# Registry verification.
#
# For every diagnostic, on real data:
#   1. does it run through the app's own code path?
#   2. does the R code the app shows the user, evaluated in a clean
#      environment, reproduce the same numbers?
#
# Check 2 is the one that matters. A reproducible-code panel that does not
# actually reproduce the panel is worse than no panel at all, and only an
# end-to-end evaluation catches a mis-deparsed argument.
#
# Run from the app directory:  Rscript tests/test-registry.R [entry_id ...]

suppressPackageStartupMessages({
    library(shiny); library(bslib); library(bsicons)
    library(scDiagnostics); library(SingleCellExperiment); library(ggplot2)
})

for (f in sort(list.files("R", pattern = "[.][Rr]$", full.names = TRUE))) {
    source(f, local = globalenv())
}

SEED <- 20260921
args_cli <- commandArgs(trailingOnly = TRUE)

# --- fixture ---------------------------------------------------------------

preset <- PRESETS$marrow
ref_desc <- preset_descriptor(preset$ref$data)
query_desc <- preset_descriptor(preset$query$data)
ctx <- build_context(
    ref = materialize(ref_desc), query = materialize(query_desc),
    ref_desc = ref_desc, query_desc = query_desc,
    ref_col = preset$ref$col, query_col = preset$query$col,
    cell_types = NULL, assay = "logcounts", pc_subset = 1:5,
    preset_id = "marrow"
)

FOCUS_TYPE <- ctx$shared_types[1]
GENE <- rownames(ctx$query)[1]
CELLS <- colnames(ctx$query)[1:4]
REF_CELLS <- colnames(ctx$ref)[1:4]

# Values for parameters that have no usable default: a gene, some cells, a
# single cell type. These stand in for what a user would pick.
OVERRIDES <- list(
    comparePCA = list(cell_types = FOCUS_TYPE),
    comparePCASubspace = list(cell_types = FOCUS_TYPE),
    plotPairwiseDistancesDensity = list(cell_type = FOCUS_TYPE),
    plotMarkerExpression = list(cell_type = FOCUS_TYPE, gene_name = GENE),
    plotGeneExpressionDimred = list(feature = GENE, method = "PCA"),
    calculateCellDistancesSimilarity = list(cell_names_query = CELLS),
    # The focus object is the query, so the cells must come from the query.
    calculateCellSimilarityPCA = list(cell_names = CELLS),
    plotQCvsAnnotation = list(qc_col = "percent_mito", score_col = "annotation_scores"),
    histQCvsAnnotation = list(qc_col = "percent_mito", score_col = "annotation_scores"),
    plotGeneSetScores = list(score_col = "gene_set_scores", method = "PCA"),
    calculateCategorizationEntropy = list(
        score_columns = c("percent_mito", "gene_set_scores", "annotation_scores")
    ),
    # Keep permutation counts low: this is a correctness test, not a benchmark.
    calculateHotellingPValue = list(n_permutation = 30),
    calculateMMDPValue = list(n_permutation = 20),
    calculateWassersteinDistance = list(n_resamples = 30),
    detectAnomaly = list(n_tree = 60),
    calculateDiscriminantSpace = list(n_tree = 60),
    calculateVarImpOverlap = list(n_tree = 60, n_top = 15),
    calculateGeneShifts = list(pc_subset = 1:2, n_top_loadings = 20),
    compareMarkers = list(n_markers = 15)
)

# calculateCategorizationEntropy needs one score per candidate cell type. The
# marrow preset carries unrelated numeric columns, so entropy over them is
# meaningless as biology but exercises the code path, which is the point here.

entry_ctx <- function(entry) {
    c <- ctx
    if (uses_focus(entry)) {
        c <- build_context(
            ref = ctx$ref, query = ctx$query, ref_desc = ctx$ref_desc,
            query_desc = ctx$query_desc, ref_col = ctx$ref_col,
            query_col = ctx$query_col, cell_types = NULL, assay = ctx$assay,
            pc_subset = ctx$pc_subset, focus_which = "query", preset_id = "marrow"
        )
    }
    c$entry_label <- entry$label
    c
}

# --- comparison -------------------------------------------------------------

flatten_num <- function(x, prefix = "", depth = 0) {
    if (depth > 6 || is.null(x)) return(numeric(0))
    if (inherits(x, c("ggplot", "gg", "ggmatrix", "igraph", "function",
                      "Heatmap", "HeatmapList"))) return(numeric(0))
    if (isS4(x)) return(numeric(0))
    if (is.data.frame(x)) {
        num <- Filter(is.numeric, as.list(x))
        if (!length(num)) return(numeric(0))
        return(unlist(lapply(names(num), function(n) {
            v <- as.numeric(num[[n]])
            stats::setNames(v, sprintf("%s.%s[%d]", prefix, n, seq_along(v)))
        })))
    }
    if (is.matrix(x) && is.numeric(x)) {
        return(stats::setNames(as.numeric(x), sprintf("%s[%d]", prefix, seq_len(length(x)))))
    }
    if (is.numeric(x) && is.null(dim(x))) {
        nm <- if (!is.null(names(x))) paste0(prefix, ".", names(x)) else sprintf("%s[%d]", prefix, seq_along(x))
        return(stats::setNames(as.numeric(x), nm))
    }
    if (is.list(x)) {
        keys <- names(x) %||% as.character(seq_along(x))
        out <- unlist(lapply(seq_along(x), function(i) {
            flatten_num(x[[i]], paste0(prefix, if (nzchar(prefix)) "." else "", keys[i]), depth + 1)
        }))
        return(if (is.null(out)) numeric(0) else out)
    }
    numeric(0)
}

compare <- function(a, b, tol = 1e-8) {
    fa <- flatten_num(a); fb <- flatten_num(b)
    if (!length(fa) && !length(fb)) {
        # Both plot-only: compare the rendered data instead where possible.
        if (inherits(a, "gg") && inherits(b, "gg")) {
            da <- tryCatch(ggplot2::ggplot_build(a)$data, error = function(e) NULL)
            db <- tryCatch(ggplot2::ggplot_build(b)$data, error = function(e) NULL)
            if (!is.null(da) && !is.null(db)) {
                fa <- flatten_num(da); fb <- flatten_num(db)
            }
        }
        if (!length(fa) && !length(fb)) return(list(ok = TRUE, detail = "no numeric output to compare"))
    }
    shared <- intersect(names(fa), names(fb))
    if (!length(shared)) {
        return(list(ok = FALSE, detail = sprintf("no shared fields (%d vs %d leaves)",
                                                 length(fa), length(fb))))
    }
    if (length(shared) != length(fa) || length(shared) != length(fb)) {
        return(list(ok = FALSE, detail = sprintf(
            "structures differ: %d shared of %d / %d", length(shared), length(fa), length(fb))))
    }
    d <- suppressWarnings(max(abs(fa[shared] - fb[shared]), na.rm = TRUE))
    if (!is.finite(d)) return(list(ok = TRUE, detail = "all-NA comparison"))
    list(ok = d <= tol, detail = sprintf("%d values, max abs diff %.2e", length(shared), d))
}

# --- run --------------------------------------------------------------------

ids <- if (length(args_cli)) args_cli else names(REGISTRY)
rows <- list()

cat(sprintf("Verifying %d diagnostics on the %s preset\n", length(ids), preset$label))
cat(sprintf("Reference %s x %s, query %s x %s\n\n",
            nrow(ctx$ref), ncol(ctx$ref), nrow(ctx$query), ncol(ctx$query)))

for (id in ids) {
    entry <- REGISTRY[[id]]
    c <- entry_ctx(entry)
    pv <- default_param_values(entry$params, c)
    pv <- utils::modifyList(pv, OVERRIDES[[id]] %||% list())

    t0 <- proc.time()[["elapsed"]]
    args <- tryCatch(build_args(entry, c, pv), error = function(e) sc_failure(conditionMessage(e)))
    if (is_failure(args)) {
        rows[[id]] <- list(id = id, run = "ERROR", repro = "-", secs = 0,
                           detail = paste("build_args:", args$message))
        cat(sprintf("%-32s ERROR  %s\n", id, args$message))
        next
    }

    set.seed(SEED)
    got <- run_entry(entry, args)
    secs <- round(proc.time()[["elapsed"]] - t0, 2)

    if (is_failure(got)) {
        rows[[id]] <- list(id = id, run = "ERROR", repro = "-", secs = secs,
                           detail = got$message)
        cat(sprintf("%-32s ERROR  (%5.2fs) %s\n", id, secs, got$message))
        next
    }

    # Now the real check: evaluate the code the app would show the user.
    code <- repro_code(entry, c, args, NULL, pv)
    # Strip the plot() call: we compare the computed object, not the device.
    code <- sub("\n\nplot\\(.*$", "", code)
    env <- new.env(parent = globalenv())
    set.seed(SEED)
    repro <- tryCatch({
        eval(parse(text = code), envir = env)
        get("result", envir = env)
    }, error = function(e) sc_failure(conditionMessage(e)))

    if (is_failure(repro)) {
        rows[[id]] <- list(id = id, run = "ok", repro = "ERROR", secs = secs,
                           detail = repro$message)
        cat(sprintf("%-32s ok     (%5.2fs)  repro FAILED: %s\n", id, secs, repro$message))
        next
    }

    cmp <- compare(got, repro)
    rows[[id]] <- list(id = id, run = "ok", repro = if (cmp$ok) "match" else "DIFFERS",
                       secs = secs, detail = cmp$detail)
    cat(sprintf("%-32s ok     (%5.2fs)  repro %-8s %s\n", id, secs,
                if (cmp$ok) "match" else "DIFFERS", cmp$detail))
}

# --- summary -----------------------------------------------------------------

df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
n_err <- sum(df$run == "ERROR")
n_repro_bad <- sum(df$repro %in% c("ERROR", "DIFFERS"))

cat("\n------------------------------------------------------------\n")
cat(sprintf("%d diagnostics: %d ran, %d failed to run\n",
            nrow(df), sum(df$run == "ok"), n_err))
cat(sprintf("reproducible code: %d match, %d differ or error\n",
            sum(df$repro == "match"), n_repro_bad))
cat(sprintf("total compute: %.1fs\n", sum(df$secs)))

if (n_err || n_repro_bad) {
    cat("\nProblems:\n")
    bad <- df[df$run == "ERROR" | df$repro %in% c("ERROR", "DIFFERS"), ]
    for (i in seq_len(nrow(bad))) {
        cat(sprintf("  %-32s run=%-6s repro=%-8s %s\n",
                    bad$id[i], bad$run[i], bad$repro[i], bad$detail[i]))
    }
    quit(status = 1)
}
cat("\nALL GOOD\n")
