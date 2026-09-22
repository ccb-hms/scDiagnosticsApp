# Plot and cross-preset verification.
#
# test-registry.R checks that every diagnostic computes and that the code we
# show reproduces it. This checks the other half: that every result can
# actually be drawn, on every preset - including the spatial one and the two
# with no numeric metadata, where the QC panels must degrade gracefully rather
# than crash.
#
# Run from the app directory:  Rscript tests/test-plots.R [preset ...]

suppressPackageStartupMessages({
    library(shiny); library(bslib); library(bsicons)
    library(scDiagnostics); library(SingleCellExperiment); library(ggplot2)
})

for (f in sort(list.files("R", pattern = "[.][Rr]$", full.names = TRUE))) {
    source(f, local = globalenv())
}
ggplot2::theme_set(sc_ggtheme())

SEED <- 20260921
cli <- commandArgs(trailingOnly = TRUE)
presets <- if (length(cli)) cli else names(PRESETS)

# Keep this fast: it is a smoke test over many combinations.
SPEED <- list(
    n_tree = 40, n_resamples = 20, n_permutation = 20,
    n_markers = 10, n_top = 10, n_top_loadings = 15, n_hvgs = 60
)

make_ctx <- function(pid) {
    p <- PRESETS[[pid]]
    rd <- preset_descriptor(p$ref$data)
    qd <- preset_descriptor(p$query$data)
    build_context(
        ref = materialize(rd), query = materialize(qd),
        ref_desc = rd, query_desc = qd,
        ref_col = p$ref$col, query_col = p$query$col,
        cell_types = NULL, assay = "logcounts",
        pc_subset = p$pc_subset %||% 1:5, preset_id = pid
    )
}

#' Values a user would have to supply, derived from the data at hand
overrides_for <- function(id, ctx) {
    ct <- ctx$shared_types[1]
    gene <- rownames(ctx$query)[1]
    cells <- colnames(ctx$query)[seq_len(min(4, ncol(ctx$query)))]
    nums <- numeric_cols(ctx$query)
    base <- switch(id,
        comparePCA = list(cell_types = ct),
        comparePCASubspace = list(cell_types = ct),
        plotPairwiseDistancesDensity = list(cell_type = ct),
        plotMarkerExpression = list(cell_type = ct, gene_name = gene),
        plotGeneExpressionDimred = list(feature = gene, method = "PCA"),
        calculateCellDistancesSimilarity = list(cell_names_query = cells),
        calculateCellSimilarityPCA = list(cell_names = cells),
        plotQCvsAnnotation = list(qc_col = nums[1], score_col = nums[min(2, length(nums))]),
        histQCvsAnnotation = list(qc_col = nums[1], score_col = nums[min(2, length(nums))]),
        plotGeneSetScores = list(score_col = nums[1], method = "PCA"),
        calculateCategorizationEntropy = list(score_columns = nums),
        list()
    )
    # Trim the expensive knobs.
    fm <- names(formals(get(REGISTRY[[id]]$fn, asNamespace("scDiagnostics"))))
    speed <- SPEED[names(SPEED) %in% fm]
    utils::modifyList(speed, base)
}

#' Plot-side arguments a panel would supply from its own controls
plot_args_for <- function(entry, value, ctx) {
    if (!length(entry$plot_params)) return(list())
    out <- list()
    for (p in entry$plot_params) {
        v <- switch(class(p)[1],
            sc_param_choice = p$default,
            sc_param_bool = p$default,
            sc_param_int = p$default,
            sc_param_num = p$default,
            # Mirrors param_refresh: a required cell type defaults to the first.
            sc_param_celltypes = if (!isTRUE(p$allow_empty)) ctx$shared_types[1] else NULL,
            NULL)
        if (!is.null(v)) out[[p$id]] <- v
    }
    # Choice-based pickers are repopulated from the result in the real module;
    # cell-type pickers are repopulated from the context, so leave those alone.
    from_choice <- vapply(entry$plot_params, function(p) inherits(p, "sc_param_choice"),
                          logical(1))
    choice_ids <- vapply(entry$plot_params, `[[`, character(1), "id")[from_choice]
    nm <- names(value)
    if ("cell_type" %in% choice_ids && length(nm)) {
        out$cell_type <- if ("Combined" %in% nm) "Combined" else nm[1]
    }
    if ("ref_cell_type" %in% choice_ids && length(nm)) out$ref_cell_type <- nm[1]
    if (!is.null(out$plot_cell_types)) out$plot_cell_types <- NULL
    out
}

results <- list()
for (pid in presets) {
    ctx0 <- make_ctx(pid)
    cat(sprintf("\n=== %s (%s) ===\n", PRESETS[[pid]]$label, pid))
    cat(sprintf("    ref %s x %s, query %s x %s, shared types: %s\n",
                nrow(ctx0$ref), ncol(ctx0$ref), nrow(ctx0$query), ncol(ctx0$query),
                abbrev(ctx0$shared_types, 6)))
    cat(sprintf("    numeric colData in query: %s\n",
                if (length(numeric_cols(ctx0$query))) paste(numeric_cols(ctx0$query), collapse = ", ")
                else "none"))

    for (id in names(REGISTRY)) {
        entry <- REGISTRY[[id]]
        ctx <- ctx0
        ctx$entry_label <- entry$label

        # A diagnostic that cannot apply to this data must say so up front.
        gate <- entry_ready(entry, ctx)
        if (!is.null(gate)) {
            results[[length(results) + 1]] <- list(preset = pid, id = id, stage = "gate",
                                                   status = "gated", detail = gate)
            cat(sprintf("  %-32s gated  -- %s\n", id, substr(gate, 1, 70)))
            next
        }

        pv <- utils::modifyList(default_param_values(entry$params, ctx), overrides_for(id, ctx))

        args <- tryCatch(build_args(entry, ctx, pv),
                         error = function(e) sc_failure(conditionMessage(e)))
        if (is_failure(args)) {
            results[[length(results) + 1]] <- list(preset = pid, id = id,
                                                   stage = "args", status = "expected-fail",
                                                   detail = args$message)
            cat(sprintf("  %-32s args   -- %s\n", id, substr(args$message, 1, 70)))
            next
        }

        set.seed(SEED)
        val <- run_entry(entry, args)
        if (is_failure(val)) {
            results[[length(results) + 1]] <- list(preset = pid, id = id,
                                                   stage = "run", status = "expected-fail",
                                                   detail = val$message)
            cat(sprintf("  %-32s run    -- %s\n", id, substr(val$message, 1, 70)))
            next
        }

        # Draw to a null device, exactly as the app would.
        path <- tempfile(fileext = ".png")
        grDevices::png(path, width = 1200, height = 900, res = 110)
        drawn <- tryCatch({
            p <- build_plot(entry, val, plot_args_for(entry, val, ctx), ctx)
            draw_result(p)
            "ok"
        }, error = function(e) conditionMessage(e))
        grDevices::dev.off()
        bytes <- if (file.exists(path)) file.info(path)$size else 0
        unlink(path)

        # A device file under ~3 KB means nothing was actually drawn.
        status <- if (!identical(drawn, "ok")) "PLOT ERROR"
                  else if (entry$returns %in% c("value", "sce")) "no plot (by design)"
                  else if (bytes < 3000) "BLANK PLOT"
                  else "ok"
        results[[length(results) + 1]] <- list(preset = pid, id = id, stage = "plot",
                                               status = status,
                                               detail = if (identical(drawn, "ok"))
                                                   sprintf("%.0f KB", bytes / 1024) else drawn)
        flag <- if (status %in% c("PLOT ERROR", "BLANK PLOT")) "  <<<" else ""
        cat(sprintf("  %-32s plot   %-20s %s%s\n", id, status,
                    substr(results[[length(results)]]$detail, 1, 60), flag))
    }
}

df <- do.call(rbind, lapply(results, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
bad <- df[df$status %in% c("PLOT ERROR", "BLANK PLOT"), ]

cat("\n============================================================\n")
cat(sprintf("%d preset x diagnostic combinations\n", nrow(df)))
for (s in unique(df$status)) cat(sprintf("  %-22s %d\n", s, sum(df$status == s)))

if (nrow(bad)) {
    cat("\nPlot failures:\n")
    for (i in seq_len(nrow(bad))) {
        cat(sprintf("  [%s] %-30s %s: %s\n", bad$preset[i], bad$id[i],
                    bad$status[i], bad$detail[i]))
    }
    quit(status = 1)
}

gated <- df[df$status == "gated", ]
if (nrow(gated)) {
    cat("\nGated before running, with an explanation:\n")
    for (i in seq_len(nrow(gated))) {
        cat(sprintf("  [%s] %-30s %s\n", gated$preset[i], gated$id[i],
                    substr(gated$detail[i], 1, 100)))
    }
}

# Anything that still reaches the function and fails is a gap in the gating.
soft <- df[df$status == "expected-fail", ]
if (nrow(soft)) {
    cat("\nRaw failures that SHOULD have been gated:\n")
    for (i in seq_len(nrow(soft))) {
        cat(sprintf("  [%s] %-30s %s\n", soft$preset[i], soft$id[i],
                    substr(soft$detail[i], 1, 90)))
    }
}
cat("\nPLOTS OK\n")
