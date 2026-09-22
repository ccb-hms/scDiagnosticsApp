# ---------------------------------------------------------------------------
# Phase 0 audit: baseline timings for the legacy app, and a parity check of the
# forked `Shiny/support/` functions against the published scDiagnostics package.
#
# Run from the repository root:
#   Rscript audit/phase0_baseline.R
#
# Writes audit/phase0-results.md and audit/phase0-results.rds.
#
# NOTE: this needs the legacy Shiny/ tree, which has since been removed. It is
# kept as the record of how phase0-results.md was produced. To re-run it,
# restore that tree first:
#
#   git checkout a52a226 -- Shiny/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(scDiagnostics)
})

set.seed(42)
options(stringsAsFactors = FALSE)

APP_DIR <- "Shiny"
OUT_MD <- "audit/phase0-results.md"
OUT_RDS <- "audit/phase0-results.rds"
dir.create("audit", showWarnings = FALSE)

results <- list()
notes <- character()

say <- function(...) cat(sprintf(...), "\n", sep = "")

# Time an expression, returning elapsed seconds and the value (or the error).
timed <- function(expr) {
    t0 <- proc.time()[["elapsed"]]
    val <- tryCatch(force(expr), error = function(e) structure(conditionMessage(e), class = "audit_error"))
    list(seconds = round(proc.time()[["elapsed"]] - t0, 2), value = val)
}

failed <- function(x) inherits(x, "audit_error")

# Peak-ish memory: size of the object graph we are holding.
mb <- function(x) round(as.numeric(object.size(x)) / 1024^2, 1)

# ---------------------------------------------------------------------------
# 1. Startup cost of the legacy app
# ---------------------------------------------------------------------------

say("== 1. Legacy startup ==")

data_files <- list.files(file.path(APP_DIR, "data"), full.names = TRUE)
disk_mb <- round(sum(file.info(data_files)$size) / 1024^2, 1)

startup <- timed({
    query_cells <- readRDS(file.path(APP_DIR, "data", "query_marrow_myeloid.rds"))
    reference_cells <- readRDS(file.path(APP_DIR, "data", "reference_marrow_myeloid.rds"))
    reference_cells_subset <- readRDS(file.path(APP_DIR, "data", "reference_subset_marrow_myeloid.rds"))
    list(query_cells, reference_cells, reference_cells_subset)
})
query_cells <- startup$value[[1]]
reference_cells <- startup$value[[2]]
reference_cells_subset <- startup$value[[3]]

results$startup <- list(
    disk_mb = disk_mb,
    load_seconds = startup$seconds,
    memory_mb = mb(startup$value),
    dims = vapply(startup$value, function(x) paste(dim(x), collapse = " x "), character(1)),
    assays = vapply(startup$value, function(x) paste(assayNames(x), collapse = ","), character(1))
)
say("  %.1f MB on disk, %.2fs to read, %.1f MB in memory",
    disk_mb, startup$seconds, results$startup$memory_mb)

# Legacy app column names on this data.
REF_COL <- "Cell Ontology Class"
QUERY_COL <- "SingleR Annotation - Promonocytes Missing in Reference"
stopifnot(REF_COL %in% colnames(colData(reference_cells)))
stopifnot(QUERY_COL %in% colnames(colData(query_cells)))

# ---------------------------------------------------------------------------
# 2. Source the forked support functions exactly as the app does
# ---------------------------------------------------------------------------

say("== 2. Sourcing forked support/ ==")

support_files <- list.files(file.path(APP_DIR, "support"), full.names = TRUE)
fork <- new.env(parent = globalenv())
for (f in support_files) sys.source(f, envir = fork)
# The forks call each other by name, so they must see one another.
parent.env(fork) <- globalenv()
for (nm in ls(fork)) {
    if (is.function(fork[[nm]])) environment(fork[[nm]]) <- fork
}
say("  %d files, %d lines, %d functions",
    length(support_files),
    sum(vapply(support_files, function(f) length(readLines(f, warn = FALSE)), integer(1))),
    length(ls(fork)))

# ---------------------------------------------------------------------------
# 3. Paired runs: fork vs package
# ---------------------------------------------------------------------------

say("== 3. Parity: fork vs package ==")

PC <- 1:5

# Recursively flatten any result object into a named vector of numeric leaves,
# so two results of unknown shape can be compared without bespoke accessors.
flatten_numeric <- function(x, prefix = "", depth = 0) {
    if (depth > 6) return(numeric(0))
    if (is.null(x)) return(numeric(0))
    if (inherits(x, c("ggplot", "gtable", "function", "igraph"))) return(numeric(0))
    if (is.data.frame(x)) {
        num <- Filter(is.numeric, as.list(x))
        if (!length(num)) return(numeric(0))
        out <- unlist(lapply(names(num), function(n) {
            v <- as.numeric(num[[n]])
            stats::setNames(v, sprintf("%s%s[%d]", paste0(prefix, "."), n, seq_along(v)))
        }))
        return(out)
    }
    if (is.matrix(x) && is.numeric(x)) {
        v <- as.numeric(x)
        return(stats::setNames(v, sprintf("%s[%d]", prefix, seq_along(v))))
    }
    if (is.numeric(x) && is.null(dim(x))) {
        nm <- if (!is.null(names(x))) paste0(prefix, ".", names(x)) else sprintf("%s[%d]", prefix, seq_along(x))
        return(stats::setNames(as.numeric(x), nm))
    }
    if (is.list(x)) {
        keys <- names(x)
        if (is.null(keys)) keys <- as.character(seq_along(x))
        out <- unlist(lapply(seq_along(x), function(i) {
            flatten_numeric(x[[i]], paste0(prefix, if (nzchar(prefix)) "." else "", keys[i]), depth + 1)
        }))
        return(if (is.null(out)) numeric(0) else out)
    }
    numeric(0)
}

# Compare two arbitrary results structurally.
compare_results <- function(a, b, tol = 1e-6) {
    fa <- tryCatch(flatten_numeric(a), error = function(e) numeric(0))
    fb <- tryCatch(flatten_numeric(b), error = function(e) numeric(0))
    if (!length(fa) && !length(fb)) {
        return(list(status = "no numeric leaves", detail = "nothing comparable (plot-only result)"))
    }
    shared <- intersect(names(fa), names(fb))
    only_fork <- setdiff(names(fa), names(fb))
    only_pkg <- setdiff(names(fb), names(fa))
    if (!length(shared)) {
        return(list(status = "INCOMPATIBLE",
                    detail = sprintf("no shared fields; fork has %d numeric leaves, package %d",
                                     length(fa), length(fb))))
    }
    d <- suppressWarnings(max(abs(fa[shared] - fb[shared]), na.rm = TRUE))
    struct <- if (length(only_fork) || length(only_pkg)) {
        sprintf("; %d fields only in fork, %d only in package", length(only_fork), length(only_pkg))
    } else ""
    if (!is.finite(d)) {
        list(status = "DIFFERENT", detail = paste0("non-finite difference", struct))
    } else if (d <= tol) {
        list(status = if (nzchar(struct)) "values agree, shape differs" else "identical",
             detail = sprintf("%d shared fields, max abs diff %.2e%s", length(shared), d, struct))
    } else {
        list(status = "DIFFERENT",
             detail = sprintf("%d shared fields, max abs diff %.3g%s", length(shared), d, struct))
    }
}

record <- function(name, fork_run, pkg_run, note = NA_character_, comparison = NULL) {
    if (is.null(comparison)) {
        comparison <- if (failed(fork_run$value) || failed(pkg_run$value)) {
            list(status = "ERROR", detail = "one side failed to run")
        } else {
            tryCatch(compare_results(fork_run$value, pkg_run$value),
                     error = function(e) list(status = "ERROR",
                                              detail = paste("comparison failed:", conditionMessage(e))))
        }
    }
    results$parity[[name]] <<- list(
        name = name,
        fork_seconds = fork_run$seconds,
        pkg_seconds = pkg_run$seconds,
        fork_error = if (failed(fork_run$value)) as.character(fork_run$value) else NA_character_,
        pkg_error = if (failed(pkg_run$value)) as.character(pkg_run$value) else NA_character_,
        status = comparison$status,
        detail = comparison$detail,
        note = note
    )
    say("  %-28s fork %6.2fs | pkg %6.2fs | %s", name,
        fork_run$seconds, pkg_run$seconds, comparison$status)
    say("  %-28s   %s", "", comparison$detail)
    if (!is.na(results$parity[[name]]$fork_error)) {
        say("  %-28s   fork error: %s", "", results$parity[[name]]$fork_error)
    }
    if (!is.na(results$parity[[name]]$pkg_error)) {
        say("  %-28s   pkg error: %s", "", results$parity[[name]]$pkg_error)
    }
}

## Each block runs the fork and the package on identical inputs under identical
## seeds, then compares every numeric leaf of the two results.

## 3a. projectPCA -----------------------------------------------------------
record("projectPCA",
       timed(fork$projectPCA(
           query_data = query_cells, reference_data = reference_cells,
           query_cell_type_col = QUERY_COL, ref_cell_type_col = REF_COL, pc_subset = PC)),
       timed(scDiagnostics::projectPCA(
           query_data = query_cells, reference_data = reference_cells,
           query_cell_type_col = QUERY_COL, ref_cell_type_col = REF_COL, pc_subset = PC)))

## 3b. detectAnomaly --------------------------------------------------------
set.seed(1)
f <- timed(fork$detectAnomaly_app(
    reference_data = reference_cells, query_data = query_cells,
    ref_cell_type_col = REF_COL, query_cell_type_col = QUERY_COL,
    pc_subset = PC, n_tree = 500, anomaly_threshold = 0.6))
set.seed(1)
p <- timed(scDiagnostics::detectAnomaly(
    reference_data = reference_cells, query_data = query_cells,
    ref_cell_type_col = REF_COL, query_cell_type_col = QUERY_COL,
    pc_subset = PC, n_tree = 500))
record("detectAnomaly", f, p,
       note = paste("fork exposes a fixed anomaly_threshold only; the package adds",
                    "threshold_method (MAD/absolute), mad_multiplier, n_hvgs and",
                    "max_cells_* downsampling"))

## 3c. calculateWassersteinDistance -----------------------------------------
set.seed(2)
f <- timed(fork$calculateWassersteinDistance_app(
    query_data = query_cells, reference_data = reference_cells,
    ref_cell_type_col = REF_COL, query_cell_type_col = QUERY_COL,
    pc_subset = PC, n_resamples = 100))
set.seed(2)
p <- timed(scDiagnostics::calculateWassersteinDistance(
    query_data = query_cells, reference_data = reference_cells,
    ref_cell_type_col = REF_COL, query_cell_type_col = QUERY_COL,
    pc_subset = PC, n_resamples = 100))
record("calculateWassersteinDistance", f, p)

## 3d. calculateDiscriminantSpace -------------------------------------------
set.seed(3)
f <- timed(fork$calculateDiscriminantSpace_app(
    reference_data = reference_cells, query_data = query_cells,
    ref_cell_type_col = REF_COL, query_cell_type_col = QUERY_COL))
set.seed(3)
p <- timed(scDiagnostics::calculateDiscriminantSpace(
    reference_data = reference_cells, query_data = query_cells,
    ref_cell_type_col = REF_COL, query_cell_type_col = QUERY_COL))
record("calculateDiscriminantSpace", f, p)

## 3e. calculateGraphIntegration --------------------------------------------
set.seed(4)
f <- timed(fork$calculateGraphIntegration_app(
    query_data = query_cells, reference_data = reference_cells,
    query_cell_type_col = QUERY_COL, ref_cell_type_col = REF_COL,
    pc_subset = 1:10, k_neighbors = 30))
set.seed(4)
p <- timed(scDiagnostics::calculateGraphIntegration(
    query_data = query_cells, reference_data = reference_cells,
    query_cell_type_col = QUERY_COL, ref_cell_type_col = REF_COL,
    pc_subset = 1:10, k_neighbors = 30))
record("calculateGraphIntegration", f, p)

## 3f. calculateVarImpOverlap -----------------------------------------------
set.seed(5)
f <- timed(fork$calculateVarImpOverlap(
    reference_data = reference_cells, query_data = query_cells,
    ref_cell_type_col = REF_COL, query_cell_type_col = QUERY_COL, n_tree = 100, n_top = 20))
set.seed(5)
p <- timed(scDiagnostics::calculateVarImpOverlap(
    reference_data = reference_cells, query_data = query_cells,
    ref_cell_type_col = REF_COL, query_cell_type_col = QUERY_COL, n_tree = 100, n_top = 20))
record("calculateVarImpOverlap", f, p,
       note = "fork uses randomForest; the package Imports ranger")

## 3g. plotCellTypePCA -------------------------------------------------------
f <- timed(fork$plotCellTypePCA_app(
    query_data = query_cells, reference_data = reference_cells,
    query_cell_type_col = QUERY_COL, ref_cell_type_col = REF_COL, pc_subset = PC))
p <- timed(scDiagnostics::plotCellTypePCA(
    query_data = query_cells, reference_data = reference_cells,
    query_cell_type_col = QUERY_COL, ref_cell_type_col = REF_COL, pc_subset = PC))
record("plotCellTypePCA", f, p,
       note = "ggplot results; the package downsamples to max_cells_* = 2000 by default, the fork does not")

## 3h. plotGeneExpressionDimred ---------------------------------------------
gene <- rownames(query_cells)[1]
f <- timed(fork$plotGeneExpressionDimred_app(
    se_object = query_cells, method = "PCA", pc_subset = 1:3,
    feature = gene, cell_type_col = QUERY_COL))
p <- timed(scDiagnostics::plotGeneExpressionDimred(
    sce_object = query_cells, method = "PCA", pc_subset = 1:3,
    feature = gene, cell_type_col = QUERY_COL))
record("plotGeneExpressionDimred", f, p,
       note = "fork argument is `se_object`, package argument is `sce_object`")


## ---------------------------------------------------------------------------
# 4. Same package functions on the package's own curated presets
# ---------------------------------------------------------------------------

say("== 4. Package functions on package preset data ==")

e <- new.env()
utils::data("reference_data", "query_data", package = "scDiagnostics", envir = e)
ref_small <- e$reference_data
qry_small <- e$query_data

preset <- list()
preset$dims <- c(ref = paste(dim(ref_small), collapse = " x "),
                 query = paste(dim(qry_small), collapse = " x "))
preset$memory_mb <- mb(list(ref_small, qry_small))

bench_small <- function(label, expr) {
    r <- timed(expr)
    preset$timings[[label]] <<- r$seconds
    say("  %-26s %6.2fs%s", label, r$seconds,
        if (failed(r$value)) paste0("  [ERROR] ", as.character(r$value)) else "")
    invisible(r)
}
set.seed(1)
bench_small("detectAnomaly", scDiagnostics::detectAnomaly(
    reference_data = ref_small, query_data = qry_small,
    ref_cell_type_col = "expert_annotation", query_cell_type_col = "SingleR_annotation",
    pc_subset = PC, n_tree = 500))
set.seed(2)
bench_small("calculateWassersteinDistance", scDiagnostics::calculateWassersteinDistance(
    query_data = qry_small, reference_data = ref_small,
    ref_cell_type_col = "expert_annotation", query_cell_type_col = "SingleR_annotation",
    pc_subset = PC, n_resamples = 100))
bench_small("plotCellTypePCA", scDiagnostics::plotCellTypePCA(
    query_data = qry_small, reference_data = ref_small,
    query_cell_type_col = "SingleR_annotation", ref_cell_type_col = "expert_annotation",
    pc_subset = PC))
set.seed(4)
bench_small("calculateGraphIntegration", scDiagnostics::calculateGraphIntegration(
    query_data = qry_small, reference_data = ref_small,
    query_cell_type_col = "SingleR_annotation", ref_cell_type_col = "expert_annotation",
    pc_subset = 1:10, k_neighbors = 30))
set.seed(6)
bench_small("calculateGeneShifts", scDiagnostics::calculateGeneShifts(
    query_data = qry_small, reference_data = ref_small,
    query_cell_type_col = "SingleR_annotation", ref_cell_type_col = "expert_annotation",
    pc_subset = 1:3))

results$preset <- preset

# ---------------------------------------------------------------------------
# 5. Report
# ---------------------------------------------------------------------------

saveRDS(results, OUT_RDS)

lines <- c(
    "# Phase 0 audit results",
    "",
    sprintf("Generated %s on R %s, scDiagnostics %s.",
            format(Sys.Date()), getRversion(), packageVersion("scDiagnostics")),
    "",
    "Timings are single runs on one machine. They are a baseline to beat, not a benchmark.",
    "",
    "## Legacy startup cost",
    "",
    sprintf("- Bundled data: **%.1f MB** on disk across %d files", disk_mb, length(data_files)),
    sprintf("- Read at process start: **%.2f s**", startup$seconds),
    sprintf("- Resident in memory: **%.1f MB**", results$startup$memory_mb),
    "",
    "| object | dimensions | assays |",
    "|---|---|---|",
    sprintf("| %s | %s | %s |",
            c("query", "reference", "reference subset"),
            results$startup$dims, results$startup$assays),
    "",
    "## Parity: forked `Shiny/support/` vs package",
    "",
    "| diagnostic | fork (s) | package (s) | agreement | detail |",
    "|---|---:|---:|---|---|",
    vapply(results$parity, function(r) {
        sprintf("| `%s` | %s | %s | %s | %s |", r$name,
                if (is.na(r$fork_error)) sprintf("%.2f", r$fork_seconds) else "error",
                if (is.na(r$pkg_error)) sprintf("%.2f", r$pkg_seconds) else "error",
                r$status, r$detail)
    }, character(1)),
    "",
    "Errors and caveats:",
    "",
    unlist(lapply(results$parity, function(r) {
        out <- character()
        if (!is.na(r$fork_error)) out <- c(out, sprintf("- `%s` fork error: %s", r$name, r$fork_error))
        if (!is.na(r$pkg_error)) out <- c(out, sprintf("- `%s` package error: %s", r$name, r$pkg_error))
        if (!is.na(r$note)) out <- c(out, sprintf("- `%s`: %s", r$name, r$note))
        out
    })),
    "",
    "## Package functions on package preset data",
    "",
    sprintf("Reference %s, query %s, %.1f MB in memory (vs %.1f MB for the legacy bundle).",
            preset$dims[["ref"]], preset$dims[["query"]], preset$memory_mb,
            results$startup$memory_mb),
    "",
    "| diagnostic | seconds |",
    "|---|---:|",
    vapply(names(preset$timings), function(k) sprintf("| `%s` | %.2f |", k, preset$timings[[k]]),
           character(1)),
    ""
)

writeLines(lines, OUT_MD)
say("\nWrote %s and %s", OUT_MD, OUT_RDS)
