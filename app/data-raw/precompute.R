# Pre-warm the shared cache.
#
# Every preset result is identical for every visitor, so computing them once at
# build time and shipping the cache means the first person to open a panel
# waits no longer than the second. Only preset data is precomputed: results
# from uploads are private and never enter this cache.
#
# Run from the app directory:  Rscript data-raw/precompute.R
#
# Writes inst/precomputed/, which the app copies into its runtime cache on
# first use (see shared_cache() in R/11-cache.R).

suppressPackageStartupMessages({
    library(shiny); library(bslib); library(bsicons)
    library(scDiagnostics); library(SingleCellExperiment); library(ggplot2)
})

OUT <- "inst/precomputed"
unlink(OUT, recursive = TRUE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(SCDIAG_CACHE_DIR = normalizePath(OUT, winslash = "/", mustWork = TRUE))

for (f in sort(list.files("R", pattern = "[.][Rr]$", full.names = TRUE))) {
    source(f, local = globalenv())
}

# Candidates: the diagnostics a visitor opens first, plus the three the guided
# audit runs.
WARM <- c("plotCellTypePCA", "detectAnomaly", "calculateGeneShifts",
          "regressPC", "calculateWassersteinDistance", "calculateGraphIntegration",
          "calculateReconstructionError", "compareMarkers", "boxplotPCA")

# What is actually kept.
#
# Cache what is slow, not what is large. A ggplot object carries the data it
# was built from, so boxplotPCA weighs 17 MB and takes a fifth of a second:
# caching it would trade a lot of disk for nothing. The Wasserstein null
# distribution is the opposite - half a minute of resampling that serialises
# to 30 KB - and is exactly what a shipped cache is for.
MIN_SECONDS <- 1.0
MAX_MB <- 3.0
worth_caching <- function(secs, mb) secs >= MIN_SECONDS && mb <= MAX_MB

SEED <- 20260921
cache <- shared_cache()
rows <- list()

for (pid in names(PRESETS)) {
    p <- PRESETS[[pid]]
    rd <- preset_descriptor(p$ref$data)
    qd <- preset_descriptor(p$query$data)
    ctx <- build_context(
        ref = materialize(rd), query = materialize(qd),
        ref_desc = rd, query_desc = qd,
        ref_col = p$ref$col, query_col = p$query$col,
        pc_subset = p$pc_subset %||% 1:5, preset_id = pid
    )
    cat(sprintf("\n== %s ==\n", p$label))

    for (id in WARM) {
        entry <- REGISTRY[[id]]
        gate <- entry_ready(entry, ctx)
        if (!is.null(gate)) {
            cat(sprintf("  %-30s skipped (%s)\n", id, substr(gate, 1, 50)))
            next
        }
        pv <- default_param_values(entry$params, ctx)
        # The audit runs gene shifts over the first three components only.
        if (id == "calculateGeneShifts") pv$pc_subset <- utils::head(ctx$pc_subset, 3)

        args <- tryCatch(build_args(entry, ctx, pv), error = function(e) sc_failure(conditionMessage(e)))
        if (is_failure(args)) {
            cat(sprintf("  %-30s skipped (%s)\n", id, substr(args$message, 1, 50)))
            next
        }
        key <- cache_key(entry$id, ctx, pv)
        t0 <- proc.time()[["elapsed"]]
        set.seed(SEED)
        value <- run_entry(entry, args)
        secs <- proc.time()[["elapsed"]] - t0
        if (is_failure(value)) {
            cat(sprintf("  %-30s FAILED  %s\n", id, substr(value$message, 1, 60)))
            next
        }
        mb <- as.numeric(utils::object.size(value)) / 1024^2
        keep <- worth_caching(secs, mb)
        if (keep) cache$set(key, value)
        rows[[length(rows) + 1]] <- data.frame(preset = pid, id = id,
                                               seconds = round(secs, 2),
                                               mb = round(mb, 2), cached = keep)
        cat(sprintf("  %-30s %6.2fs  %6.2f MB  %s\n", id, secs, mb,
                    if (keep) "cached" else "not worth caching"))
    }
}

df <- do.call(rbind, rows)
kept <- df[df$cached, , drop = FALSE]
total_mb <- round(sum(file.info(list.files(OUT, recursive = TRUE, full.names = TRUE))$size) / 1024^2, 1)

cat("\n------------------------------------------------------------\n")
cat(sprintf("%d of %d results cached: %.1f MB on disk, %.0fs of compute saved per cold start\n",
            nrow(kept), nrow(df), total_mb, sum(kept$seconds)))
if (nrow(kept)) {
    cat(sprintf("ratio: %.1f seconds saved per MB shipped\n", sum(kept$seconds) / max(total_mb, 0.1)))
}
skipped <- df[!df$cached, , drop = FALSE]
if (nrow(skipped)) {
    cat(sprintf("not cached (under %.1fs or over %.0f MB): %s\n", MIN_SECONDS, MAX_MB,
                paste(unique(skipped$id), collapse = ", ")))
}

# A warmed cache that dwarfs the app is a liability, not an optimisation.
if (total_mb > 200) {
    warning(sprintf(
        "Precomputed cache is %.1f MB. Trim WARM or drop a preset before committing.",
        total_mb), call. = FALSE)
}

writeLines(c(
    "# Precomputed cache",
    "",
    sprintf("Built %s with scDiagnostics %s.", Sys.Date(), utils::packageVersion("scDiagnostics")),
    sprintf("%d results, %.1f MB, %.0f seconds of compute saved.", nrow(kept), total_mb, sum(kept$seconds)),
    "",
    "| preset | diagnostic | seconds | MB | cached |",
    "|---|---|---:|---:|---|",
    sprintf("| %s | `%s` | %.2f | %.2f | %s |", df$preset, df$id, df$seconds, df$mb,
            ifelse(df$cached, "yes", "no"))
), file.path(OUT, "MANIFEST.md"))

cat(sprintf("Wrote %s\n", OUT))
