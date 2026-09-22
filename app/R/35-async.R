# Background execution for the slower diagnostics.
#
# Shiny is single-threaded per R process: a twenty-second graph integration
# freezes every other session sharing that process. The slow diagnostics
# therefore run in a `future` worker behind an ExtendedTask.
#
# The worker is never sent a SingleCellExperiment. It receives a descriptor -
# a preset name or a file path - and materialises the data itself, so a large
# uploaded object is read once from disk in parallel instead of being
# serialised across the process boundary on every run. Serialising the data is
# the thing that makes naive async slower than running synchronously.

ASYNC_ENABLED <- !identical(tolower(Sys.getenv("SCDIAG_ASYNC", "true")), "false")
ASYNC_WORKERS <- suppressWarnings(as.integer(Sys.getenv("SCDIAG_WORKERS", "2")))

#' Set up the worker pool once per process
init_async <- function() {
    if (!ASYNC_ENABLED) return(invisible(FALSE))
    ok <- tryCatch({
        n <- max(1L, min(ASYNC_WORKERS %|% 2L, future::availableCores()))
        future::plan(future::multisession, workers = n)
        # Results, not inputs, cross the boundary; allow room for them.
        options(future.globals.maxSize = 2 * 1024^3)
        TRUE
    }, error = function(e) {
        warning("Could not start background workers; diagnostics will run in the main process: ",
                conditionMessage(e), call. = FALSE)
        FALSE
    })
    invisible(ok)
}

#' Everything a worker needs to reproduce a run, as plain data
#'
#' Nothing here is an S4 object, so the payload is small whatever the data.
run_spec <- function(entry, ctx, params) {
    list(
        entry_id = entry$id,
        ref_desc = ctx$ref_desc,
        query_desc = ctx$query_desc,
        ref_col = ctx$ref_col,
        query_col = ctx$query_col,
        cell_types = ctx$cell_types,
        assay = ctx$assay,
        pc_subset = ctx$pc_subset,
        focus_which = ctx$focus_which,
        preset_id = ctx$preset_id,
        params = params
    )
}

#' Run a diagnostic inside a worker process
#'
#' Sources the app once per worker, then rebuilds the context from the
#' descriptor and runs exactly the code path the main process would.
worker_run <- function(app_dir, spec) {
    if (!isTRUE(getOption("sc.worker.ready", FALSE))) {
        suppressPackageStartupMessages({
            library(scDiagnostics)
            library(SingleCellExperiment)
            library(ggplot2)
        })
        for (f in sort(list.files(file.path(app_dir, "R"), pattern = "[.][Rr]$",
                                  full.names = TRUE))) {
            source(f, local = globalenv())
        }
        options(sc.worker.ready = TRUE)
    }
    ctx <- build_context(
        ref = if (is.null(spec$ref_desc)) NULL else materialize(spec$ref_desc),
        query = if (is.null(spec$query_desc)) NULL else materialize(spec$query_desc),
        ref_desc = spec$ref_desc, query_desc = spec$query_desc,
        ref_col = spec$ref_col, query_col = spec$query_col,
        cell_types = spec$cell_types, assay = spec$assay,
        pc_subset = spec$pc_subset, focus_which = spec$focus_which %||% "query",
        preset_id = spec$preset_id
    )
    entry <- REGISTRY[[spec$entry_id]]
    args <- tryCatch(build_args(entry, ctx, spec$params),
                     error = function(e) sc_failure(conditionMessage(e)))
    if (is_failure(args)) return(args)
    run_entry(entry, args)
}

#' The app's own directory, captured at load so workers can find it
APP_DIR <- normalizePath(".", winslash = "/", mustWork = FALSE)

#' Build an ExtendedTask that runs `worker_run` in the pool
new_diagnostic_task <- function() {
    shiny::ExtendedTask$new(function(app_dir, spec) {
        promises::future_promise(
            worker_run(app_dir, spec),
            seed = TRUE,
            globals = list(app_dir = app_dir, spec = spec, worker_run = worker_run)
        )
    })
}

#' Should this run go to a worker?
#'
#' Only for entries marked heavy, only when the pool started, and only when the
#' data can be reconstructed from a descriptor. A context assembled some other
#' way falls back to the main process rather than silently disagreeing with it.
use_async <- function(entry, ctx) {
    ASYNC_ENABLED &&
        isTRUE(entry$heavy) &&
        isTRUE(getOption("sc.async.ready", FALSE)) &&
        !is.null(ctx$ref_desc) &&
        (is.null(ctx$query) || !is.null(ctx$query_desc))
}
