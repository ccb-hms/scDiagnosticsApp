# Background execution.
#
# Checks the thing async is for: that a heavy diagnostic really runs in a
# separate process, that it produces the same answer as running it in the
# foreground, and that the descriptor mechanism keeps the data out of the
# message sent to the worker.
#
# Run from the app directory:  Rscript tests/test-async.R

suppressPackageStartupMessages({
    library(shiny); library(bslib); library(bsicons)
    library(scDiagnostics); library(SingleCellExperiment); library(ggplot2)
    library(promises); library(future)
})

for (f in sort(list.files("R", pattern = "[.][Rr]$", full.names = TRUE))) {
    source(f, local = globalenv())
}

failures <- character(0)
check <- function(label, cond, detail = "") {
    if (isTRUE(cond)) cat(sprintf("  ok    %s\n", label))
    else {
        failures <<- c(failures, label)
        cat(sprintf("  FAIL  %s  %s\n", label, detail))
    }
}

cat("== worker pool ==\n")
ok <- init_async()
check("workers start", isTRUE(ok))
if (!isTRUE(ok)) {
    cat("\nSkipping: no worker pool available on this machine.\n")
    quit(status = 0)
}
options(sc.async.ready = TRUE)

ctx <- local({
    p <- PRESETS$marrow
    rd <- preset_descriptor(p$ref$data); qd <- preset_descriptor(p$query$data)
    build_context(ref = materialize(rd), query = materialize(qd),
                  ref_desc = rd, query_desc = qd,
                  ref_col = p$ref$col, query_col = p$query$col,
                  pc_subset = 1:5, preset_id = "marrow")
})

entry <- REGISTRY$detectAnomaly
params <- utils::modifyList(default_param_values(entry$params, ctx), list(n_tree = 80L))

cat("\n== payload ==\n")
spec <- run_spec(entry, ctx, params)
payload_kb <- as.numeric(object.size(spec)) / 1024
check("the worker payload carries no SingleCellExperiment",
      !any(vapply(spec, function(x) inherits(x, "SummarizedExperiment"), logical(1))))
# The data itself is megabytes; the descriptor is a few hundred bytes.
data_kb <- as.numeric(object.size(list(ctx$ref, ctx$query))) / 1024
check("payload is far smaller than the data", payload_kb < data_kb / 100,
      sprintf("%.1f KB payload vs %.0f KB of data", payload_kb, data_kb))

cat("\n== same answer in a worker as in the foreground ==\n")
check("async is selected for a heavy entry", use_async(entry, ctx))

args <- build_args(entry, ctx, params)
set.seed(20260921)
here <- run_entry(entry, args)

f <- future::future(
    {
        set.seed(20260921)
        worker_run(app_dir, spec)
    },
    seed = TRUE,
    globals = list(app_dir = APP_DIR, spec = spec, worker_run = worker_run)
)
there <- future::value(f)

check("worker returns a result", !is_failure(there),
      if (is_failure(there)) there$message else "")

if (!is_failure(there) && !is_failure(here)) {
    types <- intersect(names(here), names(there))
    a <- unlist(lapply(types, function(k) here[[k]][["query_anomaly_scores"]]))
    b <- unlist(lapply(types, function(k) there[[k]][["query_anomaly_scores"]]))
    check("worker result matches the foreground result",
          length(a) == length(b) && max(abs(a - b)) < 1e-10,
          sprintf("%d vs %d values", length(a), length(b)))
}

cat("\n== worker runs in another process ==\n")
pid_here <- Sys.getpid()
pid_there <- future::value(future::future(Sys.getpid(), seed = TRUE))
check("worker has its own process id", !identical(pid_here, pid_there),
      sprintf("both %s", pid_here))

cat("\n== fallback ==\n")
upload_ctx <- ctx
upload_ctx$ref_desc <- NULL
check("a context with no descriptor falls back to the foreground",
      !use_async(entry, upload_ctx))
check("light entries never go to a worker",
      !use_async(REGISTRY$projectPCA, ctx))

future::plan(future::sequential)

cat(sprintf("\n%d failures\n", length(failures)))
if (length(failures)) {
    for (x in failures) cat("  -", x, "\n")
    quit(status = 1)
}
cat("ASYNC OK\n")
