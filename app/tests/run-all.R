# Run every test suite. From the app directory:  Rscript tests/run-all.R
#
# Suites, fastest first:
#   boot          parses and sources everything, builds the UI, checks that
#                 every registry entry's arguments match its function
#   test-modules  reactive logic via shiny::testServer
#   test-async    background workers agree with the foreground
#   test-registry every diagnostic runs, and the code shown reproduces it
#   test-plots    every result draws, on every preset
#   test-app      the whole app in a headless browser
#   test-rootrun  launching the way the README says, from the repository root
#
# Pass suite names to run a subset:  Rscript tests/run-all.R boot test-modules

SUITES <- c("boot", "test-modules", "test-async", "test-registry",
            "test-plots", "test-app", "test-rootrun")

args <- commandArgs(trailingOnly = TRUE)
suites <- if (length(args)) intersect(args, SUITES) else SUITES
if (!length(suites)) {
    cat("No matching suite. Available:", paste(SUITES, collapse = ", "), "\n")
    quit(status = 2)
}

Sys.setenv(NOT_CRAN = "true")
results <- list()

for (s in suites) {
    path <- file.path("tests", paste0(s, ".R"))
    cat(sprintf("\n=========== %s ===========\n", s))
    t0 <- proc.time()[["elapsed"]]
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(path))
    secs <- round(proc.time()[["elapsed"]] - t0, 1)
    results[[s]] <- list(status = status, seconds = secs)
    cat(sprintf("--- %s: %s in %.1fs ---\n", s,
                if (status == 0) "passed" else sprintf("FAILED (%d)", status), secs))
}

cat("\n============ summary ============\n")
for (s in names(results)) {
    cat(sprintf("  %-14s %-8s %6.1fs\n", s,
                if (results[[s]]$status == 0) "pass" else "FAIL",
                results[[s]]$seconds))
}
failed <- names(Filter(function(r) r$status != 0, results))
cat(sprintf("\n%d of %d suites passed, %.0fs total\n",
            length(results) - length(failed), length(results),
            sum(vapply(results, `[[`, numeric(1), "seconds"))))
if (length(failed)) quit(status = 1)
