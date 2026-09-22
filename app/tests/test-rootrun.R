# Launching the way the README says.
#
# Several things resolve relative to the working directory: the precomputed
# cache the app seeds itself from, and the app directory that background
# workers source. The README tells people to run `shiny::runApp("app")` from
# the repository root, so that has to work - not only a launch from inside the
# app directory, which is what every other suite does.
#
# Run from the app directory:  Rscript tests/test-rootrun.R

library(shinytest2)

# Starting the browser is the flakiest step. Chromote allows a fixed ten
# seconds for a command response and shinytest2 cannot raise it, so a busy
# machine can fail here before any check runs. Retry rather than report
# contention as a defect.
new_app_driver <- function(..., attempts = 3) {
    for (i in seq_len(attempts)) {
        app <- tryCatch(shinytest2::AppDriver$new(...), error = function(e) e)
        if (!inherits(app, "error")) return(app)
        if (i == attempts) {
            stop("Could not start the app in a browser after ", attempts,
                 " attempts. Last error: ", conditionMessage(app), call. = FALSE)
        }
        cat(sprintf("  (browser did not start, attempt %d of %d)
", i, attempts))
        Sys.sleep(5)
    }
}

setwd("..")   # the repository root, as a user would be

failures <- character(0)
check <- function(label, cond, detail = "") {
    if (isTRUE(cond)) cat(sprintf("  ok    %s\n", label))
    else {
        failures <<- c(failures, label)
        cat(sprintf("  FAIL  %s  %s\n", label, detail))
    }
}

cat("== starting from", getwd(), "==\n")
app <- new_app_driver(app_dir = "app", name = "scdiagnostics-rootrun",
                     load_timeout = 180 * 1000, timeout = 120 * 1000)
on.exit(try(app$stop(), silent = TRUE), add = TRUE)

# wait_for_idle can return during a lull on a loaded machine, so poll for the
# state we are actually waiting for rather than trusting a single settle.
wait_for <- function(pattern, secs = 300) {
    for (i in seq_len(secs)) {
        if (grepl(pattern, app$get_html("body"), fixed = TRUE)) return(TRUE)
        Sys.sleep(1)
    }
    FALSE
}

html <- app$get_html("body")
check("the app starts", grepl("Audit your cell type annotations", html, fixed = TRUE))
check("no package-version banner", !grepl("sc-context-bar-alert", html, fixed = TRUE),
      "the installed scDiagnostics is older than the app requires")

app$click("home-demo")
app$wait_for_idle(timeout = 300 * 1000)
check("a preset loads", wait_for("covid_reference_data"))
check("the guided audit produces a verdict", wait_for("sc-verdict"))
html <- app$get_html("body")

errs <- setdiff(unique(regmatches(html, gregexpr("shiny-output-error[^\"]*", html))[[1]]),
                "shiny-output-error-validation")
check("no errors", length(errs) == 0, paste(errs, collapse = ", "))

logs <- tryCatch(as.data.frame(app$get_logs()), error = function(e) NULL)
check("the shipped cache is found and seeded",
      !is.null(logs) && any(grepl("Seeded the result cache", logs$message)),
      "inst/precomputed was not located from this working directory")

app$stop()

cat(sprintf("\n%d failures\n", length(failures)))
if (length(failures)) {
    for (f in failures) cat("  -", f, "\n")
    quit(status = 1)
}
cat("ROOTRUN OK\n")
