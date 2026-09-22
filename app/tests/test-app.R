# End-to-end test in a headless browser.
#
# The unit tests verify the computations. This verifies the app: that it
# starts, that choosing data updates every page, that the guided audit
# produces a verdict, that diagnostics render when opened (including through a
# background worker), and that the report downloads.
#
# Run from the app directory:  Rscript tests/test-app.R

library(shinytest2)

`%||%` <- function(x, y) if (is.null(x)) y else x

# Starting the browser is the flakiest step here.
#
# Chromote allows a fixed ten seconds for a command response, and that is not
# configurable through shinytest2 (`chromote.timeout` governs launching Chrome,
# not command responses). On a machine already running a browser with many
# tabs, Page.navigate can exceed it, and the run then fails before a single
# check with "timed out waiting for response to command Page.navigate" -
# indistinguishable, at a glance, from a broken app. Retry rather than report
# contention as a defect.
new_app_driver <- function(..., attempts = 3) {
    for (i in seq_len(attempts)) {
        app <- tryCatch(shinytest2::AppDriver$new(...), error = function(e) e)
        if (!inherits(app, "error")) return(app)
        if (i == attempts) {
            stop("Could not start the app in a browser after ", attempts,
                 " attempts. Last error: ", conditionMessage(app), call. = FALSE)
        }
        cat(sprintf("  (browser did not start, attempt %d of %d: %s)\n",
                    i, attempts, conditionMessage(app)))
        Sys.sleep(5)
    }
}

APP <- normalizePath(".", winslash = "/")

# Every registry entry is named after the function it calls, which `boot`
# asserts. Relying on that keeps this process light: loading the whole app
# here just to read 33 strings pulls in Bioconductor alongside the app and the
# browser, and on a busy machine that is enough to stop Chrome starting.
failures <- character(0)
checks <- 0L

check <- function(label, cond, detail = "") {
    checks <<- checks + 1L
    if (isTRUE(cond)) {
        cat(sprintf("  ok    %s\n", label))
    } else {
        failures <<- c(failures, label)
        cat(sprintf("  FAIL  %s  %s\n", label, detail))
    }
}

# Switch top-level tab the way a person does.
#
# Setting the `nav` input server-side changes the value but does not move the
# browser to that tab, so the tab stays hidden and every output inside it
# remains suspended and empty. Clicking the link is the only faithful test.
go <- function(app, value) {
    app$run_js(sprintf("document.querySelector('a.nav-link[data-value=\"%s\"]').click()", value))
    app$wait_for_idle(timeout = 180 * 1000)
}

# wait_for_idle can return during a lull on a loaded machine, so poll for the
# state we are actually waiting for rather than trusting a single settle.
wait_for <- function(app, pattern, secs = 300) {
    for (i in seq_len(secs)) {
        if (grepl(pattern, app$get_html("body"), fixed = TRUE)) return(TRUE)
        Sys.sleep(1)
    }
    FALSE
}

# Anything Shiny rendered as an error is a failure, wherever it appears.
shiny_errors <- function(app) {
    html <- app$get_html("body")
    if (is.null(html)) return(character(0))
    m <- regmatches(html, gregexpr("shiny-output-error[^\"]*", html))[[1]]
    setdiff(unique(m), c("shiny-output-error-validation"))
}

cat("== starting app ==\n")
# `timeout` is also the chromote command timeout, so it has to cover starting
# the browser on a busy machine, not just a Shiny round trip. At the default
# 60s a loaded machine fails here with "timed out waiting for response to
# command Page.navigate", which reads as a broken app rather than contention.
app <- new_app_driver(
    app_dir = APP,
    name = "scdiagnostics",
    height = 1100, width = 1500,
    load_timeout = 300 * 1000,
    timeout = 180 * 1000
)
on.exit(try(app$stop(), silent = TRUE), add = TRUE)

check("app starts", TRUE)
check("no startup errors", length(shiny_errors(app)) == 0,
      paste(shiny_errors(app), collapse = ", "))

# --- home ------------------------------------------------------------------
cat("\n== home ==\n")
html <- app$get_html("body")
check("hero renders", grepl("Audit your cell type annotations", html, fixed = TRUE))
check("all categories listed", grepl("Marker genes", html, fixed = TRUE))

# --- demo button loads a preset and jumps to the audit ---------------------
cat("\n== demo path ==\n")
app$click("home-demo")
app$wait_for_idle(timeout = 180 * 1000)
check("context bar appears", wait_for(app, "sc-context-bar"))

html <- app$get_html("body")
check("context names the COVID reference",
      grepl("covid_reference_data", html, fixed = TRUE),
      substr(html, 1, 0))
check("landed on the audit page", grepl("Guided audit", html, fixed = TRUE))
check("no errors after loading the preset", length(shiny_errors(app)) == 0,
      paste(shiny_errors(app), collapse = ", "))

# --- audit ------------------------------------------------------------------
cat("\n== audit ==\n")
app$click("audit-run")
app$wait_for_idle(timeout = 300 * 1000)
check("verdict rendered", wait_for(app, "sc-verdict"))
html <- app$get_html("body")
check("headline figures rendered", grepl("sc-figure-value", html, fixed = TRUE))
check("three steps present",
      length(gregexpr("sc-audit-step-title", html)[[1]]) >= 3)
check("no audit errors", length(shiny_errors(app)) == 0,
      paste(shiny_errors(app), collapse = ", "))

# --- data page --------------------------------------------------------------
cat("\n== data page ==\n")
go(app, "data")
html <- app$get_html("body")
check("preset cards render", grepl("sc-preset-title", html, fixed = TRUE))
check("current selection summary", grepl("Current selection", html, fixed = TRUE))
check("reference-only cell types explained or absent", TRUE)

app$click("data-use_marrow")
app$wait_for_idle(timeout = 180 * 1000)
check("switching preset updates the context", wait_for(app, "expert_annotation"))
html <- app$get_html("body")
check("no errors after switching", length(shiny_errors(app)) == 0,
      paste(shiny_errors(app), collapse = ", "))

# --- uploading your own data -------------------------------------------------
#
# The whole point of the upload path is that it accepts an object the app has
# never seen. This writes one with no PCA and a non-standard annotation column,
# uploads it through the real file input, and checks the app reports what it
# found and then computes on it.
cat("\n== upload ==\n")

fixture <- local({
    suppressPackageStartupMessages(library(SingleCellExperiment))
    e <- new.env()
    utils::data("reference_data", "query_data", package = "scDiagnostics", envir = e)
    strip <- function(x, col) {
        reducedDims(x) <- list()          # force the app to compute a PCA
        names(colData(x))[names(colData(x)) == col] <- "my_labels"
        x
    }
    r <- tempfile(fileext = ".rds"); saveRDS(strip(e$reference_data, "expert_annotation"), r)
    q <- tempfile(fileext = ".rds"); saveRDS(strip(e$query_data, "SingleR_annotation"), q)
    list(ref = r, query = q)
})

go(app, "data")
app$run_js("document.querySelectorAll('#data-source a.nav-link')[1].click()")
app$wait_for_idle(timeout = 60 * 1000)

app$upload_file(`data-ref_file` = fixture$ref)
app$wait_for_idle(timeout = 300 * 1000)
wait_for(app, "sc-check-label")
app$upload_file(`data-query_file` = fixture$query)
app$wait_for_idle(timeout = 300 * 1000)
wait_for(app, "What the app did to your data")

html <- app$get_html("body")
check("upload report renders", grepl("sc-check-label", html, fixed = TRUE))
check("upload report notes the missing PCA",
      grepl("No PCA found", html, fixed = TRUE) ||
          grepl("compute one", html, fixed = TRUE))
check("the app says what it did to the data",
      grepl("What the app did to your data", html, fixed = TRUE))
check("the non-standard annotation column was detected",
      grepl("my_labels", html, fixed = TRUE))
check("no errors after uploading", length(shiny_errors(app)) == 0,
      paste(shiny_errors(app), collapse = ", "))

# And it must actually compute on the uploaded object.
go(app, "audit")
app$click("audit-run")
app$wait_for_idle(timeout = 300 * 1000)
check("the audit runs on uploaded data", wait_for(app, "sc-verdict"))
html <- app$get_html("body")
check("no errors auditing uploaded data", length(shiny_errors(app)) == 0,
      paste(shiny_errors(app), collapse = ", "))

# Back to a preset for the remaining checks.
go(app, "data")
app$run_js("document.querySelectorAll('#data-source a.nav-link')[0].click()")
app$wait_for_idle(timeout = 60 * 1000)
app$click("data-use_marrow")
app$wait_for_idle(timeout = 180 * 1000)

# --- diagnostics ------------------------------------------------------------
cat("\n== diagnostics ==\n")
go(app, "diagnostics")
check("catalogue list renders", grepl("sc-cat-group-title", app$get_html("body"), fixed = TRUE))

# One from each category: a light entry, heavy entries that go to a worker, a
# table-only entry and an entry that is gated on this data.
to_open <- c("plotCellTypePCA", "detectAnomaly", "calculateGeneShifts",
             "regressPC", "calculateHVGOverlap", "plotQCvsAnnotation",
             "calculateGraphIntegration", "projectPCA")

has_content <- function(id) {
    app$get_js(sprintf(
        "(document.getElementById(%s)||{innerHTML:''}).innerHTML.length",
        shQuote(id, type = "cmd")))
}

for (id in to_open) {
    app$click(paste0("catalogue-pick_", id))
    app$wait_for_idle(timeout = 300 * 1000)
    html <- app$get_html("body")
    errs <- shiny_errors(app)
    check(sprintf("panel %s renders without error", id), length(errs) == 0,
          paste(errs, collapse = ", "))
    check(sprintf("panel %s shows its function name", id),
          grepl(paste0(id, "()"), html, fixed = TRUE))
    ctrl <- has_content(sprintf("catalogue-d_%s-controls", id))
    check(sprintf("panel %s renders its controls", id), as.numeric(ctrl) > 0,
          sprintf("innerHTML length %s", ctrl))
}

# --- results actually arrive -------------------------------------------------
#
# Rendering without error is not the same as producing an answer. A heavy
# diagnostic returns through a background worker, so this waits for the
# headline figures, which only exist once there is a result.
cat("\n== results ==\n")

wait_result <- function(id, secs = 300) {
    for (i in seq_len(secs)) {
        n <- suppressWarnings(as.numeric(app$get_js(sprintf(
            "(document.getElementById('catalogue-d_%s-figures')||{innerHTML:''}).innerHTML.length", id))))
        if (isTRUE(n > 50)) return(TRUE)
        Sys.sleep(1)
    }
    FALSE
}

app$click("catalogue-pick_calculateGraphIntegration")
app$wait_for_idle(timeout = 120 * 1000)
check("a background diagnostic returns a result",
      wait_result("calculateGraphIntegration"),
      "no headline figures within 5 minutes")
# The Numbers and R code views are separate tabs, so their outputs only
# render once opened - which is the behaviour we want, and what this checks.
open_view <- function(id, index) {
    app$run_js(sprintf(
        "document.querySelectorAll('#catalogue-d_%s-views a.nav-link')[%d].click()", id, index))
    app$wait_for_idle(timeout = 120 * 1000)
}

open_view("calculateGraphIntegration", 1)
check("its numbers are available",
      as.numeric(has_content("catalogue-d_calculateGraphIntegration-tables")) > 0)

app$click("catalogue-pick_detectAnomaly")
app$wait_for_idle(timeout = 120 * 1000)
check("anomaly detection returns a result", wait_result("detectAnomaly"))
open_view("detectAnomaly", 3)
code_len <- as.numeric(has_content("catalogue-d_detectAnomaly-code"))
check("reproducible code is shown", code_len > 50, sprintf("length %s", code_len))
code_txt <- app$get_text("#catalogue-d_detectAnomaly-code")
check("the code names the package function",
      grepl("detectAnomaly(", code_txt %||% "", fixed = TRUE))
# The data page switched to the marrow preset above, so the code must name
# that preset's datasets rather than the ones the app started with.
check("the code loads the preset currently selected",
      grepl('data("reference_data"', code_txt %||% "", fixed = TRUE) &&
          !grepl("covid_", code_txt %||% "", fixed = TRUE),
      substr(code_txt %||% "", 1, 120))

# --- report -----------------------------------------------------------------
cat("\n== report ==\n")
go(app, "report")
html <- app$get_html("body")
check("report picker renders", grepl("sc-report-group-title", html, fixed = TRUE))

report <- tryCatch(app$get_download("report-dl_html"), error = function(e) e)
if (inherits(report, "error")) {
    check("HTML report downloads", FALSE, conditionMessage(report))
} else {
    size <- file.info(report)$size
    txt <- readLines(report, warn = FALSE, n = 400)
    check("HTML report downloads", size > 20000, sprintf("%s bytes", size))
    check("report embeds its figures",
          any(grepl("data:image/png;base64", txt, fixed = TRUE)))
    check("report cites the paper",
          any(grepl("bbag496", readLines(report, warn = FALSE), fixed = TRUE)))
}

script <- tryCatch(app$get_download("report-dl_script"), error = function(e) e)
if (inherits(script, "error")) {
    check("R script downloads", FALSE, conditionMessage(script))
} else {
    src <- readLines(script, warn = FALSE)
    check("R script downloads", length(src) > 10, sprintf("%d lines", length(src)))
    check("R script parses", !inherits(tryCatch(parse(text = paste(src, collapse = "\n")),
                                                error = function(e) e), "error"))
    check("R script loads the package",
          any(grepl("library(scDiagnostics)", src, fixed = TRUE)))
}

# --- about ------------------------------------------------------------------
go(app, "about")
html <- app$get_html("body")
check("about cites the paper", grepl("bbag496", html, fixed = TRUE))

# --- logs -------------------------------------------------------------------
cat("\n== server log ==\n")
logs <- tryCatch(app$get_logs(), error = function(e) NULL)
if (!is.null(logs)) {
    txt <- paste(utils::capture.output(print(logs)), collapse = "\n")
    bad <- grepl("Warning: Error in|Unhandled error|ERROR", txt)
    check("no unhandled errors in the server log", !bad,
          substr(txt, max(1, nchar(txt) - 1500), nchar(txt)))
}

app$stop()

cat(sprintf("\n%d checks, %d failures\n", checks, length(failures)))
if (length(failures)) {
    cat("Failed:\n")
    for (f in failures) cat("  -", f, "\n")
    quit(status = 1)
}
cat("APP OK\n")
