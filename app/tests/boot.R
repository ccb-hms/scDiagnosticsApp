# Boot check: parse every file, load the app environment, and build the UI.
# Catches syntax errors, missing functions and malformed tags without needing
# a browser. Run from the app directory:  Rscript tests/boot.R

suppressPackageStartupMessages({
    library(shiny)
    library(bslib)
    library(bsicons)
    library(scDiagnostics)
    library(SingleCellExperiment)
    library(ggplot2)
})

fail <- function(...) {
    cat("FAIL:", sprintf(...), "\n")
    quit(status = 1)
}
ok <- function(...) cat("  ok  ", sprintf(...), "\n")

files <- sort(list.files("R", pattern = "[.][Rr]$", full.names = TRUE))
cat("== parsing", length(files), "files ==\n")
for (f in files) {
    e <- tryCatch({
        parse(f)
        NULL
    }, error = function(e) e)
    if (!is.null(e)) fail("%s does not parse: %s", f, conditionMessage(e))
}
ok("all files parse")

cat("== sourcing ==\n")
for (f in files) {
    e <- tryCatch({
        source(f, local = globalenv())
        NULL
    }, error = function(e) e)
    if (!is.null(e)) fail("%s failed to source: %s", f, conditionMessage(e))
}
ok("all files source")

cat("== registry ==\n")
if (!exists("REGISTRY")) fail("REGISTRY not defined")
ok("%d entries", length(REGISTRY))

exported <- setdiff(
    Filter(function(f) is.function(tryCatch(get(f, asNamespace("scDiagnostics")),
                                            error = function(e) NULL)),
           getNamespaceExports("scDiagnostics")),
    grep("^plot[.]", getNamespaceExports("scDiagnostics"), value = TRUE)
)
covered <- vapply(REGISTRY, `[[`, character(1), "fn")
missing <- setdiff(exported, covered)
extra <- setdiff(covered, exported)
if (length(missing)) fail("not exposed in the app: %s", paste(missing, collapse = ", "))
if (length(extra)) fail("registry names a function the package does not export: %s",
                        paste(extra, collapse = ", "))
ok("all %d exported functions are covered", length(exported))

cat("== entry validity ==\n")
for (e in REGISTRY) {
    fn <- tryCatch(get(e$fn, asNamespace("scDiagnostics")), error = function(x) NULL)
    if (is.null(fn)) fail("%s: no such function %s", e$id, e$fn)
    fm <- names(formals(fn))

    bad <- setdiff(names(e$data_map), fm)
    if (length(bad)) fail("%s: data_map names arguments %s does not have: %s",
                          e$id, e$fn, paste(bad, collapse = ", "))

    # Parameters must map onto real arguments, unless a prepare() step
    # consumes them or they drive the plot rather than the call.
    pids <- vapply(e$params, `[[`, character(1), "id")
    unknown <- setdiff(pids, fm)
    if (length(unknown) && is.null(e$prepare)) {
        fail("%s: parameters %s are not arguments of %s and there is no prepare()",
             e$id, paste(unknown, collapse = ", "), e$fn)
    }

    if (length(e$plot_params)) {
        cls <- switch(e$returns, object = paste0(e$fn, "Object"), NA_character_)
        if (!is.na(cls)) {
            m <- utils::getS3method("plot", cls, optional = TRUE)
            if (is.null(m)) fail("%s: declares plot_params but no plot.%s method exists", e$id, cls)
            pm <- names(formals(m))
            if (!"..." %in% pm) {
                ppids <- setdiff(vapply(e$plot_params, `[[`, character(1), "id"),
                                 c("plot_cell_types"))
                bad2 <- setdiff(ppids, pm)
                if (length(bad2)) fail("%s: plot_params %s not accepted by plot.%s",
                                       e$id, paste(bad2, collapse = ", "), cls)
            }
        }
    }

    if (is.null(e$reading)) fail("%s: no reading guidance", e$id)
    if (!nzchar(e$blurb)) fail("%s: no blurb", e$id)
}
ok("every entry's arguments match its function")

for (e in REGISTRY) {
    # The browser suite derives the function name from the entry id rather
    # than loading the whole app into its own process, so the two must agree.
    if (!identical(e$id, e$fn)) {
        fail("entry id %s does not match its function name %s", e$id, e$fn)
    }
}
ok("every entry id matches its function name")



cat("== UI construction ==\n")
ui <- tryCatch({
    src <- readLines("app.R")
    # Evaluate app.R up to (but not including) the shinyApp call.
    stop_at <- grep("^shiny::shinyApp", src)
    eval(parse(text = paste(src[seq_len(stop_at[1] - 1)], collapse = "\n")), globalenv())
    get("ui", globalenv())
}, error = function(e) e)
if (inherits(ui, "error")) fail("app.R failed: %s", conditionMessage(ui))

html <- tryCatch(as.character(ui), error = function(e) e)
if (inherits(html, "error")) fail("UI does not render: %s", conditionMessage(html))
ok("UI renders, %s KB of HTML", format(round(nchar(html) / 1024, 1)))

for (id in names(REGISTRY)) {
    h <- tryCatch(as.character(diagnostic_ui(paste0("diag_", id), REGISTRY[[id]])),
                  error = function(e) e)
    if (inherits(h, "error")) fail("%s UI failed: %s", id, conditionMessage(h))
}
ok("all %d diagnostic panels render", length(REGISTRY))

cat("\nBOOT OK\n")
