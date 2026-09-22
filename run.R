# Launch the app locally, from the repository root:
#
#   Rscript run.R
#
# A script rather than `Rscript -e '...'` because quoting an inline expression
# differs between shells. In PowerShell the inner double quotes of
# `-e 'shiny::runApp("app")'` are stripped, so R receives a bare symbol and
# reports the unhelpful "object 'app' not found". There is nothing to quote
# here, so there is nothing to get wrong.
#
# It also checks the environment before starting. A machine with more than one
# R installation can easily run this under the wrong one, where the library is
# a different set of packages; the failure then surfaces as a missing package
# deep in startup rather than as the version problem it actually is.

if (!file.exists(file.path("app", "app.R"))) {
    stop("Cannot find app/app.R. Run this from the repository root:  Rscript run.R",
         call. = FALSE)
}

MIN_SCDIAG <- "1.7.6"

REQUIRED <- c(
    "shiny", "bslib", "bsicons", "htmltools", "DT", "cachem", "digest",
    "future", "promises", "ggplot2", "scDiagnostics", "SingleCellExperiment",
    "SummarizedExperiment", "Matrix", "GGally", "ggridges", "isotree",
    "transport", "ranger", "igraph", "scran", "scater"
)

cat(R.version.string, "\n")
cat("Rscript:  ", file.path(R.home("bin"), "Rscript"), "\n")
cat("library:  ", .libPaths()[1], "\n\n")

missing <- REQUIRED[!vapply(REQUIRED, requireNamespace, logical(1), quietly = TRUE)]
have <- tryCatch(as.character(utils::packageVersion("scDiagnostics")),
                 error = function(e) NA_character_)
too_old <- !is.na(have) && utils::compareVersion(have, MIN_SCDIAG) < 0

if (length(missing) || too_old) {
    cat("This R installation cannot run the app.\n\n")
    if (length(missing)) {
        cat("Missing packages:", paste(missing, collapse = ", "), "\n")
    }
    if (too_old) {
        cat("scDiagnostics", have, "is installed, but the app needs", MIN_SCDIAG, "or newer.\n")
    }
    # More than one R on the machine is the usual cause, so point at the
    # others rather than only listing what to install. Sibling directories of
    # this R's own home, and only those that really contain an Rscript.
    canon <- function(p) tryCatch(normalizePath(p, winslash = "/", mustWork = FALSE),
                                  error = function(e) p)
    rscript_in <- function(dir) {
        for (exe in c(file.path(dir, "bin", "x64", "Rscript.exe"),
                      file.path(dir, "bin", "Rscript.exe"),
                      file.path(dir, "bin", "Rscript"))) {
            if (file.exists(exe)) return(canon(exe))
        }
        NULL
    }
    siblings <- tryCatch(list.dirs(dirname(R.home()), recursive = FALSE),
                         error = function(e) character(0))
    others <- Filter(Negate(is.null), lapply(siblings, function(d) {
        if (identical(canon(d), canon(R.home()))) return(NULL)
        rscript_in(d)
    }))
    if (length(others)) {
        cat("\nOther R installations on this machine:\n")
        for (o in others) cat("  ", o, "\n")
        cat("\nIf one of those has the packages, run the app with it directly:\n")
        cat('  & "', others[[length(others)]], '" run.R\n', sep = "")
    }
    cat("\nOtherwise install what is missing:\n")
    cat('  install.packages("BiocManager")\n')
    cat('  BiocManager::install(c(', paste0('"', c(missing, if (too_old) "ccb-hms/scDiagnostics"), '"',
                                            collapse = ", "), '))\n', sep = "")
    quit(status = 1)
}

cat("scDiagnostics", have, "\n")
cat("Starting the app; it should open in your browser.\n")
cat("Press Ctrl+C here to stop it.\n\n")

shiny::runApp("app", launch.browser = TRUE)
