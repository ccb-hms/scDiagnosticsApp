# Deploy to Posit Connect.
#
# Run from the app directory:  Rscript deploy.R
#
# The manifest is regenerated from the installed library every time, so it can
# never drift from the code the way a hand-maintained one does. Run this on a
# machine whose scDiagnostics is the version you intend to deploy: the manifest
# records what is installed, not what you meant.

suppressPackageStartupMessages(library(rsconnect))

SERVER <- Sys.getenv("SCDIAG_CONNECT_SERVER", "ccb.connect.hms.harvard.edu")
APP_NAME <- Sys.getenv("SCDIAG_CONNECT_APP", "scDiagnosticsApp")

for (f in list.files("R", pattern = "[.][Rr]$", full.names = TRUE)) source(f)

msg <- check_package_version()
if (!is.null(msg)) stop(msg, call. = FALSE)
cat(sprintf("scDiagnostics %s\n", utils::packageVersion("scDiagnostics")))

if (!dir.exists("inst/precomputed")) {
    cat("No precomputed cache; building one (this takes a few minutes).\n")
    system2(file.path(R.home("bin"), "Rscript"), "data-raw/precompute.R")
}

cat("Regenerating manifest.json\n")
options(repos = BiocManager::repositories(),
        # The app tracks a package version newer than the current Bioconductor
        # release, which renv refuses to snapshot without this.
        renv.config.snapshot.validate = FALSE)
rsconnect::writeManifest(appDir = ".", appPrimaryDoc = "app.R")

accounts <- rsconnect::accounts()
if (!nrow(accounts)) {
    stop("No Connect account configured. Run rsconnect::connectApiUser() first.", call. = FALSE)
}

cat(sprintf("Deploying %s to %s\n", APP_NAME, SERVER))
rsconnect::deployApp(
    appDir = ".",
    appName = APP_NAME,
    appPrimaryDoc = "app.R",
    server = SERVER,
    forceUpdate = TRUE,
    # tests/ and data-raw/ are excluded by .rscignore.
    lint = FALSE
)
