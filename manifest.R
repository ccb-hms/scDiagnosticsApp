# Regenerate app/manifest.json for git-backed deployment.
#
#   Rscript manifest.R
#
# Connect clones this repository and restores packages from the manifest, so
# the manifest has to name versions Connect can actually obtain. It records
# whatever is installed on the machine that generates it, which makes one
# failure mode easy to walk into: a package installed from a local source tree
# is written down as if it came from a repository that has never heard of that
# version. The deploy then fails on the server, long after the mistake.
#
# So this regenerates the manifest and then checks every package in it is
# reachable at the version recorded, before you commit.

if (!file.exists(file.path("app", "app.R"))) {
    stop("Run this from the repository root:  Rscript manifest.R", call. = FALSE)
}
for (p in c("rsconnect", "BiocManager", "jsonlite")) {
    if (!requireNamespace(p, quietly = TRUE)) {
        stop("Missing package '", p, "'. Install it with install.packages(\"", p, "\").",
             call. = FALSE)
    }
}

cat(R.version.string, "\n")
cat("library: ", .libPaths()[1], "\n\n")

repos <- BiocManager::repositories()
options(repos = repos,
        # The app may track a package newer than the current Bioconductor
        # release, which renv's pre-flight check refuses to snapshot.
        renv.config.snapshot.validate = FALSE)

cat("Repositories the manifest will point at:\n")
for (i in seq_along(repos)) cat(sprintf("  %-14s %s\n", names(repos)[i], repos[i]))
cat("\nRegenerating app/manifest.json ...\n")
rsconnect::writeManifest(appDir = "app", appPrimaryDoc = "app.R", verbose = FALSE)

m <- jsonlite::fromJSON("app/manifest.json", simplifyVector = FALSE)
pkgs <- m$packages
cat(sprintf("\n%d packages, %d files, R %s\n", length(pkgs), length(m$files), m$platform))

# --- can Connect actually get these? ---------------------------------------
#
# Build the set of (package, version) each repository can serve, then compare.
cat("\nChecking every package is obtainable at the version recorded ...\n")
available <- list()
for (u in unique(vapply(pkgs, function(p) p$Repository %||% "", character(1)))) {
    if (!nzchar(u)) next
    available[[u]] <- tryCatch(
        available.packages(contriburl = contrib.url(u, "source")),
        error = function(e) NULL
    )
}

problems <- list()
for (nm in names(pkgs)) {
    p <- pkgs[[nm]]
    want <- p$description$Version
    src <- p$Source %||% "?"
    repo <- p$Repository %||% ""
    if (identical(src, "github") || grepl("^Remote", src, ignore.case = TRUE)) next
    ap <- available[[repo]]
    if (is.null(ap) || !nm %in% rownames(ap)) {
        problems[[nm]] <- list(want = want, src = src, repo = repo, has = "not in this repository")
        next
    }
    has <- unname(ap[nm, "Version"])
    if (utils::compareVersion(has, want) < 0) {
        problems[[nm]] <- list(want = want, src = src, repo = repo, has = has)
    }
}

if (!length(problems)) {
    cat("  all packages are obtainable.\n")
    cat("\nmanifest.json is ready. Commit it:\n")
    cat("  git add app/manifest.json && git commit -m \"Update deployment manifest\"\n")
    quit(status = 0)
}

cat("\n", strrep("-", 70), "\n", sep = "")
cat("PROBLEM: the manifest names versions the recorded repository cannot serve.\n")
cat("Connect will fail to restore these when it deploys.\n\n")
for (nm in names(problems)) {
    x <- problems[[nm]]
    cat(sprintf("  %-22s manifest wants %-10s repository has %s\n", nm, x$want, x$has))
    cat(sprintf("  %-22s   %s\n", "", x$repo))
}
cat("\nUsually this means the package was installed from a local source tree or\n")
cat("a development branch, so the version exists only on this machine.\n\n")
cat("Two ways out:\n")
cat("  1. Publish the version to the repository the manifest points at, wait for\n")
cat("     it to appear, then run this script again.\n")
cat("  2. Install it from a source Connect can reach, so the manifest records\n")
cat("     that source instead of a repository, for example:\n")
cat('       BiocManager::install("ccb-hms/scDiagnostics")\n')
cat("     then run this script again. Connect must be configured to allow\n")
cat("     package installation from GitHub for this to work.\n")
quit(status = 1)
