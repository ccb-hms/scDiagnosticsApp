# Two-tier result cache.
#
# Preset results are identical for every visitor, so they belong in a disk
# cache shared by all sessions and pre-warmed at build time. Results derived
# from uploaded data are private, so they stay in a per-session memory cache
# that dies with the session. Keeping these apart is a correctness requirement,
# not an optimisation: a single shared cache would serve one user's results to
# another whenever two keys collided on anything but the data identity.

CACHE_DIR <- Sys.getenv("SCDIAG_CACHE_DIR", unset = file.path(tempdir(), "scdiag-cache"))
CACHE_MAX_MB <- as.numeric(Sys.getenv("SCDIAG_CACHE_MAX_MB", unset = "1024"))
SESSION_CACHE_MAX_MB <- as.numeric(Sys.getenv("SCDIAG_SESSION_CACHE_MAX_MB", unset = "256"))

#' Shared, cross-session cache for preset results
shared_cache <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
            seed_cache(CACHE_DIR)
            cache <<- cachem::cache_disk(
                dir = CACHE_DIR,
                max_size = CACHE_MAX_MB * 1024^2,
                evict = "lru",
                warn_ref_objects = FALSE
            )
        }
        cache
    }
})

#' Copy the shipped, precomputed preset results into the runtime cache
#'
#' Built by data-raw/precompute.R. Copying rather than pointing at the shipped
#' directory keeps the app working where the deployment directory is read-only,
#' which it is on most servers.
seed_cache <- function(dir) {
    src <- "inst/precomputed"
    if (!dir.exists(src)) return(invisible(FALSE))
    if (length(list.files(dir, pattern = "^k"))) return(invisible(FALSE))
    files <- list.files(src, pattern = "^k", full.names = TRUE)
    if (!length(files)) return(invisible(FALSE))
    ok <- tryCatch(all(file.copy(files, dir, overwrite = FALSE)),
                   error = function(e) FALSE)
    if (isTRUE(ok)) {
        message(sprintf("Seeded the result cache with %d precomputed results.", length(files)))
    }
    invisible(ok)
}

#' Private, per-session cache for results derived from uploads
new_session_cache <- function() {
    cachem::cache_mem(max_size = SESSION_CACHE_MAX_MB * 1024^2, evict = "lru")
}

#' Build a cache key from everything that can change a result.
#'
#' Data enters the key by identity (`descriptor_id`), never by value, so the
#' key is cheap to compute even for a large object. The package version is
#' folded in via descriptor_id for presets and explicitly here for uploads, so
#' upgrading scDiagnostics invalidates every cached result.
cache_key <- function(entry_id, ctx, args, extra = NULL) {
    payload <- list(
        entry = entry_id,
        pkg = as.character(utils::packageVersion(PKG)),
        ref = if (!is.null(ctx$ref_desc)) descriptor_id(ctx$ref_desc) else NULL,
        query = if (!is.null(ctx$query_desc)) descriptor_id(ctx$query_desc) else NULL,
        ref_col = ctx$ref_col,
        query_col = ctx$query_col,
        cell_types = sort(as.character(ctx$cell_types %||% character(0))),
        assay = ctx$assay,
        pc_subset = ctx$pc_subset,
        args = args[order(names(args))],
        extra = extra
    )
    paste0("k", digest::digest(payload, algo = "xxhash64"))
}

#' Look up `key`, otherwise evaluate `expr` and store it.
#'
#' `shareable` picks the tier. Errors are never cached: a transient failure
#' should not be remembered as the answer.
with_cache <- function(key, cache, shareable, expr, enabled = TRUE) {
    if (!enabled || is.null(cache)) return(force(expr))
    hit <- cache$get(key)
    if (!cachem::is.key_missing(hit)) {
        attr(hit, "sc_cache_hit") <- TRUE
        return(hit)
    }
    val <- force(expr)
    # Only cache successful results.
    if (!inherits(val, "sc_failure")) {
        ok <- tryCatch({
            cache$set(key, val)
            TRUE
        }, error = function(e) FALSE)
        if (!ok) warning("Could not write to cache; continuing uncached.", call. = FALSE)
    }
    attr(val, "sc_cache_hit") <- FALSE
    val
}

#' Choose the right cache tier for a context
cache_for <- function(ctx, session_cache) {
    shareable <- is_shareable(ctx$ref_desc) &&
        (is.null(ctx$query_desc) || is_shareable(ctx$query_desc))
    list(cache = if (shareable) shared_cache() else session_cache, shareable = shareable)
}

#' A failed computation, carried as a value rather than thrown.
#'
#' Diagnostics fail for ordinary reasons (a cell type with too few cells, a
#' gene absent from one dataset). Those should reach the user as an explanation
#' in the right panel, not as a red bar over the whole app.
sc_failure <- function(message, call = NULL) {
    structure(list(message = message, call = call), class = "sc_failure")
}

is_failure <- function(x) inherits(x, "sc_failure")

#' Run `expr`, converting any error or warning into a value.
#'
#' Warnings are collected rather than suppressed: scDiagnostics warns about
#' skipped cell types and recomputed PCA, and those messages are exactly what a
#' user needs to see.
capture_run <- function(expr) {
    warnings <- character(0)
    value <- withCallingHandlers(
        tryCatch(force(expr), error = function(e) sc_failure(conditionMessage(e), conditionCall(e))),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        },
        message = function(m) {
            warnings <<- c(warnings, sub("\n$", "", conditionMessage(m)))
            invokeRestart("muffleMessage")
        }
    )
    attr(value, "sc_warnings") <- unique(warnings)
    value
}
