# Small helpers shared across the app. No side effects at load time.

# The oldest scDiagnostics this app is known to work against.
#
# The app calls the package directly and shows its results as the package's
# own, so running against a version whose behaviour differs would quietly
# misattribute results. 1.7.6 is the first version carrying the fix for
# detectAnomaly() on a reference cell type absent from the query, which the
# app's own presets hit.
SCDIAG_MIN_VERSION <- "1.7.6"

#' Check the installed package version, returning a message or NULL
check_package_version <- function(min_version = SCDIAG_MIN_VERSION) {
    have <- tryCatch(utils::packageVersion(PKG), error = function(e) NULL)
    if (is.null(have)) {
        return(sprintf("The %s package is not installed on this server.", PKG))
    }
    if (utils::compareVersion(as.character(have), min_version) < 0) {
        return(sprintf(paste(
            "This app needs %s %s or newer, but %s is installed. Results would",
            "not match the published package. Install a newer version from",
            "Bioconductor devel or from github.com/ccb-hms/%s."
        ), PKG, min_version, have, PKG))
    }
    NULL
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# Blank-safe: treats "", NA and character(0) as absent too.
`%|%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L) return(y)
    if (length(x) == 1L && (is.na(x) || identical(x, ""))) return(y)
    x
}

not_null <- function(x) x[!vapply(x, is.null, logical(1))]

#' Human-readable object size
size_mb <- function(x) round(as.numeric(utils::object.size(x)) / 1024^2, 1)

#' Format a count with thousands separators
n_fmt <- function(x) formatC(x, format = "d", big.mark = ",")

#' Format a proportion as a percentage string
pct <- function(x, digits = 1) {
    if (!length(x) || all(is.na(x))) return("–")
    paste0(formatC(100 * x, format = "f", digits = digits), "%")
}

#' Truncate a character vector for display, appending "+N more"
abbrev <- function(x, n = 4) {
    x <- as.character(x)
    if (length(x) <= n) return(paste(x, collapse = ", "))
    paste0(paste(x[seq_len(n)], collapse = ", "), " +", length(x) - n, " more")
}

#' Does a colData column look like a cell type annotation?
#'
#' Returns a score; higher is more likely. Used to rank candidate columns so the
#' user is offered sensible defaults rather than every column in the object.
celltype_score <- function(values, n_cells) {
    if (is.numeric(values) && !is.factor(values)) return(-1)
    v <- as.character(values)
    v <- v[!is.na(v)]
    if (!length(v)) return(-1)
    k <- length(unique(v))
    if (k < 2 || k > 200) return(-1)
    # Favour a moderate number of levels, each with several cells.
    coverage <- length(v) / n_cells
    level_penalty <- if (k > 60) 0.3 else 1
    base <- coverage * level_penalty * (1 - abs(log10(k) - log10(12)) / 3)
    max(base, 0.01)
}

#' Rank the colData columns of an SCE by how cell-type-like they look
candidate_celltype_cols <- function(sce) {
    cd <- SummarizedExperiment::colData(sce)
    if (!ncol(cd)) return(character(0))
    scores <- vapply(
        seq_len(ncol(cd)),
        function(i) celltype_score(cd[[i]], nrow(cd)),
        numeric(1)
    )
    keep <- which(scores > 0)
    if (!length(keep)) return(character(0))
    colnames(cd)[keep[order(scores[keep], decreasing = TRUE)]]
}

#' Numeric colData columns, for QC metrics and annotation scores
numeric_cols <- function(sce) {
    cd <- SummarizedExperiment::colData(sce)
    if (!ncol(cd)) return(character(0))
    colnames(cd)[vapply(seq_len(ncol(cd)), function(i) is.numeric(cd[[i]]), logical(1))]
}

#' Does this object carry a PCA that scDiagnostics can project onto?
#'
#' Mirrors the validity checks in scDiagnostics::processPCA so the app can warn
#' before a diagnostic fails rather than after.
has_valid_pca <- function(sce, assay_name = "logcounts") {
    rd <- SingleCellExperiment::reducedDims(sce)
    if (!"PCA" %in% names(rd)) return(FALSE)
    pca <- SingleCellExperiment::reducedDim(sce, "PCA")
    at <- attributes(pca)
    if (is.null(at[["rotation"]]) || is.null(at[["percentVar"]])) return(FALSE)
    if (nrow(pca) != ncol(sce)) return(FALSE)
    rot_genes <- rownames(at[["rotation"]])
    if (is.null(rot_genes)) return(FALSE)
    all(rot_genes %in% rownames(sce))
}

#' Subset cells while keeping the attributes scDiagnostics needs
#'
#' Subsetting a SingleCellExperiment by column silently drops the attributes
#' carried on each reducedDim - including the PCA `rotation` and `percentVar`
#' that every projection depends on. Any code path that subsets must restore
#' them, or downstream functions fail with an obscure error about rotation
#' matrices. This is the only safe way to subset cells in the app.
subset_cells <- function(sce, j) {
    keep <- SingleCellExperiment::reducedDims(sce)
    saved <- lapply(keep, function(m) {
        at <- attributes(m)
        at[setdiff(names(at), c("dim", "dimnames", "class"))]
    })
    out <- sce[, j]
    for (nm in names(saved)) {
        if (!length(saved[[nm]])) next
        m <- SingleCellExperiment::reducedDim(out, nm)
        for (a in names(saved[[nm]])) attr(m, a) <- saved[[nm]][[a]]
        SingleCellExperiment::reducedDim(out, nm) <- m
    }
    out
}

#' Stratified subsample of an SCE, keeping cell type proportions
#'
#' The package's own functions accept max_cells_* arguments that downsample
#' internally per call. This is for the upload path, where we want to shrink the
#' object once, up front, so that every later call is cheap and memory stays
#' bounded for the whole session.
subsample_sce <- function(sce, max_cells, cell_type_col = NULL, seed = 1L) {
    n <- ncol(sce)
    if (is.null(max_cells) || n <= max_cells) return(sce)
    withr_seed <- function(code) {
        old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
        set.seed(seed)
        on.exit({
            if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv)
        }, add = TRUE)
        code
    }
    idx <- withr_seed({
        if (!is.null(cell_type_col) && cell_type_col %in% colnames(SummarizedExperiment::colData(sce))) {
            groups <- as.character(SummarizedExperiment::colData(sce)[[cell_type_col]])
            groups[is.na(groups)] <- "__NA__"
            split_idx <- split(seq_len(n), groups)
            # Proportional allocation, but never drop a group below 2 cells if
            # it has them: small populations are exactly what diagnostics are for.
            target <- pmax(2, round(max_cells * lengths(split_idx) / n))
            target <- pmin(target, lengths(split_idx))
            # If rounding overshoots, trim the largest groups first.
            while (sum(target) > max_cells && any(target > 2)) {
                big <- which.max(target)
                target[big] <- target[big] - 1L
            }
            sort(unlist(Map(function(ix, k) if (length(ix) <= k) ix else sample(ix, k),
                            split_idx, target), use.names = FALSE))
        } else {
            sort(sample(n, max_cells))
        }
    })
    subset_cells(sce, idx)
}

#' Coerce a value to the integer vector the package expects for pc_subset
as_pc_subset <- function(x) {
    x <- suppressWarnings(as.integer(x))
    x <- x[!is.na(x) & x > 0]
    if (!length(x)) 1:5 else sort(unique(x))
}

#' A short, stable hash of anything
hash_of <- function(...) substr(digest::digest(list(...), algo = "xxhash64"), 1, 16)

#' Wrap text as a Bootstrap alert
alert <- function(..., type = c("info", "warning", "danger", "success")) {
    type <- match.arg(type)
    htmltools::div(class = paste0("alert alert-", type, " mb-3"), ...)
}

#' Standard empty-state panel
empty_state <- function(message, icon = "hourglass-split") {
    htmltools::div(
        class = "sc-empty",
        bsicons::bs_icon(icon, size = "2rem", class = "text-muted mb-2"),
        htmltools::p(class = "text-muted mb-0", message)
    )
}

#' Keep these outputs rendering even while their tab is hidden
#'
#' Shiny suspends outputs it believes are invisible and resumes them when they
#' are shown. Inside a bslib navbar page that resume does not reliably fire, so
#' a tab can be displayed with every panel in it permanently blank. Disabling
#' suspension makes rendering deterministic.
#'
#' This does not reintroduce the cost suspension exists to avoid: each panel
#' gates its own computation on an `active` reactive, so a hidden panel still
#' does no work - it just renders an empty placeholder instead of nothing.
always_render <- function(output, names) {
    for (nm in names) {
        try(shiny::outputOptions(output, nm, suspendWhenHidden = FALSE), silent = TRUE)
    }
    invisible(NULL)
}

#' Render an error from a diagnostic as a readable panel rather than a red blob
error_panel <- function(msg, context = NULL) {
    htmltools::div(
        class = "alert alert-danger",
        htmltools::tags$strong("This diagnostic could not be computed."),
        htmltools::tags$pre(class = "sc-error-detail mb-2", msg),
        if (!is.null(context)) htmltools::div(class = "small text-muted", context)
    )
}
