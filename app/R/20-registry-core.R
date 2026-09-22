# The diagnostic registry.
#
# One entry describes one scDiagnostics function completely enough that the
# generic module can build its controls, run it, cache it, plot it, tabulate
# it, explain it and print the R code that reproduces it. Adding a function to
# the app is adding an entry, not writing a tab.

CATEGORIES <- list(
    visualization = list(
        id = "visualization",
        label = "Visualization",
        icon = "bar-chart",
        blurb = "See where the query sits relative to the reference, cell type by cell type."
    ),
    alignment = list(
        id = "alignment",
        label = "Alignment & statistics",
        icon = "rulers",
        blurb = "Quantify how well the two datasets line up, and test whether differences are real."
    ),
    anomaly = list(
        id = "anomaly",
        label = "Anomaly & distances",
        icon = "exclamation-diamond",
        blurb = "Flag individual cells that do not look like their assigned reference population."
    ),
    markers = list(
        id = "markers",
        label = "Marker genes",
        icon = "activity",
        blurb = "Ask which genes differ, and whether the genes that define a cell type agree."
    ),
    qc = list(
        id = "qc",
        label = "QC & annotation scores",
        icon = "clipboard-check",
        blurb = "Relate annotation confidence to quality control metrics."
    ),
    utility = list(
        id = "utility",
        label = "Utilities",
        icon = "tools",
        blurb = "The building blocks the other diagnostics are made of."
    )
)

#' Define a diagnostic
#'
#' @param id stable identifier, used in cache keys and URLs
#' @param fn name of the scDiagnostics function, as a string, so that the call
#'   can be deparsed for the reproducible-code panel and the help looked up
#' @param data_map named list mapping the function's formal arguments to fields
#'   of the data context: "ref", "query", "ref_col", "query_col", "focus",
#'   "focus_col", "cell_types", "pc_subset", "assay"
#' @param params extra user-facing parameters
#' @param plot_params parameters passed to the result's plot() method
#' @param prepare optional function(ctx, args) returning the final argument
#'   list, for the few functions whose inputs are not SCE objects
#' @param tables function(result) returning a named list of data frames
#' @param summarise function(result) returning a list of headline figures
#' @param reading one-paragraph guidance on interpreting the output
#' @param heavy run in a background worker behind a task button
diag_entry <- function(id, fn, label, category, blurb,
                       tier = 2L,
                       question = NULL,
                       data_map = list(),
                       params = list(),
                       plot_params = list(),
                       prepare = NULL,
                       repro_prelude = NULL,
                       repro_names = NULL,
                       plot_fn = NULL,
                       tables = NULL,
                       summarise = NULL,
                       reading = NULL,
                       available = NULL,
                       needs = c("ref", "query"),
                       heavy = FALSE,
                       returns = c("object", "plot", "value", "sce"),
                       plot_height = 520) {
    returns <- match.arg(returns)
    stopifnot(category %in% names(CATEGORIES))
    structure(list(
        id = id, fn = fn, label = label, category = category, blurb = blurb,
        tier = as.integer(tier), question = question,
        data_map = data_map, params = params, plot_params = plot_params,
        prepare = prepare, repro_prelude = repro_prelude, repro_names = repro_names,
        plot_fn = plot_fn, tables = tables,
        summarise = summarise, reading = reading, available = available, needs = needs,
        heavy = heavy, returns = returns, plot_height = plot_height
    ), class = "sc_diag_entry")
}

# --- argument-map shorthands ----------------------------------------------

#' Functions taking query and reference plus both cell type columns
map_pair <- function(pcs = TRUE, cell_types = TRUE, assay = TRUE) {
    m <- list(
        query_data = "query", reference_data = "ref",
        query_cell_type_col = "query_col", ref_cell_type_col = "ref_col"
    )
    if (cell_types) m$cell_types <- "cell_types"
    if (pcs) m$pc_subset <- "pc_subset"
    if (assay) m$assay_name <- "assay"
    m
}

#' Functions whose first argument is the reference and whose query is optional
map_ref_first <- function(pcs = TRUE, cell_types = TRUE, assay = TRUE) {
    m <- list(
        reference_data = "ref", query_data = "query",
        ref_cell_type_col = "ref_col", query_cell_type_col = "query_col"
    )
    if (cell_types) m$cell_types <- "cell_types"
    if (pcs) m$pc_subset <- "pc_subset"
    if (assay) m$assay_name <- "assay"
    m
}

#' Functions operating on a single object, chosen by the user
map_single <- function(arg = "sce_object", col = "cell_type_col",
                       cell_types = TRUE, assay = TRUE) {
    m <- stats::setNames(list("focus"), arg)
    # Some single-object functions take no cell type column at all.
    if (!is.null(col)) m[[col]] <- "focus_col"
    if (cell_types) m$cell_types <- "cell_types"
    if (assay) m$assay_name <- "assay"
    m
}

#' Does this entry act on a single, user-chosen object?
#'
#' True either because the argument map references the focus, or because the
#' entry declares it in `needs` while building its inputs in prepare().
uses_focus <- function(entry) {
    any(unlist(entry$data_map) %in% c("focus", "focus_col")) ||
        "focus" %in% entry$needs
}

# --- resolving an entry into a call ----------------------------------------

#' Turn a context plus parameter values into the argument list for the call
build_args <- function(entry, ctx, param_values) {
    args <- list()
    for (arg in names(entry$data_map)) {
        field <- entry$data_map[[arg]]
        args[[arg]] <- switch(field,
            ref = ctx$ref,
            query = ctx$query,
            ref_col = ctx$ref_col,
            query_col = ctx$query_col,
            focus = ctx$focus,
            focus_col = ctx$focus_col,
            cell_types = ctx$cell_types,
            pc_subset = ctx$pc_subset,
            assay = ctx$assay,
            stop("Unknown data_map field: ", field)
        )
    }
    # Parameter values win over context defaults where ids collide, which is
    # how a panel overrides the global PC range or cell type selection.
    for (nm in names(param_values)) args[[nm]] <- param_values[[nm]]
    if (!is.null(entry$prepare)) args <- entry$prepare(ctx, args)
    # Drop NULLs the function does not default to NULL itself; keeping an
    # explicit NULL would override a real default.
    formals_fn <- names(formals(get(entry$fn, envir = asNamespace(PKG))))
    args <- args[names(args) %in% formals_fn]
    drop <- vapply(args, is.null, logical(1))
    # cell_types = NULL is meaningful ("use all"), so keep it.
    keep_null <- names(args) %in% c("cell_types")
    args[!drop | keep_null]
}

#' Run an entry, returning the result or an sc_failure
run_entry <- function(entry, args) {
    fn <- get(entry$fn, envir = asNamespace(PKG))
    capture_run(do.call(fn, args))
}

#' Which context pieces must exist before this entry can run?
entry_ready <- function(entry, ctx) {
    if (is.null(ctx$ref) && "ref" %in% entry$needs) {
        return("Choose a reference dataset first.")
    }
    if (is.null(ctx$query) && "query" %in% entry$needs) {
        return("This diagnostic compares a query against the reference. Choose a query dataset.")
    }
    if (uses_focus(entry) && is.null(ctx$focus)) {
        return("Choose a dataset first.")
    }
    # Entry-specific preconditions: a diagnostic that cannot work on this data
    # should say why, in its own terms, before the user presses anything.
    if (!is.null(entry$available)) {
        msg <- tryCatch(entry$available(ctx), error = function(e) NULL)
        if (!is.null(msg) && nzchar(msg)) return(msg)
    }
    NULL
}

#' Gene names behind an object's stored PCA rotation
rotation_genes <- function(sce) {
    if (is.null(sce)) return(NULL)
    rd <- SingleCellExperiment::reducedDims(sce)
    if (!"PCA" %in% names(rd)) return(NULL)
    rownames(attributes(SingleCellExperiment::reducedDim(sce, "PCA"))[["rotation"]])
}

#' Precondition: the object in focus carries at least `n` numeric colData columns
needs_numeric_cols <- function(n = 1, what = "a numeric metadata column") {
    force(n)
    force(what)
    function(ctx) {
        cols <- numeric_cols(ctx$focus)
        if (length(cols) >= n) return(NULL)
        sprintf(paste(
            "This diagnostic needs %s in colData, and the %s dataset has %s.",
            "The bone marrow preset carries QC metrics and annotation scores,",
            "or upload data whose colData holds them."
        ), what, if (identical(ctx$focus_which, "ref")) "reference" else "query",
           if (!length(cols)) "none" else sprintf("only %d", length(cols)))
    }
}

#' Precondition: reference and query PCA were fitted on the same genes
needs_matched_rotation <- function(ctx) {
    rg <- rotation_genes(ctx$ref)
    qg <- rotation_genes(ctx$query)
    if (is.null(rg) || is.null(qg)) return(NULL)
    if (setequal(rg, qg)) return(NULL)
    sprintf(paste(
        "This diagnostic compares the two datasets' PCA loadings gene by gene,",
        "so both must have been fitted on the same genes. Here the reference",
        "rotation covers %s genes and the query %s, sharing %s. Recompute both",
        "PCAs on the shared gene set to use it."
    ), n_fmt(length(rg)), n_fmt(length(qg)), n_fmt(length(intersect(rg, qg))))
}

# --- automatic tables ------------------------------------------------------

#' Extract displayable tables from an arbitrary result object
#'
#' Every Tier 2 diagnostic gets a numbers view without bespoke code: walk the
#' result, keep data frames and named numeric vectors, and skip the large
#' matrices that are internal state rather than output.
auto_tables <- function(x, max_tables = 6, max_rows = 5000) {
    out <- list()
    walk <- function(v, path) {
        if (length(out) >= max_tables) return(invisible(NULL))
        if (is.null(v) || inherits(v, c("ggplot", "igraph", "function"))) return(invisible(NULL))
        if (is.data.frame(v)) {
            if (nrow(v) > 0 && nrow(v) <= max_rows) out[[path]] <<- v
            return(invisible(NULL))
        }
        if (is.matrix(v) && is.numeric(v) && nrow(v) <= 200 && ncol(v) <= 60) {
            df <- as.data.frame(v)
            if (!is.null(rownames(v))) df <- cbind(row = rownames(v), df)
            out[[path]] <<- df
            return(invisible(NULL))
        }
        if (is.numeric(v) && is.null(dim(v)) && length(v) > 0 && length(v) <= max_rows) {
            df <- data.frame(name = names(v) %||% seq_along(v), value = as.numeric(v))
            out[[path]] <<- df
            return(invisible(NULL))
        }
        if (is.list(v)) {
            keys <- names(v) %||% as.character(seq_along(v))
            for (i in seq_along(v)) {
                walk(v[[i]], if (nzchar(path)) paste(path, keys[i], sep = " › ") else keys[i])
            }
        }
        invisible(NULL)
    }
    walk(x, "")
    out
}

#' Registry lookup helpers
registry_get <- function(id) {
    e <- REGISTRY[[id]]
    if (is.null(e)) stop("Unknown diagnostic: ", id)
    e
}

registry_by_category <- function(category) {
    Filter(function(e) e$category == category, REGISTRY)
}

registry_tier <- function(tier) {
    Filter(function(e) e$tier == tier, REGISTRY)
}
