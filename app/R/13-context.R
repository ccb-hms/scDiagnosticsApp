# The data context.
#
# One object carrying everything a diagnostic needs to know about the data:
# which reference, which query, which annotation columns, which assay, which
# components, which cell types. Chosen once and inherited everywhere, which is
# what removes the per-tab re-selection the previous app required.

#' Assemble a data context from its raw parts
#'
#' Derived fields (shared cell types, available components) are computed here
#' so that every consumer sees the same answer.
build_context <- function(ref = NULL, query = NULL,
                          ref_desc = NULL, query_desc = NULL,
                          ref_col = NULL, query_col = NULL,
                          cell_types = NULL, assay = "logcounts",
                          pc_subset = 1:5, focus_which = "query",
                          preset_id = NULL) {

    ref_types <- if (!is.null(ref) && !is.null(ref_col) &&
                     ref_col %in% colnames(SummarizedExperiment::colData(ref))) {
        sort(unique(stats::na.omit(as.character(ref[[ref_col]]))))
    } else character(0)

    query_types <- if (!is.null(query) && !is.null(query_col) &&
                       query_col %in% colnames(SummarizedExperiment::colData(query))) {
        sort(unique(stats::na.omit(as.character(query[[query_col]]))))
    } else character(0)

    shared <- intersect(ref_types, query_types)
    # Offer shared types first: a diagnostic comparing reference against query
    # can only act on types that exist in both. Reference-only types remain
    # selectable because some diagnostics report on them meaningfully.
    available <- unique(c(shared, setdiff(ref_types, shared), setdiff(query_types, shared)))

    n_pcs <- pca_ncol(ref)
    if (!is.null(query)) n_pcs <- min(n_pcs, pca_ncol(query), na.rm = TRUE)
    if (!is.finite(n_pcs) || n_pcs < 2) n_pcs <- 20

    focus <- if (identical(focus_which, "ref")) ref else query
    focus_col <- if (identical(focus_which, "ref")) ref_col else query_col
    if (is.null(focus)) {
        focus <- ref %||% query
        focus_col <- ref_col %||% query_col
    }

    structure(list(
        ref = ref, query = query,
        ref_desc = ref_desc, query_desc = query_desc,
        ref_col = ref_col, query_col = query_col,
        cell_types = if (length(cell_types)) as.character(cell_types) else NULL,
        assay = assay,
        pc_subset = as_pc_subset(pc_subset),
        n_pcs = n_pcs,
        focus = focus, focus_col = focus_col, focus_which = focus_which,
        ref_types = ref_types, query_types = query_types,
        shared_types = shared,
        cell_types_available = available,
        preset_id = preset_id,
        ready = !is.null(ref) && !is.null(ref_col)
    ), class = "sc_context")
}

#' Number of components stored on an object, or NA
pca_ncol <- function(sce) {
    if (is.null(sce)) return(NA_integer_)
    rd <- SingleCellExperiment::reducedDims(sce)
    if (!"PCA" %in% names(rd)) return(NA_integer_)
    ncol(SingleCellExperiment::reducedDim(sce, "PCA"))
}

#' Assays common to both objects, for the assay picker
common_assays <- function(ref, query) {
    a <- if (is.null(ref)) character(0) else SummarizedExperiment::assayNames(ref)
    if (!is.null(query)) a <- intersect(a, SummarizedExperiment::assayNames(query))
    # logcounts first: every scDiagnostics default assumes log-normalised values.
    c(intersect("logcounts", a), setdiff(a, "logcounts"))
}

#' A short signature of what the context *is*, ignoring display settings
#'
#' The context reactive fires several times as a page's inputs arrive from the
#' client, and again whenever a slider moves. Code that needs to discard
#' results when the underlying data changes must key on this instead, or it
#' will throw away answers it has only just computed.
context_signature <- function(ctx) {
    paste(
        if (is.null(ctx$ref_desc)) "-" else descriptor_id(ctx$ref_desc),
        if (is.null(ctx$query_desc)) "-" else descriptor_id(ctx$query_desc),
        ctx$ref_col %||% "-", ctx$query_col %||% "-", ctx$assay %||% "-",
        sep = "|"
    )
}

#' A one-line description of the current context, for headers and reports
context_caption <- function(ctx) {
    if (!isTRUE(ctx$ready)) return("No data selected")
    parts <- c(
        sprintf("Reference: %s (%s cells, %s)",
                ctx$ref_desc$label %||% "reference", n_fmt(ncol(ctx$ref)), ctx$ref_col)
    )
    if (!is.null(ctx$query)) {
        parts <- c(parts, sprintf("Query: %s (%s cells, %s)",
                                  ctx$query_desc$label %||% "query",
                                  n_fmt(ncol(ctx$query)), ctx$query_col))
    }
    paste(parts, collapse = "  ·  ")
}

# ---------------------------------------------------------------------------
# Upload validation
# ---------------------------------------------------------------------------

#' Inspect an uploaded object and report what the app found
#'
#' Returns a list of checks, each with a status of "ok", "warn" or "fail", so
#' the user sees a report card before anything runs rather than a stack trace
#' afterwards.
inspect_upload <- function(sce, other = NULL, role = "reference") {
    checks <- list()
    add <- function(status, label, detail) {
        checks[[length(checks) + 1]] <<- list(status = status, label = label, detail = detail)
    }

    if (!methods::is(sce, "SingleCellExperiment")) {
        add("fail", "Object type",
            sprintf(paste("This is a %s. The app needs a SingleCellExperiment",
                          "(a SpatialExperiment also works, since it extends one)."),
                    paste(class(sce), collapse = "/")))
        return(checks)
    }
    add("ok", "Object type",
        sprintf("%s, %s genes × %s cells",
                class(sce)[1], n_fmt(nrow(sce)), n_fmt(ncol(sce))))

    assays <- SummarizedExperiment::assayNames(sce)
    if (!length(assays)) {
        add("fail", "Assays", "No assays found. At least one expression matrix is required.")
    } else if ("logcounts" %in% assays) {
        add("ok", "Assays", paste(assays, collapse = ", "))
    } else {
        add("warn", "Assays",
            paste0("Found ", paste(assays, collapse = ", "),
                   " but no 'logcounts'. Diagnostics assume log-normalised values; ",
                   "run scuttle::logNormCounts() first if these are raw counts."))
    }

    cols <- candidate_celltype_cols(sce)
    if (!length(cols)) {
        add("fail", "Cell type annotation",
            paste("No colData column looks like a cell type annotation (a character",
                  "or factor column with between 2 and 200 levels)."))
    } else {
        add("ok", "Cell type annotation",
            sprintf("%d candidate column%s: %s", length(cols),
                    if (length(cols) == 1) "" else "s", abbrev(cols, 5)))
    }

    if (has_valid_pca(sce)) {
        add("ok", "PCA",
            sprintf("%d components stored, with rotation and variance explained",
                    pca_ncol(sce)))
    } else if ("PCA" %in% SingleCellExperiment::reducedDimNames(sce)) {
        add("warn", "PCA",
            paste("A PCA is present but lacks the rotation matrix or variance the",
                  "projection needs. The app will recompute it."))
    } else {
        add("warn", "PCA",
            "No PCA found. The app will compute one from highly variable genes when you continue.")
    }

    if (!is.null(other) && methods::is(other, "SingleCellExperiment")) {
        shared <- length(intersect(rownames(sce), rownames(other)))
        frac <- shared / max(1, min(nrow(sce), nrow(other)))
        # Gene overlap is the most common silent failure: the projection needs
        # the reference's rotation genes to exist in the query.
        if (frac >= 0.9) {
            add("ok", "Gene overlap with the other dataset",
                sprintf("%s genes shared (%s of the smaller object)", n_fmt(shared), pct(frac)))
        } else if (frac >= 0.5) {
            add("warn", "Gene overlap with the other dataset",
                sprintf("Only %s genes shared (%s of the smaller object). Diagnostics will run on the intersection.",
                        n_fmt(shared), pct(frac)))
        } else {
            add("fail", "Gene overlap with the other dataset",
                sprintf(paste("Just %s genes shared (%s). These two objects are unlikely to be",
                              "comparable — check that both use the same gene identifiers."),
                        n_fmt(shared), pct(frac)))
        }
    }

    checks
}

upload_blocked <- function(checks) {
    any(vapply(checks, function(c) identical(c$status, "fail"), logical(1)))
}
