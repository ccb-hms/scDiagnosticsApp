# Reproducible code generation.
#
# Every panel shows the R that produced what is on screen. This is close to
# free: the registry already holds the function name and the resolved argument
# list, so the call can simply be deparsed. It matters more than it costs,
# because an app that cannot hand back a script is a dead end for anyone who
# needs the result in a paper.

#' How a dataset is loaded, as R code
data_line <- function(desc, varname) {
    if (is.null(desc)) return(NULL)
    if (desc$kind == "preset") {
        sprintf('data("%s", package = "scDiagnostics")\n%s <- %s',
                desc$key, varname, desc$key)
    } else {
        sprintf('%s <- readRDS("%s")', varname, basename(desc$key))
    }
}

#' Deparse one argument value into R source
deparse_value <- function(v) {
    if (is.null(v)) return("NULL")
    if (inherits(v, c("SingleCellExperiment", "SummarizedExperiment"))) return("<data>")
    # Never inline a large object: a 500-gene character vector or a score
    # matrix would bury the call it belongs to. The entry's repro_prelude is
    # responsible for building these, and repro_names for naming them.
    if (!is.null(dim(v))) {
        return(sprintf("<%s, %s>", class(v)[1], paste(dim(v), collapse = " x ")))
    }
    if (length(v) > 12) {
        return(sprintf("<%s vector of length %d>", class(v)[1], length(v)))
    }
    if (is.integer(v) && length(v) > 2 && identical(v, seq(v[1], v[length(v)]))) {
        return(sprintf("%d:%d", v[1], v[length(v)]))
    }
    if (is.numeric(v) && length(v) > 2 && identical(as.integer(v), seq(v[1], v[length(v)]))) {
        return(sprintf("%d:%d", as.integer(v[1]), as.integer(v[length(v)])))
    }
    out <- paste(deparse(v, width.cutoff = 500L), collapse = " ")
    gsub("\\s+", " ", out)
}

#' Render a function call as formatted R source
call_source <- function(fn, args, assign_to = "result", data_names = list()) {
    if (!length(args)) return(sprintf("%s <- %s()", assign_to, fn))
    parts <- vapply(names(args), function(nm) {
        v <- args[[nm]]
        txt <- if (!is.null(data_names[[nm]])) data_names[[nm]] else deparse_value(v)
        sprintf("    %s = %s", nm, txt)
    }, character(1))
    sprintf("%s <- %s(\n%s\n)", assign_to, fn, paste(parts, collapse = ",\n"))
}

#' Full reproducible script for one diagnostic
#'
#' @param entry registry entry
#' @param ctx data context
#' @param args the argument list actually passed to the function
#' @param plot_args arguments passed to plot(), if any
repro_code <- function(entry, ctx, args, plot_args = NULL, param_values = list()) {
    lines <- c("library(scDiagnostics)", "library(SingleCellExperiment)")

    # Name the data objects the way the loading lines name them.
    data_names <- list()
    ref_var <- if (!is.null(ctx$ref_desc) && ctx$ref_desc$kind == "preset") {
        ctx$ref_desc$key
    } else "reference_data"
    query_var <- if (!is.null(ctx$query_desc) && ctx$query_desc$kind == "preset") {
        ctx$query_desc$key
    } else "query_data"

    has_sce_arg <- any(vapply(args, function(v) inherits(v, "SummarizedExperiment"), logical(1)))
    needs_ref <- any(vapply(args, function(v) is_same_sce(v, ctx$ref), logical(1)))
    needs_query <- any(vapply(args, function(v) is_same_sce(v, ctx$query), logical(1)))
    # A prepare() step may consume the object entirely (building a matrix or a
    # gene list from it), leaving no SCE in the call. The prelude still needs
    # the data, so load it.
    if (!has_sce_arg) {
        if (uses_focus(entry)) {
            if (identical(ctx$focus_which, "ref")) needs_ref <- TRUE else needs_query <- TRUE
        } else {
            needs_ref <- needs_ref || !is.null(ctx$ref)
            needs_query <- needs_query || !is.null(ctx$query)
        }
    }
    loads <- character(0)
    if (needs_ref) loads <- c(loads, data_line(ctx$ref_desc, ref_var))
    if (needs_query && !identical(ref_var, query_var)) {
        loads <- c(loads, data_line(ctx$query_desc, query_var))
    }
    if (!has_sce_arg && uses_focus(entry)) {
        focus_var <- if (identical(ctx$focus_which, "ref")) ref_var else query_var
        loads <- c(loads, sprintf("sce <- %s", focus_var))
    }
    if (length(loads)) lines <- c(lines, "", loads)

    for (nm in names(args)) {
        v <- args[[nm]]
        if (is_same_sce(v, ctx$ref)) data_names[[nm]] <- ref_var
        else if (is_same_sce(v, ctx$query)) data_names[[nm]] <- query_var
        # A prepare() step that subsets gets its own variable, named after the
        # argument, which the entry's repro_prelude is responsible for defining.
        else if (inherits(v, "SummarizedExperiment")) {
            data_names[[nm]] <- sub("_data$", "", nm)
        }
    }
    # Entries whose prepare() builds non-SCE inputs name them here, so the call
    # refers to the variables the prelude created rather than inlining them.
    data_names <- utils::modifyList(data_names, entry$repro_names %||% list())

    if (!is.null(entry$repro_prelude)) {
        pre <- tryCatch(
            entry$repro_prelude(ctx, args, param_values, ref_var, query_var),
            error = function(e) NULL
        )
        if (length(pre)) lines <- c(lines, "", pre)
    }

    lines <- c(lines, "", call_source(entry$fn, args, "result", data_names))

    if (entry$returns == "object") {
        pa <- Filter(Negate(is.null), plot_args %||% list())
        if (length(pa)) {
            parts <- vapply(names(pa), function(nm) {
                sprintf("    %s = %s", nm, deparse_value(pa[[nm]]))
            }, character(1))
            lines <- c(lines, "", sprintf("plot(\n    result,\n%s\n)", paste(parts, collapse = ",\n")))
        } else {
            lines <- c(lines, "", "plot(result)")
        }
    }

    paste(lines, collapse = "\n")
}

#' Are these the same SingleCellExperiment?
#'
#' Cheap structural identity, avoiding a deep comparison of large assays.
is_same_sce <- function(a, b) {
    if (!inherits(a, "SummarizedExperiment") || !inherits(b, "SummarizedExperiment")) return(FALSE)
    identical(dim(a), dim(b)) && identical(colnames(a), colnames(b))
}
