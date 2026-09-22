# Declarative parameter objects.
#
# A parameter knows three things: how to render itself as a Shiny input, how to
# read its value back, and (when its choices depend on the data) how to refresh
# itself. Registry entries therefore describe their controls instead of writing
# them, which is what lets one generic module drive every diagnostic.
#
# `advanced = TRUE` moves a control into a collapsed section. Most package
# functions have a dozen arguments; showing all of them is the difference
# between a tool and a form. Defaults are the package's own defaults, so a user
# who touches nothing gets exactly what the vignettes show.

new_param <- function(subclass, id, label, default = NULL, help = NULL,
                      advanced = FALSE, ...) {
    structure(
        list(id = id, label = label, default = default, help = help,
             advanced = advanced, ...),
        class = c(paste0("sc_param_", subclass), "sc_param")
    )
}

param_int <- function(id, label, value, min = 1, max = 1000, step = 1,
                      help = NULL, advanced = FALSE) {
    new_param("int", id, label, value, help, advanced, min = min, max = max, step = step)
}

param_num <- function(id, label, value, min = 0, max = 1, step = 0.01,
                      help = NULL, advanced = FALSE) {
    new_param("num", id, label, value, help, advanced, min = min, max = max, step = step)
}

param_choice <- function(id, label, choices, selected = choices[1], multiple = FALSE,
                         help = NULL, advanced = FALSE) {
    new_param("choice", id, label, selected, help, advanced,
              choices = choices, multiple = multiple)
}

param_bool <- function(id, label, value = FALSE, help = NULL, advanced = FALSE) {
    new_param("bool", id, label, value, help, advanced)
}

#' Which principal components to use
param_pcs <- function(id = "pc_subset", label = "Principal components",
                      value = 1:5, max = 20, help = NULL, advanced = FALSE) {
    new_param("pcs", id, label, value, help, advanced, max = max)
}

#' Cell types, populated from the current data context
param_celltypes <- function(id = "cell_types", label = "Cell types",
                            multiple = TRUE, allow_empty = TRUE,
                            help = "Leave empty to use every shared cell type.",
                            advanced = FALSE) {
    new_param("celltypes", id, label, NULL, help, advanced,
              multiple = multiple, allow_empty = allow_empty)
}

#' A single gene, chosen server-side from the row names
param_gene <- function(id = "feature", label = "Gene", source = c("query", "ref"),
                       help = NULL, advanced = FALSE) {
    new_param("gene", id, label, NULL, help, advanced, source = match.arg(source))
}

#' A numeric colData column (QC metric or annotation score)
param_numeric_col <- function(id, label, source = c("query", "ref"),
                              prefer = character(0), help = NULL, advanced = FALSE) {
    new_param("numcol", id, label, NULL, help, advanced,
              source = match.arg(source), prefer = prefer)
}

#' A dimensionality reduction present on the object
param_reduction <- function(id = "method", label = "Reduction",
                            source = c("query", "ref"), help = NULL, advanced = FALSE) {
    new_param("reduction", id, label, NULL, help, advanced, source = match.arg(source))
}

#' Specific cells, chosen server-side from the column names
param_cellnames <- function(id = "cell_names", label = "Cells",
                            source = c("query", "ref"), max_select = 20,
                            help = NULL, advanced = FALSE) {
    new_param("cellnames", id, label, NULL, help, advanced,
              source = match.arg(source), max_select = max_select)
}

# --- rendering -------------------------------------------------------------

with_help <- function(input, p) {
    if (is.null(p$help)) return(input)
    htmltools::tagList(input, htmltools::div(class = "form-text sc-help", p$help))
}

# Dispatch on the parameter class.
#
# Written out rather than left to UseMethod: Shiny sources R/ into its own
# environment, not a package namespace, and implicit S3 dispatch does not find
# methods there when the generic is invoked indirectly - through lapply, for
# instance, which is exactly how render_params calls this. Explicit lookup in
# the environment the methods were defined in always works.
param_method <- function(generic, p) {
    env <- environment(param_method)
    for (cl in class(p)) {
        nm <- paste0(generic, ".", cl)
        f <- get0(nm, envir = env, mode = "function", inherits = TRUE)
        if (!is.null(f)) return(f)
    }
    stop("No ", generic, " method for a parameter of class ", class(p)[1], call. = FALSE)
}

param_input <- function(p, ns, ctx = NULL) param_method("param_input", p)(p, ns, ctx)

param_input.sc_param_int <- function(p, ns, ctx = NULL) {
    with_help(shiny::numericInput(ns(p$id), p$label, value = p$default,
                                  min = p$min, max = p$max, step = p$step), p)
}

param_input.sc_param_num <- function(p, ns, ctx = NULL) {
    with_help(shiny::numericInput(ns(p$id), p$label, value = p$default,
                                  min = p$min, max = p$max, step = p$step), p)
}

param_input.sc_param_choice <- function(p, ns, ctx = NULL) {
    with_help(shiny::selectInput(ns(p$id), p$label, choices = p$choices,
                                 selected = p$default, multiple = p$multiple), p)
}

param_input.sc_param_bool <- function(p, ns, ctx = NULL) {
    with_help(shiny::checkboxInput(ns(p$id), p$label, value = p$default), p)
}

param_input.sc_param_pcs <- function(p, ns, ctx = NULL) {
    avail <- ctx$n_pcs %||% p$max
    hi <- max(2, min(avail, p$max))
    val <- range(as_pc_subset(p$default))
    val[2] <- min(val[2], hi)
    with_help(
        shiny::sliderInput(ns(p$id), p$label, min = 1, max = hi,
                           value = val, step = 1, ticks = FALSE),
        p
    )
}

param_input.sc_param_celltypes <- function(p, ns, ctx = NULL) {
    with_help(shiny::selectizeInput(
        ns(p$id), p$label, choices = ctx$cell_types_available %||% character(0),
        selected = NULL, multiple = p$multiple,
        options = list(placeholder = if (p$allow_empty) "All shared cell types" else "Select…")
    ), p)
}

param_input.sc_param_gene <- function(p, ns, ctx = NULL) {
    with_help(shiny::selectizeInput(
        ns(p$id), p$label, choices = NULL, multiple = FALSE,
        options = list(placeholder = "Type to search genes…", maxOptions = 200)
    ), p)
}

param_input.sc_param_cellnames <- function(p, ns, ctx = NULL) {
    with_help(shiny::selectizeInput(
        ns(p$id), p$label, choices = NULL, multiple = TRUE,
        options = list(placeholder = "Type to search cells…", maxOptions = 200,
                       maxItems = p$max_select)
    ), p)
}

param_input.sc_param_numcol <- function(p, ns, ctx = NULL) {
    # Filled at render time as well as on refresh: the panel is built lazily,
    # so an empty control here would be visible until the next flush.
    obj <- if (identical(p$source, "ref")) ctx$ref else ctx$query
    cols <- if (is.null(obj)) character(0) else numeric_cols(obj)
    sel <- c(intersect(p$prefer, cols), cols)[1]
    with_help(shiny::selectInput(ns(p$id), p$label, choices = cols, selected = sel), p)
}

param_input.sc_param_reduction <- function(p, ns, ctx = NULL) {
    obj <- if (identical(p$source, "ref")) ctx$ref else ctx$query
    rd <- if (is.null(obj)) character(0) else SingleCellExperiment::reducedDimNames(obj)
    rd <- intersect(c("TSNE", "UMAP", "PCA"), rd)
    with_help(shiny::selectInput(ns(p$id), p$label, choices = rd, selected = rd[1]), p)
}

# --- refreshing data-dependent choices -------------------------------------

param_refresh <- function(p, session, ctx) param_method("param_refresh", p)(p, session, ctx)

param_refresh.sc_param <- function(p, session, ctx) invisible(NULL)

param_refresh.sc_param_pcs <- function(p, session, ctx) {
    hi <- max(2, min(ctx$n_pcs %||% p$max, p$max))
    cur <- shiny::isolate(session$input[[p$id]])
    val <- if (is.null(cur)) range(as_pc_subset(p$default)) else c(cur[1], min(cur[2], hi))
    shiny::updateSliderInput(session, p$id, max = hi, value = val)
}

param_refresh.sc_param_celltypes <- function(p, session, ctx) {
    choices <- ctx$cell_types_available %||% character(0)
    keep <- intersect(shiny::isolate(session$input[[p$id]]) %||% character(0), choices)
    # Where a cell type is required rather than optional, pick one rather than
    # leaving the control empty: an empty required control means the panel
    # opens on an error message, which reads as a broken app.
    if (!length(keep) && !isTRUE(p$allow_empty) && length(choices)) {
        keep <- (ctx$shared_types %||% choices)[1]
    }
    shiny::updateSelectizeInput(session, p$id, choices = choices, selected = keep,
                                server = length(choices) > 200)
}

param_refresh.sc_param_gene <- function(p, session, ctx) {
    obj <- if (p$source == "ref") ctx$ref else ctx$query
    genes <- if (is.null(obj)) character(0) else rownames(obj)
    cur <- shiny::isolate(session$input[[p$id]])
    sel <- if (!is.null(cur) && nzchar(cur) && cur %in% genes) cur else genes[1]
    # Server-side selectize: a 30,000-gene client-side list locks the browser.
    shiny::updateSelectizeInput(session, p$id, choices = genes, selected = sel, server = TRUE)
}

param_refresh.sc_param_cellnames <- function(p, session, ctx) {
    obj <- if (p$source == "ref") ctx$ref else ctx$query
    cells <- if (is.null(obj)) character(0) else colnames(obj)
    cur <- intersect(shiny::isolate(session$input[[p$id]]) %||% character(0), cells)
    shiny::updateSelectizeInput(session, p$id, choices = cells, selected = cur, server = TRUE)
}

param_refresh.sc_param_numcol <- function(p, session, ctx) {
    obj <- if (p$source == "ref") ctx$ref else ctx$query
    cols <- if (is.null(obj)) character(0) else numeric_cols(obj)
    preferred <- intersect(p$prefer, cols)
    cur <- shiny::isolate(session$input[[p$id]])
    sel <- if (!is.null(cur) && cur %in% cols) cur else c(preferred, cols)[1]
    shiny::updateSelectInput(session, p$id, choices = cols, selected = sel)
}

param_refresh.sc_param_reduction <- function(p, session, ctx) {
    obj <- if (p$source == "ref") ctx$ref else ctx$query
    rd <- if (is.null(obj)) character(0) else SingleCellExperiment::reducedDimNames(obj)
    # The package accepts these three; offer only what the object actually has.
    rd <- intersect(c("TSNE", "UMAP", "PCA"), rd)
    cur <- shiny::isolate(session$input[[p$id]])
    sel <- if (!is.null(cur) && cur %in% rd) cur else rd[1]
    shiny::updateSelectInput(session, p$id, choices = rd, selected = sel)
}

# --- reading values back ---------------------------------------------------

param_value <- function(p, input) param_method("param_value", p)(p, input)

param_value.sc_param <- function(p, input) input[[p$id]] %||% p$default

param_value.sc_param_int <- function(p, input) {
    v <- suppressWarnings(as.integer(input[[p$id]]))
    if (length(v) != 1L || is.na(v)) return(as.integer(p$default))
    # Clamping with max()/min() against numeric bounds would silently return a
    # double, and several package arguments are documented as integers.
    as.integer(max(p$min, min(p$max, v)))
}

param_value.sc_param_num <- function(p, input) {
    v <- suppressWarnings(as.numeric(input[[p$id]]))
    if (length(v) != 1L || is.na(v)) p$default else max(p$min, min(p$max, v))
}

param_value.sc_param_pcs <- function(p, input) {
    v <- input[[p$id]]
    if (is.null(v) || length(v) != 2L) return(as_pc_subset(p$default))
    as_pc_subset(seq(v[1], v[2]))
}

param_value.sc_param_celltypes <- function(p, input) {
    v <- input[[p$id]]
    if (is.null(v) || !length(v)) NULL else as.character(v)
}

param_value.sc_param_gene <- function(p, input) {
    v <- input[[p$id]]
    if (is.null(v) || !nzchar(v)) NULL else as.character(v)
}

param_value.sc_param_cellnames <- function(p, input) {
    v <- input[[p$id]]
    if (is.null(v) || !length(v)) NULL else as.character(v)
}

#' Collect every parameter's value into a named list
collect_params <- function(params, input) {
    if (!length(params)) return(list())
    out <- lapply(params, param_value, input = input)
    names(out) <- vapply(params, `[[`, character(1), "id")
    out
}

#' Render a parameter list, splitting basic from advanced
render_params <- function(params, ns, ctx) {
    if (!length(params)) return(NULL)
    basic <- Filter(function(p) !isTRUE(p$advanced), params)
    adv <- Filter(function(p) isTRUE(p$advanced), params)
    htmltools::tagList(
        lapply(basic, param_input, ns = ns, ctx = ctx),
        if (length(adv)) {
            bslib::accordion(
                open = FALSE, class = "sc-advanced",
                bslib::accordion_panel(
                    "Advanced options", icon = bsicons::bs_icon("sliders"),
                    lapply(adv, param_input, ns = ns, ctx = ctx)
                )
            )
        }
    )
}
