# Data selection: presets, uploads, validation and the shared context.
#
# Everything downstream reads the context this module returns. Choosing data is
# therefore a single act, not something repeated on every panel.

MAX_CELLS_DEFAULT <- as.integer(Sys.getenv("SCDIAG_MAX_CELLS", unset = "20000"))

data_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::div(
        class = "sc-page",
        shiny::div(
            class = "sc-page-head",
            shiny::h2("Data", class = "sc-page-title"),
            shiny::p(
                class = "sc-page-lede",
                "Pick a prepared dataset or upload your own. Everything else in the app runs on what you choose here."
            )
        ),
        bslib::navset_card_underline(
            id = ns("source"),
            bslib::nav_panel(
                title = shiny::tagList(bsicons::bs_icon("collection"), " Prepared datasets"),
                value = "preset",
                shiny::p(class = "text-muted",
                         "These ship with the scDiagnostics package, so they always match the installed version."),
                shiny::uiOutput(ns("preset_cards"))
            ),
            bslib::nav_panel(
                title = shiny::tagList(bsicons::bs_icon("upload"), " Your own data"),
                value = "upload",
                shiny::div(
                    class = "row g-3",
                    shiny::div(
                        class = "col-md-6",
                        bslib::card(
                            bslib::card_header("Reference"),
                            bslib::card_body(
                                shiny::p(class = "small text-muted",
                                         "The annotated dataset you trust."),
                                shiny::fileInput(ns("ref_file"), NULL,
                                                 accept = c(".rds", ".RDS", ".h5ad"),
                                                 buttonLabel = "Choose file",
                                                 placeholder = ".rds or .h5ad"),
                                shiny::uiOutput(ns("ref_report"))
                            )
                        )
                    ),
                    shiny::div(
                        class = "col-md-6",
                        bslib::card(
                            bslib::card_header("Query"),
                            bslib::card_body(
                                shiny::p(class = "small text-muted",
                                         "The dataset whose annotations you want to audit."),
                                shiny::fileInput(ns("query_file"), NULL,
                                                 accept = c(".rds", ".RDS", ".h5ad"),
                                                 buttonLabel = "Choose file",
                                                 placeholder = ".rds or .h5ad"),
                                shiny::uiOutput(ns("query_report"))
                            )
                        )
                    )
                ),
                shiny::uiOutput(ns("upload_help"))
            )
        ),
        shiny::uiOutput(ns("configure")),
        shiny::uiOutput(ns("summary"))
    )
}

data_server <- function(id, jump_to = function(...) NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        state <- shiny::reactiveValues(
            ref = NULL, query = NULL,
            ref_desc = NULL, query_desc = NULL,
            preset_id = NULL,
            ref_checks = NULL, query_checks = NULL,
            notes = character(0)
        )

        # --- presets --------------------------------------------------------
        output$preset_cards <- shiny::renderUI({
            cards <- lapply(PRESETS, function(p) {
                active <- identical(state$preset_id, p$id)
                shiny::div(
                    class = paste("sc-preset", if (active) "sc-preset-active"),
                    shiny::div(class = "sc-preset-title", p$label),
                    shiny::div(class = "sc-preset-tagline", p$tagline),
                    shiny::p(class = "sc-preset-desc", p$description),
                    shiny::div(
                        class = "sc-preset-foot",
                        shiny::actionButton(ns(paste0("use_", p$id)),
                                            if (active) "Selected" else "Use this",
                                            class = if (active) "btn-sm btn-success" else "btn-sm btn-primary"),
                        shiny::a("Case study", href = p$vignette, target = "_blank",
                                 rel = "noopener", class = "btn btn-sm btn-link")
                    )
                )
            })
            bslib::layout_column_wrap(width = "320px", fixed_width = FALSE,
                                      class = "sc-preset-grid", !!!unname(cards))
        })

        apply_preset <- function(pid, notify = TRUE) {
            p <- PRESETS[[pid]]
            if (is.null(p)) return(invisible(NULL))
            prog <- shiny::Progress$new(session, min = 0, max = 1)
            prog$set(message = paste("Loading", p$label), value = 0.3)
            on.exit(prog$close(), add = TRUE)

            state$ref_desc <- preset_descriptor(p$ref$data)
            state$query_desc <- preset_descriptor(p$query$data)
            state$ref <- materialize(state$ref_desc)
            state$query <- materialize(state$query_desc)
            state$preset_id <- pid
            state$ref_checks <- NULL
            state$query_checks <- NULL
            state$notes <- character(0)

            prog$set(value = 0.8)
            shiny::updateSelectInput(session, "ref_col",
                                     choices = candidate_celltype_cols(state$ref),
                                     selected = p$ref$col)
            shiny::updateSelectInput(session, "query_col",
                                     choices = candidate_celltype_cols(state$query),
                                     selected = p$query$col)
            shiny::updateSliderInput(session, "pcs",
                                     value = range(p$pc_subset %||% 1:5))
            if (notify) {
                shiny::showNotification(sprintf("Loaded %s.", p$label),
                                        type = "message", duration = 4)
            }
            invisible(NULL)
        }

        # The cards are re-rendered whenever the selection changes, which
        # resets every button to zero. Without the count guard, observeEvent
        # reads that reset as a click on all of them and the page loops.
        lapply(names(PRESETS), function(pid) {
            key <- paste0("use_", pid)
            shiny::observeEvent(input[[key]], {
                if (isTRUE(input[[key]] > 0)) apply_preset(pid)
            }, ignoreInit = TRUE)
        })

        # --- uploads --------------------------------------------------------
        read_upload <- function(path, name) {
            ext <- tolower(tools::file_ext(name))
            if (ext == "h5ad") {
                if (!requireNamespace("zellkonverter", quietly = TRUE)) {
                    stop(paste("Reading .h5ad needs the zellkonverter package, which is",
                               "not installed on this server. Save your object as .rds",
                               "instead."), call. = FALSE)
                }
                return(zellkonverter::readH5AD(path))
            }
            readRDS(path)
        }

        handle_upload <- function(which, file) {
            shiny::req(file)
            prog <- shiny::Progress$new(session, min = 0, max = 1)
            prog$set(message = paste("Reading", file$name), value = 0.15)
            on.exit(prog$close(), add = TRUE)

            obj <- tryCatch(read_upload(file$datapath, file$name),
                            error = function(e) e)
            if (inherits(obj, "error")) {
                state[[paste0(which, "_checks")]] <- list(list(
                    status = "fail", label = "Could not read the file",
                    detail = conditionMessage(obj)
                ))
                return(invisible(NULL))
            }

            other <- if (which == "ref") state$query else state$ref
            checks <- inspect_upload(obj, other, which)
            state[[paste0(which, "_checks")]] <- checks
            if (upload_blocked(checks)) return(invisible(NULL))

            notes <- character(0)

            # Bound memory before anything else touches the object.
            if (ncol(obj) > MAX_CELLS_DEFAULT) {
                prog$set(message = "Subsampling", value = 0.4)
                col <- candidate_celltype_cols(obj)[1]
                before <- ncol(obj)
                obj <- subsample_sce(obj, MAX_CELLS_DEFAULT, col)
                notes <- c(notes, sprintf(
                    paste("Subsampled from %s to %s cells, keeping cell type proportions.",
                          "Diagnostics are unchanged in kind; small populations are preserved."),
                    n_fmt(before), n_fmt(ncol(obj))
                ))
            }

            if (!has_valid_pca(obj)) {
                prog$set(message = "Computing PCA", value = 0.6)
                obj <- tryCatch(
                    scDiagnostics::processPCA(obj, n_hvgs = 2000),
                    error = function(e) e
                )
                if (inherits(obj, "error")) {
                    state[[paste0(which, "_checks")]] <- c(checks, list(list(
                        status = "fail", label = "PCA could not be computed",
                        detail = conditionMessage(obj)
                    )))
                    return(invisible(NULL))
                }
                notes <- c(notes, "Computed a PCA from the top 2,000 highly variable genes.")
            }

            # Persist the processed object so background workers can read it
            # from disk rather than receive a serialised copy, and so the cache
            # key can be a content hash.
            prog$set(message = "Preparing", value = 0.85)
            path <- file.path(tempdir(), sprintf("sc-%s-%s.rds", which, hash_of(file$name, Sys.time())))
            saveRDS(obj, path)

            state[[which]] <- obj
            state[[paste0(which, "_desc")]] <- file_descriptor(path, file$name)
            state$preset_id <- NULL
            state$notes <- c(state$notes, notes)

            cols <- candidate_celltype_cols(obj)
            shiny::updateSelectInput(session, paste0(which, "_col"),
                                     choices = cols, selected = cols[1])
            shiny::showNotification(sprintf("%s loaded: %s cells.",
                                            if (which == "ref") "Reference" else "Query",
                                            n_fmt(ncol(obj))),
                                    type = "message", duration = 4)
            invisible(NULL)
        }

        shiny::observeEvent(input$ref_file, handle_upload("ref", input$ref_file))
        shiny::observeEvent(input$query_file, handle_upload("query", input$query_file))

        render_checks <- function(checks) {
            if (is.null(checks)) return(NULL)
            shiny::div(
                class = "sc-checks",
                lapply(checks, function(c) {
                    icon <- switch(c$status, ok = "check-circle-fill",
                                   warn = "exclamation-triangle-fill", "x-circle-fill")
                    shiny::div(
                        class = paste0("sc-check sc-check-", c$status),
                        bsicons::bs_icon(icon),
                        shiny::div(
                            shiny::div(class = "sc-check-label", c$label),
                            shiny::div(class = "sc-check-detail", c$detail)
                        )
                    )
                })
            )
        }

        output$ref_report <- shiny::renderUI(render_checks(state$ref_checks))
        output$query_report <- shiny::renderUI(render_checks(state$query_checks))

        output$upload_help <- shiny::renderUI({
            bslib::accordion(
                open = FALSE, class = "mt-3",
                bslib::accordion_panel(
                    "Preparing your data",
                    icon = bsicons::bs_icon("question-circle"),
                    shiny::tags$ul(
                        shiny::tags$li(shiny::HTML(paste(
                            "A <code>SingleCellExperiment</code> or",
                            "<code>SpatialExperiment</code> saved with <code>saveRDS()</code>,",
                            "or an AnnData <code>.h5ad</code> file."))),
                        shiny::tags$li(shiny::HTML(paste(
                            "A <code>logcounts</code> assay of log-normalised values. Raw",
                            "counts alone will give misleading distances — run",
                            "<code>scuttle::logNormCounts()</code> first."))),
                        shiny::tags$li(shiny::HTML("A cell type annotation in <code>colData</code>.")),
                        shiny::tags$li("A PCA is helpful but optional; the app computes one if it is missing."),
                        shiny::tags$li(shiny::HTML(sprintf(
                            "Objects larger than %s cells are subsampled, keeping cell type proportions.",
                            n_fmt(MAX_CELLS_DEFAULT))))
                    ),
                    shiny::tags$pre(class = "sc-snippet", paste(
                        "library(scuttle)",
                        "sce <- logNormCounts(sce)",
                        'saveRDS(sce, "my_reference.rds")',
                        sep = "\n"
                    )),
                    shiny::p(class = "small text-muted mb-0",
                             "Uploaded data stays in this session's memory and is not shared with other users.")
                )
            )
        })

        # --- configuration ---------------------------------------------------
        output$configure <- shiny::renderUI({
            if (is.null(state$ref)) return(NULL)
            assays <- common_assays(state$ref, state$query)
            bslib::card(
                class = "mt-3",
                bslib::card_header(bsicons::bs_icon("sliders"), " Configure"),
                bslib::card_body(
                    shiny::div(
                        class = "row g-3",
                        shiny::div(class = "col-md-4",
                                   shiny::selectInput(ns("ref_col"), "Reference cell type column",
                                                      choices = candidate_celltype_cols(state$ref))),
                        shiny::div(class = "col-md-4",
                                   if (!is.null(state$query)) {
                                       shiny::selectInput(ns("query_col"), "Query cell type column",
                                                          choices = candidate_celltype_cols(state$query))
                                   }),
                        shiny::div(class = "col-md-4",
                                   shiny::selectInput(ns("assay"), "Assay", choices = assays,
                                                      selected = assays[1]))
                    ),
                    shiny::div(
                        class = "row g-3",
                        shiny::div(class = "col-md-6",
                                   shiny::sliderInput(ns("pcs"), "Principal components",
                                                      min = 1, max = max(2, pca_ncol(state$ref) %||% 20),
                                                      value = c(1, 5), step = 1, ticks = FALSE)),
                        shiny::div(class = "col-md-6",
                                   shiny::selectizeInput(
                                       ns("cell_types"), "Cell types",
                                       choices = NULL, multiple = TRUE,
                                       options = list(placeholder = "All shared cell types")
                                   ))
                    ),
                    shiny::div(
                        class = "sc-config-actions",
                        shiny::actionButton(ns("go_audit"), "Run the guided audit",
                                            class = "btn-primary",
                                            icon = shiny::icon("wand-magic-sparkles")),
                        shiny::actionButton(ns("go_diagnostics"), "Browse all diagnostics",
                                            class = "btn-outline-primary")
                    )
                )
            )
        })

        # Which cell types are on offer, derived from the data and the chosen
        # columns only. This must not depend on the current selection: the
        # observer below writes that selection back, and a dependency here
        # would make the two chase each other indefinitely.
        available_types <- shiny::reactive({
            ref_col <- input$ref_col %|% default_col(state$ref, state$preset_id, "ref")
            query_col <- input$query_col %|% default_col(state$query, state$preset_id, "query")
            types_of <- function(sce, col) {
                if (is.null(sce) || is.null(col) ||
                    !col %in% colnames(SummarizedExperiment::colData(sce))) return(character(0))
                sort(unique(stats::na.omit(as.character(sce[[col]]))))
            }
            rt <- types_of(state$ref, ref_col)
            qt <- types_of(state$query, query_col)
            shared <- intersect(rt, qt)
            unique(c(shared, setdiff(rt, shared), setdiff(qt, shared)))
        })

        shiny::observeEvent(available_types(), {
            choices <- available_types()
            shiny::updateSelectizeInput(
                session, "cell_types", choices = choices,
                selected = intersect(shiny::isolate(input$cell_types) %||% character(0), choices)
            )
        }, ignoreNULL = FALSE)

        shiny::observeEvent(input$go_audit, jump_to("audit"))
        shiny::observeEvent(input$go_diagnostics, jump_to("diagnostics"))

        # --- the context -----------------------------------------------------
        ctx <- shiny::reactive({
            build_context(
                ref = state$ref, query = state$query,
                ref_desc = state$ref_desc, query_desc = state$query_desc,
                ref_col = input$ref_col %|% default_col(state$ref, state$preset_id, "ref"),
                query_col = input$query_col %|% default_col(state$query, state$preset_id, "query"),
                cell_types = input$cell_types,
                assay = input$assay %|% "logcounts",
                pc_subset = if (is.null(input$pcs)) 1:5 else seq(input$pcs[1], input$pcs[2]),
                preset_id = state$preset_id
            )
        })

        output$summary <- shiny::renderUI({
            c <- ctx()
            if (!isTRUE(c$ready)) return(NULL)
            shiny::tagList(
                if (length(state$notes)) {
                    alert(type = "info",
                          shiny::tags$strong("What the app did to your data"),
                          shiny::tags$ul(lapply(unique(state$notes), shiny::tags$li)))
                },
                bslib::card(
                    class = "mt-3",
                    bslib::card_header(bsicons::bs_icon("info-circle"), " Current selection"),
                    bslib::card_body(shiny::div(class = "sc-summary", context_summary(c)))
                )
            )
        })

        always_render(output, c("preset_cards", "ref_report", "query_report",
                                "upload_help", "configure", "summary"))

        list(ctx = ctx, load_preset = apply_preset)
    })
}

#' Preset default column, used before the configure panel has rendered
default_col <- function(sce, preset_id, role) {
    if (is.null(sce)) return(NULL)
    if (!is.null(preset_id)) {
        p <- PRESETS[[preset_id]]
        col <- if (role == "ref") p$ref$col else p$query$col
        if (!is.null(col)) return(col)
    }
    candidate_celltype_cols(sce)[1]
}

#' Rich summary block for the chosen data
context_summary <- function(c) {
    row <- function(label, value) {
        shiny::div(class = "sc-summary-row",
                   shiny::div(class = "sc-summary-label", label),
                   shiny::div(class = "sc-summary-value", value))
    }
    only_ref <- setdiff(c$ref_types, c$shared_types)
    only_query <- setdiff(c$query_types, c$shared_types)
    shiny::tagList(
        row("Reference", sprintf("%s · %s genes × %s cells · %s",
                                 c$ref_desc$label %||% "–", n_fmt(nrow(c$ref)),
                                 n_fmt(ncol(c$ref)), c$ref_col)),
        if (!is.null(c$query)) {
            row("Query", sprintf("%s · %s genes × %s cells · %s",
                                 c$query_desc$label %||% "–", n_fmt(nrow(c$query)),
                                 n_fmt(ncol(c$query)), c$query_col))
        },
        row("Shared cell types", if (length(c$shared_types)) abbrev(c$shared_types, 10) else "none"),
        if (length(only_ref)) {
            row("Reference only",
                shiny::tagList(abbrev(only_ref, 8),
                               shiny::span(class = "sc-note",
                                           " — no query cells carry these labels")))
        },
        if (length(only_query)) {
            row("Query only",
                shiny::tagList(abbrev(only_query, 8),
                               shiny::span(class = "sc-note",
                                           " — absent from the reference, so nothing to compare against")))
        },
        row("Assay", c$assay),
        row("Components", sprintf("PC%d–PC%d of %d available",
                                  min(c$pc_subset), max(c$pc_subset), c$n_pcs))
    )
}
