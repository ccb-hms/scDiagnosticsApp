# The generic diagnostic module.
#
# Given a registry entry and a data context, this builds the whole panel:
# controls, headline figures, plot, tables, interpretation, reproducible code
# and downloads. Every diagnostic in the app runs through it.

FOCUS_ID <- "__focus"

diagnostic_ui <- function(id, entry) {
    ns <- shiny::NS(id)
    bslib::layout_sidebar(
        sidebar = bslib::sidebar(
            width = 320, title = "Settings", open = "desktop",
            shiny::uiOutput(ns("controls")),
            shiny::uiOutput(ns("run_control")),
            shiny::div(
                class = "sc-sidebar-foot",
                shiny::downloadButton(ns("dl_plot"), "Plot (PNG)",
                                      class = "btn-sm btn-outline-secondary"),
                shiny::downloadButton(ns("dl_data"), "Table (CSV)",
                                      class = "btn-sm btn-outline-secondary")
            )
        ),
        shiny::div(
            class = "sc-panel",
            shiny::div(
                class = "sc-panel-head",
                shiny::h2(entry$label, class = "sc-panel-title"),
                if (!is.null(entry$question)) {
                    shiny::p(class = "sc-panel-question", entry$question)
                },
                shiny::p(class = "sc-panel-blurb", entry$blurb),
                shiny::div(
                    class = "sc-panel-meta",
                    shiny::tags$code(paste0(entry$fn, "()")),
                    shiny::tags$span(class = "sc-badge", CATEGORIES[[entry$category]]$label),
                    if (entry$heavy) shiny::tags$span(class = "sc-badge sc-badge-slow", "slower")
                )
            ),
            shiny::uiOutput(ns("gate")),
            shiny::uiOutput(ns("figures")),
            shiny::uiOutput(ns("messages")),
            bslib::navset_card_underline(
                id = ns("views"),
                bslib::nav_panel(
                    "Plot",
                    shiny::div(
                        class = "sc-plot-wrap",
                        shiny::plotOutput(ns("plot"), height = paste0(entry$plot_height, "px"))
                    ),
                    shiny::uiOutput(ns("plot_controls"))
                ),
                bslib::nav_panel("Numbers", shiny::uiOutput(ns("tables"))),
                bslib::nav_panel(
                    "How to read this",
                    shiny::div(class = "sc-reading",
                               shiny::p(entry$reading %||% entry$blurb),
                               shiny::uiOutput(ns("docs_link")))
                ),
                bslib::nav_panel(
                    "R code",
                    shiny::div(
                        class = "sc-code-wrap",
                        shiny::p(class = "text-muted small",
                                 "This reproduces the panel outside the app."),
                        shiny::verbatimTextOutput(ns("code")),
                        shiny::downloadButton(ns("dl_code"), "Download script",
                                              class = "btn-sm btn-outline-secondary")
                    )
                )
            )
        )
    )
}

diagnostic_server <- function(id, entry, ctx_r, session_cache, autorun = TRUE,
                              active = shiny::reactive(TRUE)) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        async_on <- isTRUE(entry$heavy) && isTRUE(getOption("sc.async.ready", FALSE))

        # --- controls ------------------------------------------------------
        all_params <- entry$params
        if (uses_focus(entry)) {
            all_params <- c(
                list(param_choice(FOCUS_ID, "Dataset",
                                  c(Query = "query", Reference = "ref"), "query",
                                  help = "Which object to run this on.")),
                all_params
            )
        }

        output$controls <- shiny::renderUI({
            ctx <- ctx_r()
            render_params(all_params, ns, ctx)
        })

        output$plot_controls <- shiny::renderUI({
            if (!length(entry$plot_params)) return(NULL)
            ctx <- ctx_r()
            shiny::div(class = "sc-plot-controls",
                       shiny::div(class = "sc-plot-controls-title", "Display"),
                       render_params(entry$plot_params, ns, ctx))
        })

        refresh_all <- function(ctx) {
            for (p in c(all_params, entry$plot_params)) param_refresh(p, session, ctx)
        }

        # Refresh data-dependent choices whenever the context changes.
        shiny::observeEvent(ctx_r(), refresh_all(ctx_r()), ignoreNULL = FALSE)

        # The panel's UI is built lazily when it is opened, so the controls do
        # not exist on the client until that flush completes. Server-side
        # selectize inputs (genes, cell names) must be populated after it,
        # otherwise the update message arrives before the element does.
        shiny::observeEvent(active(), {
            shiny::req(active())
            ctx <- shiny::isolate(ctx_r())
            session$onFlushed(function() shiny::isolate(refresh_all(ctx)), once = TRUE)
        }, ignoreInit = FALSE)

        # The effective context for this panel: the global one, with the
        # panel's own focus choice applied.
        panel_ctx <- shiny::reactive({
            ctx <- ctx_r()
            if (uses_focus(entry)) {
                which <- input[[FOCUS_ID]] %||% "query"
                ctx <- build_context(
                    ref = ctx$ref, query = ctx$query,
                    ref_desc = ctx$ref_desc, query_desc = ctx$query_desc,
                    ref_col = ctx$ref_col, query_col = ctx$query_col,
                    cell_types = ctx$cell_types, assay = ctx$assay,
                    pc_subset = ctx$pc_subset, focus_which = which,
                    preset_id = ctx$preset_id
                )
            }
            ctx$entry_label <- entry$label
            ctx
        })

        # --- gating --------------------------------------------------------
        gate_message <- shiny::reactive({
            ctx <- panel_ctx()
            if (!isTRUE(ctx$ready)) return("Choose a dataset on the Data page to begin.")
            entry_ready(entry, ctx)
        })

        output$gate <- shiny::renderUI({
            msg <- gate_message()
            if (is.null(msg)) return(NULL)
            alert(msg, type = "info")
        })

        # --- resolving the call --------------------------------------------
        resolved <- shiny::reactive({
            shiny::req(is.null(gate_message()))
            ctx <- panel_ctx()
            pv <- collect_params(entry$params, input)
            args <- tryCatch(build_args(entry, ctx, pv),
                             error = function(e) sc_failure(conditionMessage(e)))
            list(ctx = ctx, params = pv, args = args)
        })

        # --- running --------------------------------------------------------
        # Light diagnostics recompute as soon as a control changes. Heavy ones
        # wait for an explicit click, because silently starting a 20-second job
        # on every keystroke is how the previous app felt slow.
        run_trigger <- shiny::reactiveVal(0)

        output$run_control <- shiny::renderUI({
            if (!entry$heavy) {
                return(shiny::div(class = "sc-autorun small text-muted",
                                  bsicons::bs_icon("lightning-charge"),
                                  " Updates automatically"))
            }
            shiny::tagList(
                bslib::input_task_button(ns("run"), "Run diagnostic",
                                         icon = shiny::icon("play"),
                                         class = "btn-primary w-100"),
                if (async_on) {
                    shiny::div(class = "sc-autorun small text-muted mt-1",
                               bsicons::bs_icon("cpu"), " Runs in the background")
                }
            )
        })

        # The button is recreated whenever the panel's markup is rebuilt, which
        # resets its count to zero. Only a real press should start a run.
        shiny::observeEvent(input$run, {
            if (isTRUE(input$run > 0)) shiny::isolate(run_trigger(run_trigger() + 1))
        }, ignoreInit = TRUE)

        # Heavy panels run once when first opened so the user sees something
        # immediately, then wait for the button after that.
        first_run_done <- shiny::reactiveVal(FALSE)
        shiny::observe({
            shiny::req(autorun, active(), entry$heavy, !first_run_done(),
                       is.null(gate_message()))
            first_run_done(TRUE)
            shiny::isolate(run_trigger(run_trigger() + 1))
        })

        # Changing the data invalidates a heavy panel's result. Keyed on the
        # data signature, so moving a slider does not throw away a result the
        # user then has to recompute by hand.
        shiny::observeEvent(context_signature(ctx_r()), {
            if (entry$heavy) {
                first_run_done(FALSE)
                result_store(NULL)
            }
        }, ignoreInit = TRUE)

        result_store <- shiny::reactiveVal(NULL)
        pending <- shiny::reactiveVal(NULL)

        # A heavy diagnostic runs in a worker process when one is available, so
        # it does not block other sessions sharing this R process.
        task <- if (entry$heavy) new_diagnostic_task() else NULL
        if (!is.null(task)) bslib::bind_task_button(task, "run")

        cache_hit <- function(res) {
            tier <- cache_for(res$ctx, session_cache)
            key <- cache_key(entry$id, res$ctx, res$params)
            hit <- tier$cache$get(key)
            if (cachem::is.key_missing(hit)) NULL else hit
        }

        cache_store <- function(res, value) {
            if (inherits(value, "sc_failure")) return(invisible(NULL))
            tier <- cache_for(res$ctx, session_cache)
            key <- cache_key(entry$id, res$ctx, res$params)
            tryCatch(tier$cache$set(key, value), error = function(e) NULL)
            invisible(NULL)
        }

        shiny::observe({
            # Only the visible panel computes; the rest stay idle.
            shiny::req(active())
            if (entry$heavy) {
                shiny::req(run_trigger() > 0)
                res <- shiny::isolate(resolved())
            } else {
                res <- resolved()
            }
            shiny::req(res)

            if (is_failure(res$args)) {
                result_store(list(value = res$args, ctx = res$ctx,
                                  params = res$params, args = res$args))
                return(invisible(NULL))
            }

            hit <- cache_hit(res)
            if (!is.null(hit)) {
                attr(hit, "sc_cache_hit") <- TRUE
                result_store(list(value = hit, ctx = res$ctx,
                                  params = res$params, args = res$args))
                return(invisible(NULL))
            }

            if (!is.null(task) && use_async(entry, res$ctx)) {
                pending(res)
                task$invoke(APP_DIR, run_spec(entry, res$ctx, res$params))
                return(invisible(NULL))
            }

            prog <- shiny::Progress$new(session, min = 0, max = 1)
            prog$set(message = paste("Running", entry$label), value = 0.35)
            on.exit(prog$close(), add = TRUE)
            out <- run_entry(entry, res$args)
            cache_store(res, out)
            result_store(list(value = out, ctx = res$ctx, params = res$params, args = res$args))
        })

        if (!is.null(task)) {
            shiny::observe({
                st <- task$status()
                shiny::req(st %in% c("success", "error"))
                res <- shiny::isolate(pending())
                shiny::req(res)
                out <- if (identical(st, "success")) {
                    tryCatch(task$result(),
                             error = function(e) sc_failure(conditionMessage(e)))
                } else {
                    # A worker that dies should not take the panel with it.
                    val <- tryCatch({
                        task$result()
                        NULL
                    }, error = function(e) sc_failure(conditionMessage(e)))
                    val %||% sc_failure("The background job failed without a message.")
                }
                cache_store(res, out)
                result_store(list(value = out, ctx = res$ctx,
                                  params = res$params, args = res$args))
                pending(NULL)
            })
        }

        # --- plot-method arguments -----------------------------------------
        # Cell type pickers on the plot side are filled from the result, not
        # from the context, because a diagnostic may have dropped cell types
        # that were too small to analyse.
        shiny::observeEvent(result_store(), {
            st <- result_store()
            if (is.null(st) || is_failure(st$value)) return()
            names_in_result <- names(st$value)
            for (p in entry$plot_params) {
                if (p$id %in% c("cell_type", "ref_cell_type") &&
                    inherits(p, "sc_param_choice")) {
                    ch <- names_in_result
                    if (!length(ch)) next
                    sel <- if (!is.null(input[[p$id]]) && input[[p$id]] %in% ch) {
                        input[[p$id]]
                    } else if ("Combined" %in% ch) "Combined" else ch[1]
                    shiny::updateSelectInput(session, p$id, choices = ch, selected = sel)
                }
            }
        })

        plot_args <- shiny::reactive({
            if (!length(entry$plot_params)) return(list())
            pv <- collect_params(entry$plot_params, input)
            # A plot-side cell type picker named plot_cell_types feeds the
            # method's `cell_types` argument.
            if (!is.null(pv$plot_cell_types)) {
                pv$cell_types <- pv$plot_cell_types
                pv$plot_cell_types <- NULL
            }
            Filter(function(v) !is.null(v), pv)
        })

        # --- outputs --------------------------------------------------------
        output$figures <- shiny::renderUI({
            st <- result_store()
            if (is.null(st) || is_failure(st$value) || is.null(entry$summarise)) return(NULL)
            figs <- tryCatch(entry$summarise(st$value), error = function(e) list())
            if (!length(figs)) return(NULL)
            bslib::layout_column_wrap(
                width = "220px", fixed_width = FALSE, class = "sc-figures",
                !!!lapply(figs, function(f) {
                    shiny::div(
                        class = paste("sc-figure", tone_class(f$tone)),
                        shiny::div(class = "sc-figure-label", f$label),
                        shiny::div(class = "sc-figure-value", f$value),
                        if (!is.null(f$caption)) {
                            shiny::div(class = "sc-figure-caption", f$caption)
                        }
                    )
                })
            )
        })

        output$messages <- shiny::renderUI({
            st <- result_store()
            if (is.null(st)) return(NULL)
            if (is_failure(st$value)) {
                return(error_panel(st$value$message, explain_failure(entry, st)))
            }
            w <- attr(st$value, "sc_warnings")
            if (!length(w)) return(NULL)
            alert(type = "warning",
                  shiny::tags$strong("Notes from the computation"),
                  shiny::tags$ul(lapply(unique(w), shiny::tags$li)))
        })

        current_plot <- shiny::reactive({
            st <- result_store()
            shiny::req(st, !is_failure(st$value))
            build_plot(entry, st$value, plot_args(), st$ctx)
        })

        output$plot <- shiny::renderPlot({
            st <- result_store()
            if (is.null(st)) {
                # A background job leaves the panel empty for a while, so say
                # what is happening rather than inviting a second click.
                if (!is.null(task) && identical(task$status(), "running")) {
                    return(blank_note("Running in the background…"))
                }
                return(blank_note(if (entry$heavy) "Press Run to compute this diagnostic."
                                  else "Waiting for data…"))
            }
            if (is_failure(st$value)) return(blank_note("See the message above."))
            p <- tryCatch(current_plot(), error = function(e) e)
            if (inherits(p, "error")) return(blank_note(conditionMessage(p)))
            draw_result(p)
        }, res = 100)

        result_tables <- shiny::reactive({
            st <- result_store()
            shiny::req(st, !is_failure(st$value))
            tabs <- if (!is.null(entry$tables)) {
                tryCatch(entry$tables(st$value), error = function(e) NULL)
            } else NULL
            if (is.null(tabs) || !length(tabs)) tabs <- auto_tables(st$value)
            Filter(function(d) !is.null(d) && NROW(d) > 0, tabs)
        })

        output$tables <- shiny::renderUI({
            st <- result_store()
            if (is.null(st)) return(empty_state("Nothing computed yet."))
            if (is_failure(st$value)) return(empty_state("See the message above.", "exclamation-circle"))
            tabs <- result_tables()
            if (!length(tabs)) return(empty_state("This diagnostic returns a figure only.", "bar-chart"))
            panels <- Map(function(nm, df) {
                bslib::nav_panel(nm, DT::DTOutput(ns(table_out_id(nm))))
            }, names(tabs), tabs)
            do.call(bslib::navset_pill, c(list(id = ns("tabletabs")), unname(panels)))
        })

        # Table outputs are created lazily and by name, so any number of
        # tables can be rendered without declaring them up front.
        shiny::observe({
            tabs <- result_tables()
            for (nm in names(tabs)) {
                local({
                    this <- nm
                    output[[table_out_id(this)]] <- DT::renderDT({
                        df <- result_tables()[[this]]
                        shiny::req(df)
                        nice_table(df)
                    })
                })
            }
        })

        output$docs_link <- shiny::renderUI({
            shiny::tagList(
                shiny::p(class = "small text-muted mt-3",
                         "Full argument documentation for this function:"),
                shiny::a(
                    href = sprintf("https://ccb-hms.github.io/scDiagnostics/reference/%s.html", entry$fn),
                    target = "_blank", rel = "noopener",
                    class = "btn btn-sm btn-outline-primary",
                    bsicons::bs_icon("book"), sprintf(" ?%s", entry$fn)
                )
            )
        })

        code_text <- shiny::reactive({
            st <- result_store()
            if (is.null(st) || is_failure(st$args)) {
                res <- tryCatch(resolved(), error = function(e) NULL)
                if (is.null(res) || is_failure(res$args)) return("# Configure the panel to see the code.")
                return(repro_code(entry, res$ctx, res$args, plot_args(), res$params))
            }
            repro_code(entry, st$ctx, st$args, plot_args(), st$params)
        })

        output$code <- shiny::renderText(code_text())

        # --- downloads ------------------------------------------------------
        output$dl_plot <- shiny::downloadHandler(
            filename = function() sprintf("%s_%s.png", entry$id, Sys.Date()),
            content = function(file) {
                grDevices::png(file, width = 2400, height = 1800, res = 220)
                on.exit(grDevices::dev.off(), add = TRUE)
                p <- tryCatch(current_plot(), error = function(e) NULL)
                if (is.null(p)) blank_note("No plot available.") else draw_result(p)
            }
        )

        # The sidebar button exports the primary table. Every other table has
        # its own CSV button in the Numbers view, so nothing is unreachable.
        output$dl_data <- shiny::downloadHandler(
            filename = function() {
                tabs <- tryCatch(result_tables(), error = function(e) list())
                slug <- if (length(tabs)) gsub("[^A-Za-z0-9]+", "-", tolower(names(tabs)[1])) else "result"
                sprintf("%s_%s_%s.csv", entry$id, slug, Sys.Date())
            },
            content = function(file) {
                tabs <- tryCatch(result_tables(), error = function(e) list())
                if (!length(tabs)) {
                    utils::write.csv(data.frame(note = "This diagnostic returns a figure only."),
                                     file, row.names = FALSE)
                    return(invisible(NULL))
                }
                utils::write.csv(tabs[[1]], file, row.names = FALSE)
            }
        )

        output$dl_code <- shiny::downloadHandler(
            filename = function() sprintf("%s_%s.R", entry$id, Sys.Date()),
            content = function(file) writeLines(code_text(), file)
        )

        # No always_render here. The panel's markup is inserted only while the
        # Diagnostics tab is on screen, so these outputs are visible from the
        # moment they are bound and Shiny's normal suspension behaviour is
        # exactly what we want for the panels the user is not looking at.

        # Expose the panel's state so the Audit and Report pages can reuse it.
        shiny::reactive(result_store())
    })
}

# --- rendering helpers -----------------------------------------------------

table_out_id <- function(nm) paste0("tbl_", gsub("[^A-Za-z0-9]+", "_", nm))

#' Build the plot for a result, whatever form the entry returns
build_plot <- function(entry, value, plot_args, ctx) {
    if (!is.null(entry$plot_fn)) return(entry$plot_fn(value, plot_args, ctx))
    if (entry$returns == "plot") return(value)
    if (entry$returns %in% c("value", "sce")) return(NULL)
    args <- c(list(value), plot_args)
    # Drop arguments the method does not accept, so a shared control (for
    # example a cell type picker) can be reused across methods safely.
    m <- utils::getS3method("plot", class(value)[length(class(value))], optional = TRUE)
    if (!is.null(m)) {
        ok <- names(formals(m))
        if (!"..." %in% ok) args <- args[names(args) == "" | names(args) %in% ok]
    }
    do.call(graphics::plot, args)
}

#' Draw whatever a plot builder returned
#'
#' Package plots come back as ggplot objects, GGally matrices, ComplexHeatmap
#' objects, or as NULL after drawing straight to the device with base graphics.
draw_result <- function(p) {
    if (is.null(p)) return(invisible(NULL))
    if (inherits(p, c("gg", "ggplot", "ggmatrix", "patchwork"))) return(print(p))
    if (inherits(p, c("Heatmap", "HeatmapList"))) return(ComplexHeatmap::draw(p))
    if (inherits(p, "grob") || inherits(p, "gtable")) {
        grid::grid.newpage()
        return(grid::grid.draw(p))
    }
    if (inherits(p, "recordedplot")) return(grDevices::replayPlot(p))
    invisible(NULL)
}

#' A blank plot carrying a short message, used instead of an empty device
blank_note <- function(msg) {
    ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = msg,
                          size = 4.2, colour = SC_COLORS$secondary) +
        ggplot2::theme_void()
}

#' Consistent DT configuration
nice_table <- function(df) {
    num <- vapply(df, is.numeric, logical(1))
    dt <- DT::datatable(
        df,
        rownames = FALSE,
        extensions = "Buttons",
        class = "compact stripe hover",
        options = list(
            pageLength = 15,
            scrollX = TRUE,
            dom = "Bfrtip",
            buttons = list(list(extend = "csv", text = "CSV"),
                           list(extend = "copy", text = "Copy")),
            columnDefs = list(list(className = "dt-right", targets = which(num) - 1))
        )
    )
    if (any(num)) {
        cols <- names(df)[num]
        # Integers stay integers; everything else gets a sane number of digits.
        is_int <- vapply(df[cols], function(v) all(is.na(v) | v == round(v)), logical(1))
        if (any(!is_int)) dt <- DT::formatSignif(dt, cols[!is_int], digits = 4)
    }
    dt
}

#' Turn a raw error into something the user can act on
explain_failure <- function(entry, st) {
    msg <- st$value$message
    hints <- c(
        "only one cell type" = paste(
            "Select a single cell type in the settings: this diagnostic",
            "compares one population at a time."),
        "not found in" = paste(
            "The chosen annotation column is missing from one of the datasets.",
            "Set the columns on the Data page."),
        "pre-computed PCA" = paste(
            "This object has no usable PCA. Open the Data page and let the app",
            "compute one."),
        "Not enough" = paste(
            "Too few cells in at least one cell type. Try selecting fewer cell",
            "types, or lower the minimum cell thresholds."),
        "subscript out of bounds" = paste(
            "A cell type or gene in the settings is absent from the data.",
            "Reset the selection and try again.")
    )
    hit <- names(hints)[vapply(names(hints), function(k) grepl(k, msg, fixed = TRUE), logical(1))]
    if (length(hit)) return(hints[[hit[1]]])
    sprintf("Raised by %s(). Check the settings in the sidebar, or see the R code tab for the exact call.", entry$fn)
}
