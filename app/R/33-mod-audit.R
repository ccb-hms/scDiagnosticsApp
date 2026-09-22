# The guided audit.
#
# The package's three-step workflow - project, detect, characterize - run end
# to end with sensible defaults, summarised in plain language. This is the page
# a first-time visitor should land on: three clicks from opening the app to a
# defensible statement about their annotations.

AUDIT_STEPS <- list(
    list(id = "project", n = "1", title = "Project",
         entry = "plotCellTypePCA",
         question = "Does the query sit where the reference says it should?"),
    list(id = "detect", n = "2", title = "Detect",
         entry = "detectAnomaly",
         question = "Which cells do not look like their assigned type?"),
    list(id = "characterize", n = "3", title = "Characterize",
         entry = "calculateGeneShifts",
         question = "Which genes explain the difference?")
)

audit_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::div(
        class = "sc-page sc-audit",
        shiny::div(
            class = "sc-page-head",
            shiny::h2("Guided audit", class = "sc-page-title"),
            shiny::p(class = "sc-page-lede",
                     paste("The three steps the scDiagnostics case studies follow, run on",
                           "your current data with the package defaults."))
        ),
        shiny::uiOutput(ns("gate")),
        shiny::uiOutput(ns("controls")),
        shiny::uiOutput(ns("verdict")),
        shiny::uiOutput(ns("steps"))
    )
}

audit_server <- function(id, ctx_r, session_cache, jump_to = function(...) NULL,
                         open_diagnostic = function(...) NULL,
                         active = shiny::reactive(TRUE)) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        results <- shiny::reactiveValues(project = NULL, detect = NULL, characterize = NULL)
        ran <- shiny::reactiveVal(FALSE)

        output$gate <- shiny::renderUI({
            ctx <- ctx_r()
            if (isTRUE(ctx$ready) && !is.null(ctx$query)) return(NULL)
            shiny::div(
                alert(type = "info",
                      "The audit compares a query against a reference. Choose both on the Data page."),
                shiny::actionButton(ns("go_data"), "Go to Data", class = "btn-primary")
            )
        })
        shiny::observeEvent(input$go_data, jump_to("data"))

        output$controls <- shiny::renderUI({
            ctx <- ctx_r()
            shiny::req(isTRUE(ctx$ready), !is.null(ctx$query))
            focus <- audit_default_type(ctx)
            bslib::card(
                class = "sc-audit-controls",
                bslib::card_body(
                    shiny::div(
                        class = "row g-3 align-items-end",
                        shiny::div(
                            class = "col-md-4",
                            shiny::selectInput(
                                ns("cell_type"), "Focus cell type",
                                choices = ctx$shared_types, selected = focus
                            )
                        ),
                        shiny::div(
                            class = "col-md-4",
                            shiny::sliderInput(ns("pcs"), "Principal components",
                                               min = 1, max = max(3, ctx$n_pcs),
                                               value = range(ctx$pc_subset), step = 1,
                                               ticks = FALSE)
                        ),
                        shiny::div(
                            class = "col-md-4",
                            shiny::actionButton(ns("run"), "Run the audit",
                                                class = "btn-primary btn-lg w-100",
                                                icon = shiny::icon("play"))
                        )
                    )
                )
            )
        })

        audit_ctx <- shiny::reactive({
            ctx <- ctx_r()
            ctx$pc_subset <- if (is.null(input$pcs)) ctx$pc_subset else seq(input$pcs[1], input$pcs[2])
            ctx
        })

        run_step <- function(entry_id, extra = list(), ctx = NULL) {
            entry <- registry_get(entry_id)
            ctx <- ctx %||% audit_ctx()
            pv <- default_param_values(entry$params, ctx)
            pv <- utils::modifyList(pv, extra)
            args <- tryCatch(build_args(entry, ctx, pv),
                             error = function(e) sc_failure(conditionMessage(e)))
            if (is_failure(args)) return(list(value = args, ctx = ctx, params = pv, args = args))
            tier <- cache_for(ctx, session_cache)
            key <- cache_key(entry$id, ctx, pv)
            value <- with_cache(key, tier$cache, tier$shareable, run_entry(entry, args))
            list(value = value, ctx = ctx, params = pv, args = args, entry = entry)
        }

        shiny::observeEvent(input$run, {
            ctx <- audit_ctx()
            shiny::req(isTRUE(ctx$ready), !is.null(ctx$query))
            run_all(ctx)
            ran(TRUE)
        })

        run_all <- function(ctx) {
            prog <- shiny::Progress$new(session, min = 0, max = 3)
            on.exit(prog$close(), add = TRUE)
            prog$set(message = "Step 1 of 3: projecting onto the reference", value = 0)
            results$project <- run_step("plotCellTypePCA", ctx = ctx)
            prog$set(message = "Step 2 of 3: detecting anomalous cells", value = 1)
            results$detect <- run_step("detectAnomaly", ctx = ctx)
            prog$set(message = "Step 3 of 3: testing for gene shifts", value = 2)
            results$characterize <- run_step(
                "calculateGeneShifts",
                extra = list(pc_subset = utils::head(ctx$pc_subset, 3)), ctx = ctx)
            prog$set(value = 3)
            invisible(NULL)
        }

        # Run once automatically on first visit with data loaded, so the page
        # is never an empty form waiting to be told to do the obvious thing.
        shiny::observe({
            shiny::req(active(), !ran())
            ctx <- ctx_r()
            shiny::req(isTRUE(ctx$ready), !is.null(ctx$query), !is.null(input$cell_type))
            ran(TRUE)
            shiny::isolate(run_all(audit_ctx()))
        })

        # Recomputing is the user's call, but stale results must not linger
        # after the data changes. Keyed on the data signature rather than the
        # whole context: the context reactive also fires as inputs arrive and
        # when a slider moves, which would discard results just computed.
        shiny::observeEvent(context_signature(ctx_r()), {
            ran(FALSE)
            results$project <- NULL
            results$detect <- NULL
            results$characterize <- NULL
        }, ignoreInit = TRUE)

        # --- verdict --------------------------------------------------------
        output$verdict <- shiny::renderUI({
            det <- results$detect
            if (is.null(det) || is_failure(det$value)) return(NULL)
            ct <- input$cell_type
            v <- audit_verdict(det$value, results$characterize, ct)
            shiny::div(
                class = paste("sc-verdict", paste0("sc-verdict-", v$tone)),
                shiny::div(class = "sc-verdict-head",
                           bsicons::bs_icon(v$icon, size = "1.3rem"),
                           shiny::span(v$headline)),
                shiny::p(class = "sc-verdict-body", v$detail),
                if (length(v$bullets)) shiny::tags$ul(class = "sc-verdict-list",
                                                      lapply(v$bullets, shiny::tags$li))
            )
        })

        # --- steps -----------------------------------------------------------
        output$steps <- shiny::renderUI({
            ctx <- ctx_r()
            shiny::req(isTRUE(ctx$ready), !is.null(ctx$query))
            shiny::tagList(lapply(AUDIT_STEPS, function(s) {
                shiny::div(
                    class = "sc-audit-step",
                    shiny::div(
                        class = "sc-audit-step-head",
                        shiny::div(class = "sc-step-num", s$n),
                        shiny::div(
                            shiny::h3(class = "sc-audit-step-title", s$title),
                            shiny::p(class = "sc-audit-step-q", s$question)
                        )
                    ),
                    shiny::uiOutput(ns(paste0("figs_", s$id))),
                    shiny::uiOutput(ns(paste0("body_", s$id))),
                    shiny::uiOutput(ns(paste0("note_", s$id))),
                    shiny::div(
                        class = "sc-audit-step-foot",
                        shiny::actionLink(ns(paste0("open_", s$id)),
                                          shiny::tagList("Open the full ",
                                                         registry_get(s$entry)$label,
                                                         " panel ",
                                                         bsicons::bs_icon("arrow-right-short")))
                    )
                )
            }))
        })

        lapply(AUDIT_STEPS, function(s) {
            local({
                step <- s
                res_of <- function() results[[step$id]]

                # A full-height empty plot for each of three steps makes the page
                # mostly whitespace before the audit runs, so the placeholder is
                # compact and the plot area only appears once there is a plot.
                output[[paste0("body_", step$id)]] <- shiny::renderUI({
                    r <- res_of()
                    if (is.null(r)) {
                        return(empty_state("Press Run the audit to compute this step.", "play-circle"))
                    }
                    if (is_failure(r$value)) return(NULL)
                    shiny::div(class = "sc-plot-wrap",
                               shiny::plotOutput(session$ns(paste0("plot_", step$id)),
                                                 height = "520px"))
                })

                output[[paste0("plot_", step$id)]] <- shiny::renderPlot({
                    r <- res_of()
                    if (is.null(r)) return(blank_note("Press Run the audit."))
                    if (is_failure(r$value)) return(blank_note(r$value$message))
                    entry <- registry_get(step$entry)
                    pa <- audit_plot_args(step$id, r$value, input$cell_type)
                    p <- tryCatch(build_plot(entry, r$value, pa, r$ctx),
                                  error = function(e) e)
                    if (inherits(p, "error")) return(blank_note(conditionMessage(p)))
                    draw_result(p)
                }, res = 100)

                output[[paste0("figs_", step$id)]] <- shiny::renderUI({
                    r <- res_of()
                    if (is.null(r) || is_failure(r$value)) return(NULL)
                    entry <- registry_get(step$entry)
                    if (is.null(entry$summarise)) return(NULL)
                    figs <- tryCatch(entry$summarise(r$value), error = function(e) list())
                    if (!length(figs)) return(NULL)
                    bslib::layout_column_wrap(
                        width = "220px", fixed_width = FALSE, class = "sc-figures",
                        !!!lapply(figs, function(f) {
                            shiny::div(class = paste("sc-figure", tone_class(f$tone)),
                                       shiny::div(class = "sc-figure-label", f$label),
                                       shiny::div(class = "sc-figure-value", f$value),
                                       if (!is.null(f$caption)) {
                                           shiny::div(class = "sc-figure-caption", f$caption)
                                       })
                        })
                    )
                })

                output[[paste0("note_", step$id)]] <- shiny::renderUI({
                    r <- res_of()
                    if (is.null(r)) return(NULL)
                    if (is_failure(r$value)) {
                        return(error_panel(r$value$message,
                                           explain_failure(registry_get(step$entry), r)))
                    }
                    w <- attr(r$value, "sc_warnings")
                    if (!length(w)) return(NULL)
                    alert(type = "warning", shiny::tags$ul(lapply(unique(w), shiny::tags$li)))
                })

                # Guarded like every other link built inside a renderUI: a
                # re-render resets the counter to zero, which would otherwise
                # read as a click.
                open_key <- paste0("open_", step$id)
                shiny::observeEvent(input[[open_key]], {
                    if (!isTRUE(input[[open_key]] > 0)) return()
                    open_diagnostic(step$entry)
                    jump_to("diagnostics")
                }, ignoreInit = TRUE)
            })
        })

        always_render(output, c("gate", "controls", "verdict", "steps",
                                paste0("body_", vapply(AUDIT_STEPS, `[[`, character(1), "id")),
                                paste0("plot_", vapply(AUDIT_STEPS, `[[`, character(1), "id")),
                                paste0("figs_", vapply(AUDIT_STEPS, `[[`, character(1), "id")),
                                paste0("note_", vapply(AUDIT_STEPS, `[[`, character(1), "id"))))

        shiny::reactive(list(
            project = results$project,
            detect = results$detect,
            characterize = results$characterize,
            cell_type = input$cell_type
        ))
    })
}

# --- helpers ----------------------------------------------------------------

#' Default values for a set of parameters, given a context
default_param_values <- function(params, ctx) {
    out <- list()
    for (p in params) {
        out[[p$id]] <- switch(
            class(p)[1],
            sc_param_pcs = as_pc_subset(ctx$pc_subset),
            sc_param_celltypes = ctx$cell_types,
            p$default
        )
    }
    out[!vapply(out, is.null, logical(1))]
}

#' Which cell type should the audit focus on?
#'
#' The preset can name one (the COVID study is about CD14 monocytes); otherwise
#' take the largest shared population, which gives the most stable statistics.
audit_default_type <- function(ctx) {
    if (!length(ctx$shared_types)) return(NULL)
    if (!is.null(ctx$preset_id)) {
        want <- PRESETS[[ctx$preset_id]]$focus_cell_type
        if (!is.null(want) && want %in% ctx$shared_types) return(want)
    }
    counts <- table(as.character(ctx$query[[ctx$query_col]]))
    counts <- counts[names(counts) %in% ctx$shared_types]
    if (!length(counts)) return(ctx$shared_types[1])
    names(counts)[which.max(counts)]
}

audit_plot_args <- function(step_id, value, cell_type) {
    switch(step_id,
        detect = list(cell_type = if (!is.null(cell_type) && cell_type %in% names(value)) {
            cell_type
        } else "Combined", data_type = "query"),
        characterize = list(cell_type = cell_type, plot_type = "heatmap", n_genes = 12),
        list()
    )
}

#' Turn the audit results into a plain-language statement
#'
#' Deliberately hedged: these diagnostics flag candidates for inspection, and
#' the app should not imply that a flagged cell is a proven misannotation.
audit_verdict <- function(detect, characterize, cell_type) {
    tb <- anomaly_table(detect, "query_anomaly", "query_anomaly_scores")
    total <- sum(tb$flagged, na.rm = TRUE)
    n <- sum(tb$query_cells, na.rm = TRUE)
    rate <- if (n) total / n else NA_real_
    w <- worst_type(tb)

    bullets <- character(0)
    if (!is.null(w)) {
        bullets <- c(bullets, sprintf(
            "%s is the most affected cell type: %s of its %s query cells were flagged.",
            w$cell_type, pct(w$flagged_pct / 100), n_fmt(w$query_cells)))
    }
    if (!is.null(cell_type) && cell_type %in% tb$cell_type) {
        row <- tb[tb$cell_type == cell_type, ][1, ]
        bullets <- c(bullets, sprintf(
            "For the focus type %s, %s of %s query cells were flagged.",
            cell_type, n_fmt(row$flagged), n_fmt(row$query_cells)))
    }
    if (!is.null(characterize) && !is_failure(characterize$value)) {
        pcs <- grep("^PC[0-9]+$", names(characterize$value), value = TRUE)
        if (length(pcs)) {
            all_df <- do.call(rbind, lapply(pcs, function(p) characterize$value[[p]]))
            sig <- all_df[!is.na(all_df$significant) & all_df$significant, , drop = FALSE]
            if (nrow(sig)) {
                genes <- unique(sig[order(sig$p_adjusted), "gene"])
                top <- utils::head(genes, 5)
                bullets <- c(bullets, sprintf(
                    "Genes shifting most strongly between reference and query: %s.",
                    paste(top, collapse = ", ")))
            }
        }
    }

    if (!is.finite(rate)) {
        return(list(tone = "neutral", icon = "question-circle",
                    headline = "No verdict",
                    detail = "Anomaly detection returned no query scores.",
                    bullets = bullets))
    }
    if (rate < 0.05) {
        list(tone = "good", icon = "check-circle",
             headline = sprintf("Annotations look consistent with the reference (%s flagged)", pct(rate)),
             detail = paste(
                 "Few query cells are unusual relative to the reference cells sharing their label.",
                 "That is what a well-aligned annotation transfer looks like. Worth confirming with",
                 "variance attribution, which separates biological signal from batch effects."),
             bullets = bullets)
    } else if (rate < 0.15) {
        list(tone = "watch", icon = "exclamation-triangle",
             headline = sprintf("Some cells warrant a second look (%s flagged)", pct(rate)),
             detail = paste(
                 "A minority of query cells sit outside the reference distribution for their",
                 "assigned type. This is common and not automatically a problem: it can mean a",
                 "genuine cell state the reference lacks. The gene shifts below say which it is."),
             bullets = bullets)
    } else {
        list(tone = "concern", icon = "exclamation-octagon",
             headline = sprintf("A large share of the query looks unlike the reference (%s flagged)", pct(rate)),
             detail = paste(
                 "When flagging is this widespread, check whether it is concentrated in one cell",
                 "type or spread across all of them. Concentrated flagging points at that",
                 "annotation; uniform flagging usually means a batch effect, which variance",
                 "attribution will confirm."),
             bullets = bullets)
    }
}
