# The diagnostics catalogue.
#
# A grouped list on the left, one panel on the right. The panel is rendered
# lazily for the selected diagnostic only.
#
# Rendering all 33 panels into hidden tabs is the obvious alternative and it
# does not work: outputs inside a hidden tab are suspended, and a nested tab
# set does not reliably resume them when its parent becomes visible, so panels
# arrive empty. Rendering one panel at a time also keeps the initial HTML
# small.

catalogue_ui <- function(id) {
    ns <- shiny::NS(id)
    bslib::layout_sidebar(
        sidebar = bslib::sidebar(
            width = 280, title = "Diagnostics", open = "desktop",
            shiny::div(
                class = "sc-catalogue-search",
                shiny::textInput(ns("filter"), NULL, placeholder = "Filter…")
            ),
            shiny::uiOutput(ns("list"))
        ),
        shiny::uiOutput(ns("panel"))
    )
}

catalogue_server <- function(id, ctx_r, session_cache, nav_active) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        selected <- shiny::reactiveVal(names(REGISTRY)[1])

        # Which entries are unavailable on the current data. Depending on this
        # rather than on the whole context keeps the list from re-rendering
        # every time an unrelated setting changes.
        gated_ids <- shiny::reactive({
            ctx <- ctx_r()
            names(Filter(function(e) !is.null(entry_ready(e, ctx)), REGISTRY))
        })

        output$list <- shiny::renderUI({
            q <- tolower(input$filter %|% "")
            gated <- gated_ids()
            active_id <- selected()
            sections <- lapply(names(CATEGORIES), function(ck) {
                entries <- registry_by_category(ck)
                if (nzchar(q)) {
                    entries <- Filter(function(e) {
                        grepl(q, tolower(paste(e$label, e$fn, e$blurb, e$question)), fixed = TRUE)
                    }, entries)
                }
                if (!length(entries)) return(NULL)
                shiny::div(
                    class = "sc-cat-group",
                    shiny::div(class = "sc-cat-group-title",
                               bsicons::bs_icon(CATEGORIES[[ck]]$icon),
                               " ", CATEGORIES[[ck]]$label),
                    lapply(entries, function(e) {
                        shiny::actionLink(
                            ns(paste0("pick_", e$id)),
                            shiny::tagList(
                                shiny::span(class = "sc-cat-item-label", e$label),
                                if (e$heavy) shiny::span(class = "sc-cat-item-slow", "●"),
                                if (e$id %in% gated) {
                                    shiny::span(class = "sc-cat-item-gated",
                                                title = "Not applicable to the current data",
                                                bsicons::bs_icon("slash-circle"))
                                }
                            ),
                            class = paste("sc-cat-item",
                                          if (identical(active_id, e$id)) "sc-cat-item-active")
                        )
                    })
                )
            })
            sections <- not_null(sections)
            if (!length(sections)) return(empty_state("Nothing matches that filter.", "search"))
            shiny::div(class = "sc-cat-list", sections)
        })

        # One observer per entry, created once.
        #
        # The guard on the click count matters: re-rendering the list resets
        # every actionLink to zero, which observeEvent sees as a change and
        # would otherwise treat as 33 simultaneous clicks.
        for (eid in names(REGISTRY)) {
            local({
                this <- eid
                key <- paste0("pick_", this)
                shiny::observeEvent(input[[key]], {
                    if (isTRUE(input[[key]] > 0)) selected(this)
                }, ignoreInit = TRUE)
            })
        }

        # Build the panel only once the tab has actually been opened, and keep
        # it thereafter. Inserting the markup while the tab is visible means
        # its outputs bind in a visible state, which is what makes them render
        # at all: a nested output created inside a hidden bslib tab is
        # suspended and is not reliably resumed when the tab is shown.
        opened <- shiny::reactiveVal(FALSE)
        shiny::observeEvent(nav_active(), {
            if (isTRUE(nav_active())) opened(TRUE)
        }, ignoreInit = FALSE)

        output$panel <- shiny::renderUI({
            shiny::req(opened())
            entry <- REGISTRY[[selected()]]
            diagnostic_ui(ns(paste0("d_", entry$id)), entry)
        })

        # Module servers for every diagnostic, each idle unless it is the one
        # on screen. Creating them up front costs nothing; they hold only
        # reactives until activated.
        for (eid in names(REGISTRY)) {
            local({
                this <- eid
                diagnostic_server(
                    paste0("d_", this), REGISTRY[[this]], ctx_r, session_cache,
                    active = shiny::reactive(nav_active() && identical(selected(), this))
                )
            })
        }

        always_render(output, c("list", "panel"))

        # Let other pages jump straight to a diagnostic.
        list(select = function(entry_id) {
            if (entry_id %in% names(REGISTRY)) selected(entry_id)
        })
    })
}
