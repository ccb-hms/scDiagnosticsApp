# scDiagnostics interactive app
#
# Companion to Christidis et al. (2026), Briefings in Bioinformatics 27(5),
# bbag496. The app calls the published scDiagnostics package directly; it holds
# no copy of any statistical method.
#
# Files in R/ are sourced automatically by Shiny, in alphabetical order.

suppressPackageStartupMessages({
    library(shiny)
    library(bslib)
    library(bsicons)
    library(scDiagnostics)
    library(SingleCellExperiment)
    library(ggplot2)
    library(promises)
    library(future)
})

options(
    shiny.maxRequestSize = as.numeric(Sys.getenv("SCDIAG_MAX_UPLOAD_MB", "500")) * 1024^2,
    shiny.sanitize.errors = FALSE
)

# The app theme applies to package plots too, without touching the package.
ggplot2::theme_set(sc_ggtheme())

# Background workers, started once per R process. If they cannot start, the
# slow diagnostics simply run in the foreground.
options(sc.async.ready = init_async())

# Refuse to start against a package version whose results would differ from
# what the app claims to be showing. Better a clear message at startup than
# figures that quietly disagree with the paper.
PKG_VERSION_PROBLEM <- check_package_version()
if (!is.null(PKG_VERSION_PROBLEM)) warning(PKG_VERSION_PROBLEM, call. = FALSE)

# --- UI --------------------------------------------------------------------

ui <- bslib::page_navbar(
    title = shiny::tagList(
        shiny::span(class = "sc-brand-mark", "scD"),
        shiny::span(class = "sc-brand-name", "scDiagnostics")
    ),
    id = "nav",
    theme = sc_theme(),
    fillable = FALSE,
    header = shiny::tagList(
        shiny::tags$head(
            shiny::tags$link(rel = "stylesheet", href = "styles.css"),
            shiny::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1")
        ),
        shiny::uiOutput("context_bar")
    ),
    bslib::nav_panel("Home", value = "home", home_ui("home")),
    bslib::nav_panel("Data", value = "data", data_ui("data")),
    bslib::nav_panel("Audit", value = "audit", audit_ui("audit")),
    bslib::nav_panel("Diagnostics", value = "diagnostics", catalogue_ui("catalogue")),
    bslib::nav_panel("Report", value = "report", report_ui("report")),
    bslib::nav_spacer(),
    bslib::nav_item(
        shiny::a(bsicons::bs_icon("github"), " Source",
                 href = "https://github.com/ccb-hms/scDiagnosticsApp",
                 target = "_blank", rel = "noopener", class = "nav-link")
    ),
    bslib::nav_panel("About", value = "about", about_ui())
)

# --- server ----------------------------------------------------------------

server <- function(input, output, session) {

    # Results from uploaded data never touch the shared disk cache.
    session_cache <- new_session_cache()

    jump_to <- function(target) {
        bslib::nav_select("nav", target, session = session)
    }

    data_mod <- data_server("data", jump_to = jump_to)
    ctx_r <- data_mod$ctx

    home_server("home", jump_to = jump_to, load_preset = data_mod$load_preset)

    # A persistent bar naming the current data, so no panel is ambiguous about
    # what it is showing.
    output$context_bar <- shiny::renderUI({
        if (!is.null(PKG_VERSION_PROBLEM)) {
            return(shiny::div(class = "sc-context-bar sc-context-bar-alert",
                              bsicons::bs_icon("exclamation-triangle-fill"),
                              shiny::span(class = "sc-context-text", PKG_VERSION_PROBLEM)))
        }
        ctx <- ctx_r()
        if (!isTRUE(ctx$ready)) return(NULL)
        shiny::div(
            class = "sc-context-bar",
            bsicons::bs_icon("database"),
            shiny::span(class = "sc-context-text", context_caption(ctx)),
            shiny::actionLink("change_data", "Change", class = "sc-context-link")
        )
    })
    shiny::observeEvent(input$change_data, jump_to("data"))

    # Each diagnostic computes only while its own panel is on screen.
    catalogue <- catalogue_server(
        "catalogue", ctx_r, session_cache,
        nav_active = shiny::reactive(identical(input$nav, "diagnostics"))
    )

    audit_server("audit", ctx_r, session_cache, jump_to = jump_to,
                 open_diagnostic = catalogue$select,
                 active = shiny::reactive(identical(input$nav, "audit")))

    report_server("report", ctx_r, session_cache,
                  active = shiny::reactive(identical(input$nav, "report")))
}

shiny::shinyApp(ui, server)
