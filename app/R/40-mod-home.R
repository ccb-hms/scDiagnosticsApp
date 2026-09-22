# Landing page.
#
# One sentence on what the app is for, one obvious way to start, and a short
# account of the workflow the package is built around.

home_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::div(
        class = "sc-page sc-home",
        shiny::div(
            class = "sc-hero",
            shiny::h1(class = "sc-hero-title", "Audit your cell type annotations"),
            shiny::p(
                class = "sc-hero-lede",
                "Automated annotation is only as trustworthy as the alignment between",
                "your reference and your query. This app runs the",
                shiny::tags$b("scDiagnostics"), "toolkit over both and shows you where",
                "the labels hold up and where they do not."
            ),
            shiny::div(
                class = "sc-hero-actions",
                shiny::actionButton(ns("demo"), "Try it on the COVID-19 case study",
                                    class = "btn-lg btn-primary",
                                    icon = shiny::icon("play")),
                shiny::actionButton(ns("own"), "Use my own data",
                                    class = "btn-lg btn-outline-primary",
                                    icon = shiny::icon("upload"))
            ),
            shiny::p(class = "sc-hero-note",
                     "Nothing is uploaded anywhere: your data stays in this session.")
        ),
        shiny::div(
            class = "sc-steps",
            shiny::h2("The workflow", class = "sc-section-title"),
            shiny::p(class = "sc-section-lede",
                     "The same three steps recur across every case study in the paper."),
            bslib::layout_column_wrap(
                width = "300px", fixed_width = FALSE, class = "sc-step-grid",
                step_card("1", "Project",
                          paste("Put the query into the reference's own PCA space and look",
                                "at each cell type. Alignment is a geometric question first.")),
                step_card("2", "Detect",
                          paste("Score every query cell against the reference cells sharing",
                                "its label, and flag the ones that do not belong.")),
                step_card("3", "Characterize",
                          paste("Find the genes that make the flagged cells different, so a",
                                "geometric anomaly becomes a biological statement."))
            ),
            shiny::div(
                class = "text-center mt-3",
                shiny::actionButton(ns("audit"), "Run these three steps now",
                                    class = "btn-primary", icon = shiny::icon("wand-magic-sparkles"))
            )
        ),
        shiny::div(
            class = "sc-catalogue",
            shiny::h2("Everything else", class = "sc-section-title"),
            shiny::p(class = "sc-section-lede",
                     paste("All", length(REGISTRY), "diagnostics in the package are available,",
                           "grouped the way the package documentation groups them.")),
            bslib::layout_column_wrap(
                width = "260px", fixed_width = FALSE, class = "sc-cat-grid",
                !!!lapply(setdiff(names(CATEGORIES), "utility"), function(k) {
                    cat_info <- CATEGORIES[[k]]
                    n <- length(registry_by_category(k))
                    shiny::div(
                        class = "sc-cat-card",
                        bsicons::bs_icon(cat_info$icon, size = "1.4rem"),
                        shiny::div(class = "sc-cat-title", cat_info$label),
                        shiny::div(class = "sc-cat-count", sprintf("%d diagnostic%s", n,
                                                                   if (n == 1) "" else "s")),
                        shiny::div(class = "sc-cat-blurb", cat_info$blurb)
                    )
                })
            )
        ),
        shiny::div(
            class = "sc-cite",
            shiny::h2("Citation", class = "sc-section-title"),
            shiny::p(
                "Christidis A, Ghazi A, Chawla S, Turaga N, Gentleman R, Geistlinger L (2026).",
                shiny::tags$em(paste("scDiagnostics: systematic assessment of cell type",
                                             "annotation in single-cell transcriptomics data.")),
                "Briefings in Bioinformatics, 27(5), bbag496.",
                shiny::a("doi:10.1093/bib/bbag496", href = "https://doi.org/10.1093/bib/bbag496",
                         target = "_blank", rel = "noopener")
            )
        )
    )
}

step_card <- function(n, title, text) {
    shiny::div(
        class = "sc-step",
        shiny::div(class = "sc-step-num", n),
        shiny::div(class = "sc-step-body",
                   shiny::div(class = "sc-step-title", title),
                   shiny::p(class = "sc-step-text", text))
    )
}

home_server <- function(id, jump_to, load_preset) {
    shiny::moduleServer(id, function(input, output, session) {
        shiny::observeEvent(input$demo, {
            load_preset("covid")
            jump_to("audit")
        })
        shiny::observeEvent(input$own, jump_to("data"))
        shiny::observeEvent(input$audit, jump_to("audit"))
    })
}
