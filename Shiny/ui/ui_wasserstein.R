# Wasserstein Distance Tab UI
wasserstein_tab <- tabItem(
    tabName = "wasserstein",
    fluidRow(
        box(
            title = "Wasserstein Distance Parameters", status = "primary", solidHeader = TRUE, width = 4,
            selectInput("wass_ref_cell_type_col", "Reference Cell Type Column:",
                        choices = NULL),
            selectInput("wass_query_cell_type_col", "Query Cell Type Column:",
                        choices = NULL),
            selectInput("wass_assay_name", "Assay:",
                        choices = c("logcounts", "counts"),
                        selected = "logcounts"),
            checkboxGroupInput("wass_pc_subset", "Principal Components:",
                               choices = setNames(1:15, paste0("PC", 1:15)),
                               selected = 1:10),
            selectInput("wass_cell_types", "Cell Types (leave empty for all):",
                        choices = NULL,
                        multiple = TRUE),
            br(),
            h4("Analysis Parameters:"),
            numericInput("n_resamples", "Number of Resamples:", 
                         value = 200, min = 50, max = 1000, step = 50),
            br(),
            h4("Visualization Options:"),
            selectInput("wass_plot_cell_types", "Cell Types to Plot:",
                        choices = NULL,
                        multiple = TRUE),
            helpText("Leave empty to plot all analyzed cell types"),
            br(),
            actionButton("generate_wasserstein", "Generate Wasserstein Analysis", 
                         class = "btn-primary", width = "100%"),
            br(), br(),
            downloadButton("download_wasserstein_plot", "Download Plot", 
                           class = "btn-success", width = "100%")
        ),
        box(
            title = "Interactive Wasserstein Distance Analysis", status = "success", solidHeader = TRUE, width = 8,
            conditionalPanel(
                condition = "input.generate_wasserstein == 0",
                div(style = "text-align: center; padding: 50px;",
                    h3("Configure parameters and click 'Generate Wasserstein Analysis' to start", 
                       style = "color: #999;"))
            ),
            conditionalPanel(
                condition = "input.generate_wasserstein > 0",
                withSpinner(plotOutput("wasserstein_plot", height = "700px"), color = "#0dc5c1")
            )
        )
    ),
    fluidRow(
        box(
            title = "Probability of Superiority", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_wasserstein > 0",
                h4("Wasserstein Distance Results"),
                verbatimTextOutput("wasserstein_stats")
            )
        ),
        box(
            title = "Analysis Summary", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_wasserstein > 0",
                verbatimTextOutput("wasserstein_analysis_summary")
            )
        )
    )
)