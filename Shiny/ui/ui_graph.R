# Graph Integration Tab UI
graph_tab <- tabItem(
    tabName = "graph",
    fluidRow(
        box(
            title = "Graph Integration Parameters", status = "primary", solidHeader = TRUE, width = 4,
            selectInput("graph_ref_cell_type_col", "Reference Cell Type Column:",
                        choices = NULL),
            selectInput("graph_query_cell_type_col", "Query Cell Type Column:",
                        choices = NULL),
            selectInput("graph_assay_name", "Assay:",
                        choices = c("logcounts", "counts"),
                        selected = "logcounts"),
            checkboxGroupInput("graph_pc_subset", "Principal Components:",
                               choices = setNames(1:15, paste0("PC", 1:15)),
                               selected = 1:10),
            selectInput("graph_cell_types", "Cell Types (leave empty for all):",
                        choices = NULL,
                        multiple = TRUE),
            br(),
            h4("Graph Parameters:"),
            numericInput("k_neighbors", "K Neighbors:", value = 30, min = 5, max = 100, step = 5),
            numericInput("resolution", "Clustering Resolution:", value = 0.1, min = 0.01, max = 1, step = 0.01),
            numericInput("min_cells_per_community", "Min Cells per Community:", value = 10, min = 5, max = 50, step = 5),
            br(),
            h4("Threshold Parameters:"),
            numericInput("high_query_prop_threshold", "Query-Only Threshold:", value = 0.9, min = 0.5, max = 1, step = 0.05),
            numericInput("cross_type_threshold", "Cross-Type Threshold:", value = 0.15, min = 0.05, max = 0.5, step = 0.05),
            numericInput("local_consistency_threshold", "Local Consistency Threshold:", value = 0.6, min = 0.1, max = 1, step = 0.05),
            br(),
            actionButton("generate_graph", "Generate Graph Analysis", 
                         class = "btn-primary", width = "100%"),
            br(), br(),
            downloadButton("download_graph_plot", "Download Plot", 
                           class = "btn-success", width = "100%")
        ),
        box(
            title = "Interactive Graph Integration Analysis", status = "success", solidHeader = TRUE, width = 8,
            conditionalPanel(
                condition = "input.generate_graph == 0",
                div(style = "text-align: center; padding: 50px;",
                    h3("Configure parameters and click 'Generate Graph Analysis' to start", 
                       style = "color: #999;"))
            ),
            conditionalPanel(
                condition = "input.generate_graph > 0",
                fluidRow(
                    column(4,
                           selectInput("graph_plot_type", "Plot Type:",
                                       choices = c("Community Network" = "community_network",
                                                   "Summary" = "summary", 
                                                   "Community Data" = "community_data"),
                                       selected = "community_network")
                    ),
                    column(4,
                           conditionalPanel(
                               condition = "input.graph_plot_type == 'community_network'",
                               selectInput("graph_color_by", "Color By:",
                                           choices = c("Cell Type" = "cell_type",
                                                       "Community Type" = "community_type"),
                                           selected = "cell_type")
                           )
                    ),
                    column(4,
                           checkboxInput("exclude_reference_only", "Exclude Reference-Only Communities", value = FALSE)
                    )
                ),
                withSpinner(plotOutput("graph_plot", height = "600px"), color = "#0dc5c1")
            )
        )
    ),
    fluidRow(
        box(
            title = "Community Details", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_graph > 0",
                h4("Community Analysis Results"),
                verbatimTextOutput("community_details")
            )
        ),
        box(
            title = "Analysis Summary", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_graph > 0",
                verbatimTextOutput("graph_analysis_summary")
            )
        )
    )
)