# Marker Expression Tab UI
marker_tab <- tabItem(
    tabName = "marker",
    fluidRow(
        box(
            title = "Marker Expression Parameters", status = "primary", solidHeader = TRUE, width = 4,
            selectInput("marker_ref_cell_type_col", "Reference Cell Type Column:",
                        choices = NULL),
            selectInput("marker_query_cell_type_col", "Query Cell Type Column:",
                        choices = NULL),
            selectInput("marker_assay_name", "Assay:",
                        choices = c("logcounts", "counts"),
                        selected = "logcounts"),
            br(),
            h4("Gene and Cell Type Selection:"),
            div(
                style = "margin-bottom: 10px;",
                selectizeInput("gene_name", "Gene Name:", 
                               choices = NULL,
                               options = list(
                                   placeholder = "Start typing to search for genes...",
                                   maxOptions = 50,
                                   searchField = c('value', 'text')
                               ))
            ),
            helpText("Search and select from available genes in both datasets"),
            
            selectInput("marker_cell_type", "Cell Type:",
                        choices = NULL,
                        multiple = FALSE),
            helpText("Select one cell type for cell type-specific analysis"),
            br(),
            h4("Normalization Options:"),
            selectInput("normalization", "Expression Normalization:",
                        choices = list(
                            "Z-Score (Standard)" = "z_score",
                            "Min-Max (0-1 scale)" = "min_max", 
                            "Quantile Rank (0-100)" = "rank",
                            "None (Original values)" = "none"
                        ),
                        selected = "z_score"),
            helpText("Z-score: Centers and scales within each dataset"),
            helpText("Min-Max: Scales to 0-1 range within each dataset"),
            helpText("Rank: Converts to percentile ranks (0-100)"),
            helpText("None: Preserves original expression scale"),
            br(),
            actionButton("generate_marker", "Generate Expression Plot", 
                         class = "btn-primary", width = "100%"),
            br(), br(),
            downloadButton("download_marker_plot", "Download Plot", 
                           class = "btn-success", width = "100%")
        ),
        box(
            title = "Interactive Marker Expression Analysis", status = "success", solidHeader = TRUE, width = 8,
            conditionalPanel(
                condition = "input.generate_marker == 0",
                div(style = "text-align: center; padding: 50px;",
                    h3("Select a gene and cell type, then click 'Generate Expression Plot'", 
                       style = "color: #999;"))
            ),
            conditionalPanel(
                condition = "input.generate_marker > 0",
                withSpinner(plotOutput("marker_plot", height = "600px"), color = "#0dc5c1")
            )
        )
    ),
    fluidRow(
        box(
            title = "Expression Statistics", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_marker > 0",
                h4("Gene Expression Summary"),
                verbatimTextOutput("marker_stats")
            )
        ),
        box(
            title = "Analysis Summary", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_marker > 0",
                verbatimTextOutput("marker_analysis_summary")
            )
        )
    )
)