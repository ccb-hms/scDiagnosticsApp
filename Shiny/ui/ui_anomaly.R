# Anomaly Detection Tab UI
anomaly_tab <- tabItem(
    tabName = "anomaly",
    fluidRow(
        box(
            title = "Anomaly Detection Parameters", status = "primary", solidHeader = TRUE, width = 4,
            selectInput("anomaly_ref_cell_type_col", "Reference Cell Type Column:",
                        choices = NULL),
            selectInput("anomaly_query_cell_type_col", "Query Cell Type Column:",
                        choices = NULL),
            selectInput("anomaly_assay_name", "Assay:",
                        choices = c("logcounts", "counts"),
                        selected = "logcounts"),
            checkboxGroupInput("anomaly_pc_subset", "Principal Components:",
                               choices = setNames(1:10, paste0("PC", 1:10)),
                               selected = 1:5),
            selectInput("anomaly_cell_types", "Cell Types (leave empty for all):",
                        choices = NULL,
                        multiple = TRUE),
            br(),
            h4("Isolation Forest Parameters:"),
            numericInput("n_tree", "Number of Trees:", value = 500, min = 100, max = 2000, step = 100),
            numericInput("anomaly_threshold", "Anomaly Threshold:", value = 0.6, min = 0.1, max = 0.9, step = 0.05),
            br(),
            h4("Visualization Options:"),
            selectInput("anomaly_cell_type_plot", "Cell Type to Plot:",
                        choices = NULL),
            selectInput("anomaly_data_type", "Data Type:",
                        choices = c("Query" = "query", "Reference" = "reference"),
                        selected = "query"),
            selectInput("anomaly_upper_facet", "Upper Panels:",
                        choices = c("blank", "contour", "ellipse"),
                        selected = "blank"),
            selectInput("anomaly_diagonal_facet", "Diagonal Panels:",
                        choices = c("density", "ridge", "boxplot", "blank"),
                        selected = "density"),
            br(),
            actionButton("generate_anomaly", "Generate Anomaly Analysis", 
                         class = "btn-primary", width = "100%"),
            br(), br(),
            downloadButton("download_anomaly_plot", "Download Plot", 
                           class = "btn-success", width = "100%")
        ),
        box(
            title = "Interactive Anomaly Detection Analysis", status = "success", solidHeader = TRUE, width = 8,
            conditionalPanel(
                condition = "input.generate_anomaly == 0",
                div(style = "text-align: center; padding: 50px;",
                    h3("Configure parameters and click 'Generate Anomaly Analysis' to start", 
                       style = "color: #999;"))
            ),
            conditionalPanel(
                condition = "input.generate_anomaly > 0",
                withSpinner(plotOutput("anomaly_plot", height = "700px"), color = "#0dc5c1")
            )
        )
    ),
    fluidRow(
        box(
            title = "Anomaly Statistics", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_anomaly > 0",
                h4("Anomaly Detection Results"),
                verbatimTextOutput("anomaly_stats")
            )
        ),
        box(
            title = "Analysis Summary", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_anomaly > 0",
                verbatimTextOutput("anomaly_analysis_summary")
            )
        )
    )
)