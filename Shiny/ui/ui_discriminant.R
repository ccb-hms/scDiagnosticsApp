# Discriminant Projection Tab UI
discriminant_tab <- tabItem(
    tabName = "discriminant",
    fluidRow(
        box(
            title = "Discriminant Projection Parameters", status = "primary", solidHeader = TRUE, width = 4,
            selectInput("disc_ref_cell_type_col", "Reference Cell Type Column:",
                        choices = NULL),
            selectInput("disc_query_cell_type_col", "Query Cell Type Column:",
                        choices = NULL),
            selectInput("disc_assay_name", "Assay:",
                        choices = c("logcounts", "counts"),
                        selected = "logcounts"),
            selectInput("disc_cell_types", "Cell Types (leave empty for all):",
                        choices = NULL,
                        multiple = TRUE),
            br(),
            h4("Random Forest Parameters:"),
            numericInput("disc_n_tree", "Number of Trees:", value = 500, min = 100, max = 2000, step = 100),
            numericInput("disc_n_top", "Top Variables per Comparison:", value = 20, min = 5, max = 100, step = 5),
            br(),
            h4("Discriminant Analysis Parameters:"),
            numericInput("eigen_threshold", "Eigenvalue Threshold:", value = 0.1, min = 0.001, max = 1, step = 0.01),
            checkboxInput("calculate_metrics", "Calculate Similarity Metrics", value = FALSE),
            conditionalPanel(
                condition = "input.calculate_metrics == true",
                numericInput("alpha", "Significance Level (α):", value = 0.01, min = 0.001, max = 0.1, step = 0.001)
            ),
            br(),
            h4("Visualization Options:"),
            checkboxGroupInput("dv_subset", "Discriminant Vectors:",
                               choices = setNames(1:10, paste0("DV", 1:10)),
                               selected = 1:3),
            selectInput("disc_lower_facet", "Lower Panels:",
                        choices = c("scatter", "contour", "ellipse", "blank"),
                        selected = "scatter"),
            selectInput("disc_diagonal_facet", "Diagonal Panels:",
                        choices = c("ridge", "density", "boxplot", "blank"),
                        selected = "ridge"),
            selectInput("disc_upper_facet", "Upper Panels:",
                        choices = c("blank", "scatter", "contour", "ellipse"),
                        selected = "blank"),
            br(),
            actionButton("generate_discriminant", "Generate Discriminant Analysis", 
                         class = "btn-primary", width = "100%"),
            br(), br(),
            downloadButton("download_discriminant_plot", "Download Plot", 
                           class = "btn-success", width = "100%")
        ),
        box(
            title = "Interactive Discriminant Projection Analysis", status = "success", solidHeader = TRUE, width = 8,
            conditionalPanel(
                condition = "input.generate_discriminant == 0",
                div(style = "text-align: center; padding: 50px;",
                    h3("Configure parameters and click 'Generate Discriminant Analysis' to start", 
                       style = "color: #999;"))
            ),
            conditionalPanel(
                condition = "input.generate_discriminant > 0",
                withSpinner(plotOutput("discriminant_plot", height = "700px"), color = "#0dc5c1")
            )
        )
    ),
    fluidRow(
        box(
            title = "Discriminant Vectors & Eigenvalues", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_discriminant > 0",
                h4("Discriminant Analysis Results"),
                verbatimTextOutput("discriminant_stats")
            )
        ),
        box(
            title = "Analysis Summary", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_discriminant > 0",
                verbatimTextOutput("discriminant_analysis_summary")
            )
        )
    )
)