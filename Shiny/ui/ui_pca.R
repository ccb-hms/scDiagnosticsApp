# PCA Projection Tab UI
pca_tab <- tabItem(
    tabName = "pca",
    fluidRow(
        box(
            title = "PCA Parameters", status = "primary", solidHeader = TRUE, width = 4,
            selectInput("ref_cell_type_col", "Reference Cell Type Column:",
                        choices = NULL),
            selectInput("query_cell_type_col", "Query Cell Type Column:",
                        choices = NULL),
            selectInput("assay_name", "Assay:",
                        choices = c("logcounts", "counts"),
                        selected = "logcounts"),
            checkboxGroupInput("pc_subset", "Principal Components:",
                               choices = setNames(1:10, paste0("PC", 1:10)),
                               selected = 1:5),
            selectInput("cell_types", "Cell Types (leave empty for all):",
                        choices = NULL,
                        multiple = TRUE),
            br(),
            h4("Plot Options:"),
            selectInput("lower_facet", "Lower Panels:",
                        choices = c("scatter", "contour", "ellipse", "blank"),
                        selected = "scatter"),
            selectInput("diagonal_facet", "Diagonal Panels:",
                        choices = c("ridge", "density", "boxplot"),
                        selected = "ridge"),
            selectInput("upper_facet", "Upper Panels:",
                        choices = c("blank", "scatter", "contour", "ellipse"),
                        selected = "blank"),
            br(),
            actionButton("generate_plot", "Generate PCA Plot", 
                         class = "btn-primary", width = "100%"),
            br(), br(),
            downloadButton("download_pca_plot", "Download Plot", 
                           class = "btn-success", width = "100%")
        ),
        box(
            title = "Interactive PCA Projection", status = "success", solidHeader = TRUE, width = 8,
            conditionalPanel(
                condition = "input.generate_plot == 0",
                div(style = "text-align: center; padding: 50px;",
                    h3("Configure parameters and click 'Generate PCA Plot' to start", 
                       style = "color: #999;"))
            ),
            conditionalPanel(
                condition = "input.generate_plot > 0",
                withSpinner(plotOutput("pca_plot", height = "700px"), 
                            color = "#0dc5c1")
            )
        )
    ),
    fluidRow(
        box(
            title = "Plot Information & Settings Summary", status = "info", solidHeader = TRUE, width = 12,
            verbatimTextOutput("plot_info")
        )
    )
)