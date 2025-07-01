# Data Overview Tab UI
overview_tab <- tabItem(
    tabName = "overview",
    fluidRow(
        box(
            title = "Dataset Selection", status = "primary", solidHeader = TRUE, width = 12,
            h4("Choose Your Datasets for Analysis"),
            p("Select the reference and query datasets you want to use for all diagnostic analyses. Built-in datasets are always available, and any uploaded datasets will appear here as well."),
            br(),
            
            fluidRow(
                column(6,
                       wellPanel(
                           h4("Reference Dataset", style = "color: #5cb85c;"),
                           uiOutput("reference_dataset_selection"),
                           helpText("The reference dataset serves as the baseline for comparisons")
                       )
                ),
                column(6,
                       wellPanel(
                           h4("Query Dataset", style = "color: #f0ad4e;"),
                           uiOutput("query_dataset_selection"),
                           helpText("The query dataset will be compared against the reference")
                       )
                )
            )
        )
    ),
    
    fluidRow(
        box(
            title = "Reference Dataset Summary", status = "success", solidHeader = TRUE, width = 6,
            verbatimTextOutput("ref_summary")
        ),
        box(
            title = "Query Dataset Summary", status = "warning", solidHeader = TRUE, width = 6,
            verbatimTextOutput("query_summary")
        )
    ),
    
    fluidRow(
        box(
            title = "Reference Cell Type Columns", status = "info", solidHeader = TRUE, width = 6,
            verbatimTextOutput("ref_cell_type_cols")
        ),
        box(
            title = "Query Cell Type Columns", status = "info", solidHeader = TRUE, width = 6,
            verbatimTextOutput("query_cell_type_cols")
        )
    ),
    
    fluidRow(
        box(
            title = "Dataset Compatibility Analysis", status = "primary", solidHeader = TRUE, width = 12,
            verbatimTextOutput("compatibility_check")
        )
    )
)