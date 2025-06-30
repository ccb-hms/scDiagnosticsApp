# Data Overview Tab UI
overview_tab <- tabItem(
    tabName = "overview",
    fluidRow(
        box(
            title = "Dataset Selection", status = "primary", solidHeader = TRUE, width = 12,
            fluidRow(
                column(6,
                       selectInput("reference_dataset", "Reference Dataset:",
                                   choices = list("Reference Marrow Myeloid Cells" = "reference_cells",
                                                  "Reference Marrow Myeloid Cells (Promonocytes Removed)" = "reference_cells_subset"),
                                   selected = "reference_cells")
                ),
                column(6,
                       selectInput("query_dataset", "Query Dataset:",
                                   choices = list("Query Marrow Myeloid Cells" = "query_cells"),
                                   selected = "query_cells")
                )
            )
        )
    ),
    fluidRow(
        box(
            title = "Reference Dataset Overview", status = "info", solidHeader = TRUE, width = 6,
            verbatimTextOutput("ref_summary"),
            br(),
            h4("Available Cell Type Columns:"),
            verbatimTextOutput("ref_cell_type_cols")
        ),
        box(
            title = "Query Dataset Overview", status = "success", solidHeader = TRUE, width = 6,
            verbatimTextOutput("query_summary"),
            br(),
            h4("Available Cell Type Columns:"),
            verbatimTextOutput("query_cell_type_cols")
        )
    )
)