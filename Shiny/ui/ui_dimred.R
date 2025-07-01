# Gene Expression Dimensional Reduction Tab UI
dimred_tab <- tabItem(
    tabName = "dimred",
    fluidRow(
        box(
            title = "Gene Expression Dimensional Reduction Parameters", status = "primary", solidHeader = TRUE, width = 4,
            h4("Dataset Selection:"),
            selectInput("dimred_dataset", "Dataset:",
                        choices = list(
                            "Query Marrow Myeloid Cells" = "query",
                            "Reference Marrow Myeloid Cells" = "reference_cells",
                            "Reference Marrow Myeloid Cells (Promonocytes Removed)" = "reference_cells_subset"
                        ),
                        selected = "query"),
            
            selectInput("dimred_cell_type_col", "Cell Type Column:",
                        choices = NULL),
            selectInput("dimred_assay_name", "Assay:",
                        choices = c("logcounts", "counts"),
                        selected = "logcounts"),
            br(),
            h4("Visualization Method:"),
            selectInput("dimred_method", "Reduction Method:",
                        choices = list(
                            "t-SNE" = "TSNE",
                            "UMAP" = "UMAP", 
                            "PCA" = "PCA"
                        ),
                        selected = "UMAP"),
            
            conditionalPanel(
                condition = "input.dimred_method == 'PCA'",
                checkboxGroupInput("dimred_pc_subset", "Principal Components:",
                                   choices = setNames(1:10, paste0("PC", 1:10)),
                                   selected = 1:5)
            ),
            br(),
            h4("Gene and Cell Type Selection:"),
            selectizeInput("dimred_gene_name", "Gene Name:", 
                           choices = NULL,
                           options = list(
                               placeholder = "Start typing to search for genes...",
                               maxOptions = 50,
                               searchField = c('value', 'text')
                           )),
            helpText("Search and select from available genes in the dataset"),
            
            selectInput("dimred_cell_types", "Cell Types (leave empty for all):",
                        choices = NULL,
                        multiple = TRUE),
            helpText("Leave empty to show all cell types, or select specific types"),
            br(),
            actionButton("generate_dimred", "Generate Dimensional Reduction Plot", 
                         class = "btn-primary", width = "100%"),
            br(), br(),
            downloadButton("download_dimred_plot", "Download Plot", 
                           class = "btn-success", width = "100%")
        ),
        box(
            title = "Interactive Gene Expression Dimensional Reduction", status = "success", solidHeader = TRUE, width = 8,
            conditionalPanel(
                condition = "input.generate_dimred == 0",
                div(style = "text-align: center; padding: 50px;",
                    h3("Select a dataset, gene, and visualization method, then click 'Generate Dimensional Reduction Plot'", 
                       style = "color: #999;"))
            ),
            conditionalPanel(
                condition = "input.generate_dimred > 0",
                withSpinner(plotOutput("dimred_plot", height = "700px"), color = "#0dc5c1")
            )
        )
    ),
    fluidRow(
        box(
            title = "Expression Summary", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_dimred > 0",
                h4("Gene Expression Statistics"),
                verbatimTextOutput("dimred_stats")
            )
        ),
        box(
            title = "Analysis Summary", status = "info", solidHeader = TRUE, width = 6,
            conditionalPanel(
                condition = "input.generate_dimred > 0",
                verbatimTextOutput("dimred_analysis_summary")
            )
        )
    )
)