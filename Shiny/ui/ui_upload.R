# Data Upload Tab UI
upload_tab <- tabItem(
    tabName = "upload",
    fluidRow(
        box(
            title = "Upload SingleCellExperiment Objects", status = "primary", solidHeader = TRUE, width = 12,
            h4("Upload Your Own Data", style = "color: #337ab7;"),
            p("Upload SingleCellExperiment objects saved as .rds files. Both reference and query datasets are optional - you can upload one or both."),
            br(),
            
            fluidRow(
                column(6,
                       wellPanel(
                           h4("Reference Dataset", style = "color: #5cb85c;"),
                           fileInput("upload_reference", "Choose Reference .rds File:",
                                     accept = ".rds",
                                     placeholder = "No file selected"),
                           textInput("reference_name", "Dataset Name:",
                                     placeholder = "e.g., My Reference Dataset"),
                           conditionalPanel(
                               condition = "output.reference_uploaded",
                               div(style = "margin-top: 10px; padding: 10px; background-color: #dff0d8; border: 1px solid #d6e9c6; border-radius: 4px;",
                                   h5("✓ Reference Dataset Loaded", style = "color: #3c763d; margin: 0;"),
                                   verbatimTextOutput("reference_summary_upload")
                               )
                           )
                       )
                ),
                column(6,
                       wellPanel(
                           h4("Query Dataset", style = "color: #f0ad4e;"),
                           fileInput("upload_query", "Choose Query .rds File:",
                                     accept = ".rds",
                                     placeholder = "No file selected"),
                           textInput("query_name", "Dataset Name:",
                                     placeholder = "e.g., My Query Dataset"),
                           conditionalPanel(
                               condition = "output.query_uploaded",
                               div(style = "margin-top: 10px; padding: 10px; background-color: #fcf8e3; border: 1px solid #faebcc; border-radius: 4px;",
                                   h5("✓ Query Dataset Loaded", style = "color: #8a6d3b; margin: 0;"),
                                   verbatimTextOutput("query_summary_upload")
                               )
                           )
                       )
                )
            )
        )
    ),
    
    fluidRow(
        box(
            title = "Cell Type Column Configuration", status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE,
            conditionalPanel(
                condition = "output.reference_uploaded || output.query_uploaded",
                h4("Configure Cell Type Columns"),
                p("Select the appropriate cell type columns for your uploaded datasets. These will be used across all diagnostic functions."),
                br(),
                
                fluidRow(
                    conditionalPanel(
                        condition = "output.reference_uploaded",
                        column(6,
                               wellPanel(
                                   h5("Reference Dataset Cell Type Column:"),
                                   selectInput("uploaded_ref_cell_type_col", "Select Column:",
                                               choices = NULL),
                                   verbatimTextOutput("ref_cell_types_preview")
                               )
                        )
                    ),
                    conditionalPanel(
                        condition = "output.query_uploaded",
                        column(6,
                               wellPanel(
                                   h5("Query Dataset Cell Type Column:"),
                                   selectInput("uploaded_query_cell_type_col", "Select Column:",
                                               choices = NULL),
                                   verbatimTextOutput("query_cell_types_preview")
                               )
                        )
                    )
                )
            )
        )
    ),
    
    fluidRow(
        box(
            title = "Data Upload Instructions", status = "warning", solidHeader = TRUE, width = 12, collapsible = TRUE, collapsed = TRUE,
            h4("How to Prepare Your Data"),
            tags$ol(
                tags$li("Save your SingleCellExperiment objects as .rds files using saveRDS()"),
                tags$li("Ensure your data contains:"),
                tags$ul(
                    tags$li("Expression assays (logcounts recommended)"),
                    tags$li("Cell metadata in colData with cell type annotations"),
                    tags$li("Reduced dimension embeddings (PCA, TSNE, UMAP) if using dimensional reduction plots")
                ),
                tags$li("Upload files and configure cell type columns"),
                tags$li("Your datasets will appear in the main dataset selection dropdown")
            ),
            br(),
            h4("Example R Code:"),
            tags$pre(
                "# Save your SingleCellExperiment object\n",
                "saveRDS(my_sce_object, 'my_reference_data.rds')\n",
                "saveRDS(my_query_object, 'my_query_data.rds')"
            ),
            br(),
            h4("Data Requirements:"),
            tags$ul(
                tags$li("Files must be in .rds format"),
                tags$li("Objects must be SingleCellExperiment class"),
                tags$li("Maximum file size: 100MB per file"),
                tags$li("Cell type columns should contain character/factor annotations")
            )
        )
    )
)