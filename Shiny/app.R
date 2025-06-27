# Complete library list for rsconnect deployment
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(BiocGenerics)
library(S4Vectors)
library(ggplot2)
library(GGally)
library(dplyr)
library(shinycssloaders)
library(ggridges)
library(MASS)
library(rlang)
library(igraph)
library(FNN)
library(isotree)
library(scales)
library(methods)
library(stats)
library(utils)
library(grDevices)
library(graphics)

# Source support functions
source("support/argumentCheck.R")
source("support/projectPCA.R")
source("support/generateColors.R")
source("support/plotCellTypePCA_app.R")
source("support/calculateGraphIntegration_app.R")
source("support/plotGraphIntegration_app.R")
source("support/detectAnomaly_app.R")
source("support/plotAnomaly_app.R")

# Load the datasets
query_cells <- readRDS("data/query_marrow_myeloid.rds")
reference_cells <- readRDS("data/reference_marrow_myeloid.rds")
reference_cells_subset <- readRDS("data/reference_subset_marrow_myeloid.rds")

# Define UI
ui <- dashboardPage(
    dashboardHeader(title = "scDiagnostics Interactive App"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Data Overview", tabName = "overview", icon = icon("table")),
            menuItem("PCA Projection", tabName = "pca", icon = icon("project-diagram")),
            menuItem("Graph Integration", tabName = "graph", icon = icon("network-wired")),
            menuItem("Anomaly Detection", tabName = "anomaly", icon = icon("exclamation-triangle")),
            menuItem("Cell Type Distribution", tabName = "distribution", icon = icon("chart-pie")),
            menuItem("Marker Gene Analysis", tabName = "markers", icon = icon("dna")),
            menuItem("Quality Metrics", tabName = "quality", icon = icon("chart-line"))
        )
    ),
    
    dashboardBody(
        tabItems(
            # Data Overview Tab
            tabItem(
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
            ),
            
            # PCA Projection Tab
            tabItem(
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
            ),
            
            # Graph Integration Tab
            tabItem(
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
            ),
            
            # Anomaly Detection Tab
            tabItem(
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
            ),
            
            # Placeholder tabs for future diagnostics
            tabItem(
                tabName = "distribution", 
                fluidRow(
                    box(
                        title = "Cell Type Distribution Analysis", status = "primary", solidHeader = TRUE, width = 12,
                        div(style = "text-align: center; padding: 100px;",
                            h2("Cell Type Distribution Analysis - Coming Soon", style = "color: #999;"),
                            p("This section will include cell type proportion analysis, distribution comparisons, and statistical tests.", style = "color: #666;")
                        )
                    )
                )
            ),
            tabItem(
                tabName = "markers", 
                fluidRow(
                    box(
                        title = "Marker Gene Analysis", status = "primary", solidHeader = TRUE, width = 12,
                        div(style = "text-align: center; padding: 100px;",
                            h2("Marker Gene Analysis - Coming Soon", style = "color: #999;"),
                            p("This section will include differential expression analysis, marker gene identification, and gene set enrichment.", style = "color: #666;")
                        )
                    )
                )
            ),
            tabItem(
                tabName = "quality", 
                fluidRow(
                    box(
                        title = "Quality Metrics", status = "primary", solidHeader = TRUE, width = 12,
                        div(style = "text-align: center; padding: 100px;",
                            h2("Quality Metrics - Coming Soon", style = "color: #999;"),
                            p("This section will include annotation confidence scores, quality control metrics, and batch effect analysis.", style = "color: #666;")
                        )
                    )
                )
            )
        )
    )
)

# Define Server
server <- function(input, output, session) {
    
    # Reactive values for datasets
    ref_data <- reactive({
        switch(input$reference_dataset,
               "reference_cells" = reference_cells,
               "reference_cells_subset" = reference_cells_subset)
    })
    
    query_data <- reactive({
        query_cells
    })
    
    # Update cell type column choices based on selected datasets (for PCA)
    observe({
        ref_cols <- colnames(colData(ref_data()))
        query_cols <- colnames(colData(query_data()))
        
        # Filter for likely cell type columns
        ref_cell_type_cols <- ref_cols
        query_cell_type_cols <- query_cols
        
        updateSelectInput(session, "ref_cell_type_col",
                          choices = ref_cell_type_cols,
                          selected = if(length(ref_cell_type_cols) > 0) ref_cell_type_cols[1] else NULL)
        
        updateSelectInput(session, "query_cell_type_col", 
                          choices = query_cell_type_cols,
                          selected = if(length(query_cell_type_cols) > 0) query_cell_type_cols[1] else NULL)
    })
    
    # Update cell type column choices for graph analysis
    observe({
        ref_cols <- colnames(colData(ref_data()))
        query_cols <- colnames(colData(query_data()))
        
        ref_cell_type_cols <- ref_cols
        query_cell_type_cols <- query_cols
        
        updateSelectInput(session, "graph_ref_cell_type_col",
                          choices = ref_cell_type_cols,
                          selected = if(length(ref_cell_type_cols) > 0) ref_cell_type_cols[1] else NULL)
        
        updateSelectInput(session, "graph_query_cell_type_col", 
                          choices = query_cell_type_cols,
                          selected = if(length(query_cell_type_cols) > 0) query_cell_type_cols[1] else NULL)
    })
    
    # Update cell type column choices for anomaly detection
    observe({
        ref_cols <- colnames(colData(ref_data()))
        query_cols <- colnames(colData(query_data()))
        
        ref_cell_type_cols <- ref_cols
        query_cell_type_cols <- query_cols
        
        updateSelectInput(session, "anomaly_ref_cell_type_col",
                          choices = ref_cell_type_cols,
                          selected = if(length(ref_cell_type_cols) > 0) ref_cell_type_cols[1] else NULL)
        
        updateSelectInput(session, "anomaly_query_cell_type_col", 
                          choices = query_cell_type_cols,
                          selected = if(length(query_cell_type_cols) > 0) query_cell_type_cols[1] else NULL)
    })
    
    # Update available cell types (for PCA)
    observe({
        req(input$ref_cell_type_col, input$query_cell_type_col)
        
        ref_types <- unique(colData(ref_data())[[input$ref_cell_type_col]])
        query_types <- unique(colData(query_data())[[input$query_cell_type_col]])
        all_types <- sort(unique(c(ref_types, query_types)))
        all_types <- all_types[!is.na(all_types)]
        
        updateSelectInput(session, "cell_types", choices = all_types)
    })
    
    # Update available cell types for graph analysis
    observe({
        req(input$graph_ref_cell_type_col, input$graph_query_cell_type_col)
        
        ref_types <- unique(colData(ref_data())[[input$graph_ref_cell_type_col]])
        query_types <- unique(colData(query_data())[[input$graph_query_cell_type_col]])
        all_types <- sort(unique(c(ref_types, query_types)))
        all_types <- all_types[!is.na(all_types)]
        
        updateSelectInput(session, "graph_cell_types", choices = all_types)
    })
    
    # Update available cell types for anomaly detection
    observe({
        req(input$anomaly_ref_cell_type_col, input$anomaly_query_cell_type_col)
        
        ref_types <- unique(colData(ref_data())[[input$anomaly_ref_cell_type_col]])
        query_types <- unique(colData(query_data())[[input$anomaly_query_cell_type_col]])
        all_types <- sort(unique(c(ref_types, query_types)))
        all_types <- all_types[!is.na(all_types)]
        
        updateSelectInput(session, "anomaly_cell_types", choices = all_types)
        updateSelectInput(session, "anomaly_cell_type_plot", choices = c("Combined", all_types), selected = "Combined")
    })
    
    # Dataset summaries with updated names
    output$ref_summary <- renderText({
        data <- ref_data()
        dataset_name <- switch(input$reference_dataset,
                               "reference_cells" = "Reference Marrow Myeloid Cells",
                               "reference_cells_subset" = "Reference Marrow Myeloid Cells (Promonocytes Removed)")
        paste(
            paste("Dataset:", dataset_name),
            paste("Dimensions:", paste(dim(data), collapse = " x ")),
            paste("Assays:", paste(assayNames(data), collapse = ", ")),
            paste("Reduced Dims:", paste(reducedDimNames(data), collapse = ", ")),
            paste("Cells:", ncol(data)),
            paste("Genes:", nrow(data)),
            sep = "\n"
        )
    })
    
    output$query_summary <- renderText({
        data <- query_data()
        paste(
            paste("Dataset: Query Marrow Myeloid Cells"),
            paste("Dimensions:", paste(dim(data), collapse = " x ")),
            paste("Assays:", paste(assayNames(data), collapse = ", ")),
            paste("Reduced Dims:", paste(reducedDimNames(data), collapse = ", ")),
            paste("Cells:", ncol(data)),
            paste("Genes:", nrow(data)),
            sep = "\n"
        )
    })
    
    output$ref_cell_type_cols <- renderText({
        cols <- colnames(colData(ref_data()))
        cell_type_cols <- cols[grepl("cell|type|label|cluster", cols, ignore.case = TRUE)]
        if(length(cell_type_cols) > 0) {
            paste(cell_type_cols, collapse = ", ")
        } else {
            "No obvious cell type columns found"
        }
    })
    
    output$query_cell_type_cols <- renderText({
        cols <- colnames(colData(query_data()))
        cell_type_cols <- cols[grepl("cell|type|label|cluster", cols, ignore.case = TRUE)]
        if(length(cell_type_cols) > 0) {
            paste(cell_type_cols, collapse = ", ")
        } else {
            "No obvious cell type columns found"
        }
    })
    
    # ============ PCA ANALYSIS SECTION ============
    
    # Reactive for PCA plot generation
    pca_plot_reactive <- eventReactive(input$generate_plot, {
        req(input$ref_cell_type_col, input$query_cell_type_col, input$pc_subset)
        
        # Convert PC subset to numeric
        pc_nums <- as.numeric(input$pc_subset)
        
        # Get selected cell types
        selected_cell_types <- input$cell_types
        if(is.null(selected_cell_types) || length(selected_cell_types) == 0) {
            selected_cell_types <- NULL  # Will use all cell types in function
        }
        
        # Call the adapted PCA function
        tryCatch({
            plotCellTypePCA_app(
                query_data = query_data(),
                reference_data = ref_data(),
                query_cell_type_col = input$query_cell_type_col,
                ref_cell_type_col = input$ref_cell_type_col,
                cell_types = selected_cell_types,
                pc_subset = pc_nums,
                assay_name = input$assay_name,
                lower_facet = input$lower_facet,
                diagonal_facet = input$diagonal_facet,
                upper_facet = input$upper_facet
            )
        }, error = function(e) {
            showNotification(paste("Error generating PCA plot:", e$message), type = "error")
            NULL
        })
    })
    
    # Render PCA plot
    output$pca_plot <- renderPlot({
        pca_plot_reactive()
    })
    
    # PCA Plot information
    output$plot_info <- renderText({
        if(input$generate_plot > 0) {
            pc_nums <- as.numeric(input$pc_subset)
            selected_cell_types <- input$cell_types
            if(is.null(selected_cell_types) || length(selected_cell_types) == 0) {
                ref_types <- unique(colData(ref_data())[[input$ref_cell_type_col]])
                query_types <- unique(colData(query_data())[[input$query_cell_type_col]])
                selected_cell_types <- sort(unique(c(ref_types, query_types)))
                selected_cell_types <- selected_cell_types[!is.na(selected_cell_types)]
            }
            
            dataset_name <- switch(input$reference_dataset,
                                   "reference_cells" = "Reference Marrow Myeloid Cells",
                                   "reference_cells_subset" = "Reference Marrow Myeloid Cells (Promonocytes Removed)")
            
            paste(
                paste("Reference Dataset:", dataset_name),
                paste("Query Dataset: Query Marrow Myeloid Cells"),
                paste("Selected PCs:", paste(pc_nums, collapse = ", ")),
                paste("Cell Types:", paste(selected_cell_types, collapse = ", ")),
                paste("Reference Column:", input$ref_cell_type_col),
                paste("Query Column:", input$query_cell_type_col),
                paste("Assay:", input$assay_name),
                paste("Lower Facet:", input$lower_facet),
                paste("Diagonal Facet:", input$diagonal_facet),
                paste("Upper Facet:", input$upper_facet),
                sep = "\n"
            )
        }
    })
    
    # Download handler for PCA plot
    output$download_pca_plot <- downloadHandler(
        filename = function() {
            paste0("PCA_projection_", Sys.Date(), ".png")
        },
        content = function(file) {
            plot_obj <- pca_plot_reactive()
            if(!is.null(plot_obj)) {
                ggsave(file, plot_obj, width = 12, height = 10, dpi = 300)
            }
        }
    )
    
    # ============ GRAPH INTEGRATION SECTION ============
    
    # Reactive for Graph Integration analysis
    graph_result_reactive <- eventReactive(input$generate_graph, {
        req(input$graph_ref_cell_type_col, input$graph_query_cell_type_col, input$graph_pc_subset)
        
        # Convert PC subset to numeric
        pc_nums <- as.numeric(input$graph_pc_subset)
        
        # Get selected cell types
        selected_cell_types <- input$graph_cell_types
        if(is.null(selected_cell_types) || length(selected_cell_types) == 0) {
            selected_cell_types <- NULL
        }
        
        # Call the graph integration function
        tryCatch({
            calculateGraphIntegration_app(
                query_data = query_data(),
                reference_data = ref_data(),
                query_cell_type_col = input$graph_query_cell_type_col,
                ref_cell_type_col = input$graph_ref_cell_type_col,
                cell_types = selected_cell_types,
                pc_subset = pc_nums,
                k_neighbors = input$k_neighbors,
                assay_name = input$graph_assay_name,
                resolution = input$resolution,
                min_cells_per_community = input$min_cells_per_community,
                high_query_prop_threshold = input$high_query_prop_threshold,
                cross_type_threshold = input$cross_type_threshold,
                local_consistency_threshold = input$local_consistency_threshold
            )
        }, error = function(e) {
            showNotification(paste("Error in graph analysis:", e$message), type = "error")
            NULL
        })
    })
    
    # Render Graph Integration plot
    output$graph_plot <- renderPlot({
        req(graph_result_reactive())
        
        result <- graph_result_reactive()
        
        tryCatch({
            plotGraphIntegration_app(
                x = result,
                plot_type = input$graph_plot_type,
                color_by = input$graph_color_by,
                exclude_reference_only = input$exclude_reference_only
            )
        }, error = function(e) {
            showNotification(paste("Error generating graph plot:", e$message), type = "error")
            NULL
        })
    })
    
    # Community details output
    output$community_details <- renderText({
        req(graph_result_reactive())
        
        result <- graph_result_reactive()
        community_comp <- result$community_composition
        
        if(nrow(community_comp) > 0) {
            paste(
                "=== COMMUNITY ANALYSIS RESULTS ===",
                paste("Total Communities:", nrow(community_comp)),
                paste("Query-Only Communities:", result$overall_metrics$high_query_prop_communities),
                paste("Cross-Type Communities:", result$overall_metrics$true_cross_type_communities),
                paste("Local Inconsistencies:", result$overall_metrics$total_locally_inconsistent_cells),
                "",
                "=== AFFECTED CELLS ===",
                paste("Query Cells in Query-Only Communities:", result$overall_metrics$total_high_query_prop_cells),
                paste("Query Cells in Cross-Type Communities:", result$overall_metrics$total_true_cross_mixing_cells),
                paste("Locally Inconsistent Query Cells:", result$overall_metrics$total_locally_inconsistent_cells),
                "",
                "=== RATES ===",
                paste("Mean Query Isolation Rate:", round(result$overall_metrics$mean_query_isolation_rate, 3)),
                paste("Mean Cross-Mixing Rate:", round(result$overall_metrics$mean_true_cross_mixing_rate, 3)),
                paste("Mean Local Inconsistency Rate:", round(result$overall_metrics$mean_local_inconsistency_rate, 3)),
                sep = "\n"
            )
        } else {
            "No community data available"
        }
    })
    
    # Graph analysis summary
    output$graph_analysis_summary <- renderText({
        req(graph_result_reactive())
        
        result <- graph_result_reactive()
        dataset_name <- switch(input$reference_dataset,
                               "reference_cells" = "Reference Marrow Myeloid Cells",
                               "reference_cells_subset" = "Reference Marrow Myeloid Cells (Promonocytes Removed)")
        
        paste(
            "=== ANALYSIS PARAMETERS ===",
            paste("Reference Dataset:", dataset_name),
            paste("Query Dataset: Query Marrow Myeloid Cells"),
            paste("Reference Cell Type Column:", input$graph_ref_cell_type_col),
            paste("Query Cell Type Column:", input$graph_query_cell_type_col),
            paste("Assay Used:", input$graph_assay_name),
            paste("Principal Components:", paste(input$graph_pc_subset, collapse = ", ")),
            "",
            "=== GRAPH PARAMETERS ===",
            paste("K-Neighbors:", input$k_neighbors),
            paste("Clustering Resolution:", input$resolution),
            paste("Min Cells per Community:", input$min_cells_per_community),
            paste("Graph Modularity:", round(result$overall_metrics$modularity, 3)),
            "",
            "=== THRESHOLDS ===",
            paste("Query-Only Threshold:", input$high_query_prop_threshold),
            paste("Cross-Type Threshold:", input$cross_type_threshold),
            paste("Local Consistency Threshold:", input$local_consistency_threshold),
            "",
            "=== CELL TYPES ANALYZED ===",
            paste(result$parameters$cell_types_analyzed, collapse = ", "),
            sep = "\n"
        )
    })
    
    # Download handler for graph plot
    output$download_graph_plot <- downloadHandler(
        filename = function() {
            paste0("Graph_integration_", input$graph_plot_type, "_", Sys.Date(), ".png")
        },
        content = function(file) {
            result <- graph_result_reactive()
            if(!is.null(result)) {
                plot_obj <- plotGraphIntegration_app(
                    x = result,
                    plot_type = input$graph_plot_type,
                    color_by = input$graph_color_by,
                    exclude_reference_only = input$exclude_reference_only
                )
                if(!is.null(plot_obj)) {
                    ggsave(file, plot_obj, width = 12, height = 10, dpi = 300)
                }
            }
        }
    )
    
    # ============ ANOMALY DETECTION SECTION ============
    
    # Reactive for Anomaly Detection analysis
    anomaly_result_reactive <- eventReactive(input$generate_anomaly, {
        req(input$anomaly_ref_cell_type_col, input$anomaly_query_cell_type_col, input$anomaly_pc_subset)
        
        # Convert PC subset to numeric
        pc_nums <- as.numeric(input$anomaly_pc_subset)
        
        # Get selected cell types
        selected_cell_types <- input$anomaly_cell_types
        if(is.null(selected_cell_types) || length(selected_cell_types) == 0) {
            selected_cell_types <- NULL
        }
        
        # Call the anomaly detection function
        tryCatch({
            detectAnomaly_app(
                reference_data = ref_data(),
                query_data = query_data(),
                ref_cell_type_col = input$anomaly_ref_cell_type_col,
                query_cell_type_col = input$anomaly_query_cell_type_col,
                cell_types = selected_cell_types,
                pc_subset = pc_nums,
                n_tree = input$n_tree,
                anomaly_threshold = input$anomaly_threshold,
                assay_name = input$anomaly_assay_name
            )
        }, error = function(e) {
            showNotification(paste("Error in anomaly detection:", e$message), type = "error")
            NULL
        })
    })
    
    # Render Anomaly Detection plot
    output$anomaly_plot <- renderPlot({
        req(anomaly_result_reactive())
        
        result <- anomaly_result_reactive()
        
        tryCatch({
            plotAnomaly_app(
                x = result,
                cell_type = input$anomaly_cell_type_plot,
                pc_subset = as.numeric(input$anomaly_pc_subset),
                data_type = input$anomaly_data_type,
                n_tree = input$n_tree,
                upper_facet = input$anomaly_upper_facet,
                diagonal_facet = input$anomaly_diagonal_facet
            )
        }, error = function(e) {
            showNotification(paste("Error generating anomaly plot:", e$message), type = "error")
            NULL
        })
    })
    
    # Anomaly statistics output
    output$anomaly_stats <- renderText({
        req(anomaly_result_reactive())
        
        result <- anomaly_result_reactive()
        cell_type <- input$anomaly_cell_type_plot
        data_type <- input$anomaly_data_type
        
        if(cell_type %in% names(result)) {
            cell_result <- result[[cell_type]]
            
            if(data_type == "query" && !is.null(cell_result[["query_anomaly_scores"]])) {
                anomaly_scores <- cell_result[["query_anomaly_scores"]]
                anomaly_flags <- cell_result[["query_anomaly"]]
                data_name <- "Query"
            } else if(data_type == "reference") {
                anomaly_scores <- cell_result[["reference_anomaly_scores"]]
                anomaly_flags <- cell_result[["reference_anomaly"]]
                data_name <- "Reference"
            } else {
                return("No data available for selected combination")
            }
            
            n_total <- length(anomaly_scores)
            n_anomalies <- sum(anomaly_flags)
            anomaly_rate <- n_anomalies / n_total
            mean_score <- mean(anomaly_scores, na.rm = TRUE)
            median_score <- median(anomaly_scores, na.rm = TRUE)
            
            paste(
                paste("=== ANOMALY STATISTICS ==="),
                paste("Cell Type:", cell_type),
                paste("Data Type:", data_name),
                paste("Anomaly Threshold:", input$anomaly_threshold),
                "",
                paste("=== RESULTS ==="),
                paste("Total Cells:", n_total),
                paste("Anomalous Cells:", n_anomalies),
                paste("Anomaly Rate:", sprintf("%.2f%%", anomaly_rate * 100)),
                "",
                paste("=== ANOMALY SCORES ==="),
                paste("Mean Score:", round(mean_score, 3)),
                paste("Median Score:", round(median_score, 3)),
                paste("Min Score:", round(min(anomaly_scores, na.rm = TRUE), 3)),
                paste("Max Score:", round(max(anomaly_scores, na.rm = TRUE), 3)),
                sep = "\n"
            )
        } else {
            "No results available for selected cell type"
        }
    })
    
    # Anomaly analysis summary
    output$anomaly_analysis_summary <- renderText({
        req(anomaly_result_reactive())
        
        result <- anomaly_result_reactive()
        dataset_name <- switch(input$reference_dataset,
                               "reference_cells" = "Reference Marrow Myeloid Cells",
                               "reference_cells_subset" = "Reference Marrow Myeloid Cells (Promonocytes Removed)")
        
        # Calculate overall statistics
        all_cell_types <- names(result)
        total_query_cells <- 0
        total_anomalous_query <- 0
        total_ref_cells <- 0
        total_anomalous_ref <- 0
        
        for(ct in all_cell_types) {
            if(!is.null(result[[ct]][["query_anomaly"]])) {
                total_query_cells <- total_query_cells + length(result[[ct]][["query_anomaly"]])
                total_anomalous_query <- total_anomalous_query + sum(result[[ct]][["query_anomaly"]])
            }
            if(!is.null(result[[ct]][["reference_anomaly"]])) {
                total_ref_cells <- total_ref_cells + length(result[[ct]][["reference_anomaly"]])
                total_anomalous_ref <- total_anomalous_ref + sum(result[[ct]][["reference_anomaly"]])
            }
        }
        
        paste(
            "=== ANALYSIS PARAMETERS ===",
            paste("Reference Dataset:", dataset_name),
            paste("Query Dataset: Query Marrow Myeloid Cells"),
            paste("Reference Cell Type Column:", input$anomaly_ref_cell_type_col),
            paste("Query Cell Type Column:", input$anomaly_query_cell_type_col),
            paste("Assay Used:", input$anomaly_assay_name),
            paste("Principal Components:", paste(input$anomaly_pc_subset, collapse = ", ")),
            "",
            "=== ISOLATION FOREST PARAMETERS ===",
            paste("Number of Trees:", input$n_tree),
            paste("Anomaly Threshold:", input$anomaly_threshold),
            "",
            "=== OVERALL RESULTS ===",
            paste("Cell Types Analyzed:", length(all_cell_types)),
            if(total_query_cells > 0) paste("Query Anomaly Rate:", sprintf("%.2f%%", (total_anomalous_query/total_query_cells) * 100)) else "",
            if(total_ref_cells > 0) paste("Reference Anomaly Rate:", sprintf("%.2f%%", (total_anomalous_ref/total_ref_cells) * 100)) else "",
            paste("Total Query Cells:", total_query_cells),
            paste("Total Anomalous Query Cells:", total_anomalous_query),
            "",
            "=== VISUALIZATION SETTINGS ===",
            paste("Plot Cell Type:", input$anomaly_cell_type_plot),
            paste("Plot Data Type:", input$anomaly_data_type),
            paste("Upper Facet:", input$anomaly_upper_facet),
            paste("Diagonal Facet:", input$anomaly_diagonal_facet),
            sep = "\n"
        )
    })
    
    # Download handler for anomaly plot
    output$download_anomaly_plot <- downloadHandler(
        filename = function() {
            paste0("Anomaly_detection_", input$anomaly_cell_type_plot, "_", input$anomaly_data_type, "_", Sys.Date(), ".png")
        },
        content = function(file) {
            result <- anomaly_result_reactive()
            if(!is.null(result)) {
                plot_obj <- plotAnomaly_app(
                    x = result,
                    cell_type = input$anomaly_cell_type_plot,
                    pc_subset = as.numeric(input$anomaly_pc_subset),
                    data_type = input$anomaly_data_type,
                    n_tree = input$n_tree,
                    upper_facet = input$anomaly_upper_facet,
                    diagonal_facet = input$anomaly_diagonal_facet
                )
                if(!is.null(plot_obj)) {
                    ggsave(file, plot_obj, width = 12, height = 10, dpi = 300)
                }
            }
        }
    )
}

# Run the application
shinyApp(ui = ui, server = server)
