# Gene expression dimensional reduction server functions
setup_dimred_outputs <- function(input, output, session, ref_data, query_data) {
    
    # Reactive for selected dataset - with proper error handling
    selected_dataset <- reactive({
        req(input$dimred_dataset)
        
        tryCatch({
            switch(input$dimred_dataset,
                   "query" = query_data(),
                   "reference_cells" = reference_cells,
                   "reference_cells_subset" = reference_cells_subset)
        }, error = function(e) {
            NULL
        })
    })
    
    # Get available genes for selected dataset - with proper error handling
    available_genes <- reactive({
        req(input$dimred_dataset)
        
        tryCatch({
            dataset <- selected_dataset()
            if(is.null(dataset)) {
                return(character(0))
            }
            
            genes <- rownames(dataset)
            sort(genes)  # Sort alphabetically
        }, error = function(e) {
            return(character(0))
        })
    })
    
    # Update gene choices when dataset changes
    observe({
        req(input$dimred_dataset)
        
        genes <- available_genes()
        if(length(genes) > 0) {
            updateSelectizeInput(session, "dimred_gene_name",
                                 choices = setNames(genes, genes),
                                 server = TRUE)
        }
    })
    
    # Update cell type column choices based on selected dataset
    observe({
        req(input$dimred_dataset)
        
        dataset <- selected_dataset()
        if(!is.null(dataset)) {
            cols <- colnames(colData(dataset))
            
            updateSelectInput(session, "dimred_cell_type_col",
                              choices = cols,
                              selected = if(length(cols) > 0) cols[1] else NULL)
        }
    })
    
    # Update available cell types for selected dataset
    observe({
        req(input$dimred_cell_type_col, input$dimred_dataset)
        
        dataset <- selected_dataset()
        if(!is.null(dataset) && input$dimred_cell_type_col %in% colnames(colData(dataset))) {
            cell_types <- unique(colData(dataset)[[input$dimred_cell_type_col]])
            cell_types <- sort(cell_types[!is.na(cell_types)])
            
            updateSelectInput(session, "dimred_cell_types", choices = cell_types)
        }
    })
    
    # Reactive for Dimensional Reduction plot
    dimred_plot_reactive <- eventReactive(input$generate_dimred, {
        req(input$dimred_cell_type_col, input$dimred_gene_name, input$dimred_method, input$dimred_dataset)
        
        # Validate gene name
        if(input$dimred_gene_name == "" || is.null(input$dimred_gene_name)) {
            showNotification("Please select a gene name", type = "warning")
            return(NULL)
        }
        
        # Get selected cell types (can be NULL for all)
        selected_cell_types <- input$dimred_cell_types
        if(is.null(selected_cell_types) || length(selected_cell_types) == 0) {
            selected_cell_types <- NULL
        }
        
        # Get PC subset for PCA
        pc_subset <- if(input$dimred_method == "PCA") {
            as.numeric(input$dimred_pc_subset)
        } else {
            1:5  # Default, not used for TSNE/UMAP
        }
        
        # Call the dimensional reduction function
        tryCatch({
            plotGeneExpressionDimred_app(
                se_object = selected_dataset(),
                method = input$dimred_method,
                pc_subset = pc_subset,
                feature = input$dimred_gene_name,
                cell_type_col = input$dimred_cell_type_col,
                cell_types = selected_cell_types,
                assay_name = input$dimred_assay_name
            )
        }, error = function(e) {
            showNotification(paste("Error generating dimensional reduction plot:", e$message), type = "error")
            NULL
        })
    })
    
    # Render Dimensional Reduction plot
    output$dimred_plot <- renderPlot({
        dimred_plot_reactive()
    })
    
    # Calculate expression statistics
    dimred_stats_reactive <- reactive({
        req(input$generate_dimred > 0, input$dimred_gene_name, input$dimred_dataset)
        
        tryCatch({
            dataset <- selected_dataset()
            if(is.null(dataset)) return(NULL)
            
            # Get gene expression data
            gene_expr <- assay(dataset, input$dimred_assay_name)[input$dimred_gene_name, ]
            
            # Get cell type specific expression if cell types selected
            if(!is.null(input$dimred_cell_types) && length(input$dimred_cell_types) > 0) {
                ct_indices <- which(colData(dataset)[[input$dimred_cell_type_col]] %in% input$dimred_cell_types)
                ct_expr <- gene_expr[ct_indices]
            } else {
                ct_expr <- gene_expr
            }
            
            list(
                overall = gene_expr,
                filtered = ct_expr,
                n_total = length(gene_expr),
                n_filtered = length(ct_expr)
            )
        }, error = function(e) {
            NULL
        })
    })
    
    # Dimensional reduction statistics output
    output$dimred_stats <- renderText({
        stats_data <- dimred_stats_reactive()
        if(is.null(stats_data)) return("No statistics available")
        
        # Calculate summary statistics for overall expression
        overall_stats <- c(
            mean = mean(stats_data$overall, na.rm = TRUE),
            median = median(stats_data$overall, na.rm = TRUE),
            sd = sd(stats_data$overall, na.rm = TRUE),
            min = min(stats_data$overall, na.rm = TRUE),
            max = max(stats_data$overall, na.rm = TRUE),
            zeros = sum(stats_data$overall == 0, na.rm = TRUE)
        )
        
        # Calculate summary statistics for filtered expression
        filtered_stats <- c(
            mean = mean(stats_data$filtered, na.rm = TRUE),
            median = median(stats_data$filtered, na.rm = TRUE),
            sd = sd(stats_data$filtered, na.rm = TRUE),
            min = min(stats_data$filtered, na.rm = TRUE),
            max = max(stats_data$filtered, na.rm = TRUE),
            zeros = sum(stats_data$filtered == 0, na.rm = TRUE)
        )
        
        dataset_name <- switch(input$dimred_dataset,
                               "query" = "Query Marrow Myeloid Cells",
                               "reference_cells" = "Reference Marrow Myeloid Cells",
                               "reference_cells_subset" = "Reference Marrow Myeloid Cells (Promonocytes Removed)")
        
        paste(
            paste("=== GENE EXPRESSION STATISTICS ==="),
            paste("Dataset:", dataset_name),
            paste("Gene:", input$dimred_gene_name),
            paste("Assay:", input$dimred_assay_name),
            "",
            paste("=== OVERALL DATASET ==="),
            paste("Total Cells:", stats_data$n_total),
            paste("Mean Expression:", round(overall_stats[["mean"]], 3)),
            paste("Median Expression:", round(overall_stats[["median"]], 3)),
            paste("SD:", round(overall_stats[["sd"]], 3)),
            paste("Range:", paste(round(overall_stats[["min"]], 3), "-", round(overall_stats[["max"]], 3))),
            paste("Zero Expression Cells:", overall_stats[["zeros"]], 
                  paste0("(", round(100 * overall_stats[["zeros"]] / stats_data$n_total, 1), "%)")),
            "",
            if(!is.null(input$dimred_cell_types) && length(input$dimred_cell_types) > 0) {
                paste(
                    paste("=== SELECTED CELL TYPES ==="),
                    paste("Cell Types:", paste(input$dimred_cell_types, collapse = ", ")),
                    paste("Filtered Cells:", stats_data$n_filtered),
                    paste("Mean Expression:", round(filtered_stats[["mean"]], 3)),
                    paste("Median Expression:", round(filtered_stats[["median"]], 3)),
                    paste("SD:", round(filtered_stats[["sd"]], 3)),
                    paste("Range:", paste(round(filtered_stats[["min"]], 3), "-", round(filtered_stats[["max"]], 3))),
                    paste("Zero Expression Cells:", filtered_stats[["zeros"]], 
                          paste0("(", round(100 * filtered_stats[["zeros"]] / stats_data$n_filtered, 1), "%)")),
                    sep = "\n"
                )
            } else {
                paste("=== CELL TYPE FILTERING ===", "All cell types included", sep = "\n")
            },
            sep = "\n"
        )
    })
    
    # Dimensional reduction analysis summary
    output$dimred_analysis_summary <- renderText({
        req(input$generate_dimred > 0, input$dimred_dataset)
        
        dataset_name <- switch(input$dimred_dataset,
                               "query" = "Query Marrow Myeloid Cells",
                               "reference_cells" = "Reference Marrow Myeloid Cells",
                               "reference_cells_subset" = "Reference Marrow Myeloid Cells (Promonocytes Removed)")
        
        method_desc <- switch(input$dimred_method,
                              "TSNE" = "t-SNE: t-distributed Stochastic Neighbor Embedding",
                              "UMAP" = "UMAP: Uniform Manifold Approximation and Projection",
                              "PCA" = "PCA: Principal Component Analysis")
        
        paste(
            "=== ANALYSIS PARAMETERS ===",
            paste("Dataset:", dataset_name),
            paste("Target Gene:", input$dimred_gene_name),
            paste("Cell Type Column:", input$dimred_cell_type_col),
            paste("Assay Used:", input$dimred_assay_name),
            if(!is.null(input$dimred_cell_types) && length(input$dimred_cell_types) > 0) {
                paste("Selected Cell Types:", paste(input$dimred_cell_types, collapse = ", "))
            } else {
                "Cell Type Filter: All cell types included"
            },
            "",
            "=== VISUALIZATION METHOD ===",
            paste("Method:", input$dimred_method),
            method_desc,
            if(input$dimred_method == "PCA") {
                paste("Principal Components:", paste(input$dimred_pc_subset, collapse = ", "))
            } else "",
            "",
            "=== PLOT INTERPRETATION ===",
            "??? Color intensity represents gene expression level",
            "??? Gray/light colors: Low or no expression",
            "??? Red/blue colors: High expression",
            "??? Spatial clustering may indicate cell type-specific expression",
            "",
            "=== ANALYSIS TIPS ===",
            "??? Look for expression patterns that correlate with cell clustering",
            "??? Compare expression across different cell types",
            "??? Consider trying different reduction methods for comparison",
            "??? Zero expression (gray dots) is common for many genes",
            if(input$dimred_method == "PCA") {
                "??? PCA shows multiple PC projections simultaneously"
            } else "",
            sep = "\n"
        )
    })
    
    # Download handler for dimensional reduction plot
    output$download_dimred_plot <- downloadHandler(
        filename = function() {
            dataset_short <- switch(input$dimred_dataset,
                                    "query" = "Query",
                                    "reference_cells" = "Ref",
                                    "reference_cells_subset" = "RefSub")
            paste0("DimRed_", input$dimred_method, "_", 
                   gsub("[^A-Za-z0-9]", "_", input$dimred_gene_name), "_", 
                   dataset_short, "_", Sys.Date(), ".png")
        },
        content = function(file) {
            plot_obj <- dimred_plot_reactive()
            if(!is.null(plot_obj)) {
                # Adjust size based on method
                if(input$dimred_method == "PCA") {
                    ggsave(file, plot_obj, width = 14, height = 12, dpi = 300)
                } else {
                    ggsave(file, plot_obj, width = 10, height = 8, dpi = 300)
                }
            }
        }
    )
}