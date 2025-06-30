# Anomaly detection server functions
setup_anomaly_outputs <- function(input, output, session, ref_data, query_data) {
    
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