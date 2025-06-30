# Graph integration server functions
setup_graph_outputs <- function(input, output, session, ref_data, query_data) {
    
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
}