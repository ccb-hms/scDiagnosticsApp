# Wasserstein distance server functions
setup_wasserstein_outputs <- function(input, output, session, ref_data, query_data) {
    
    # Reactive for Wasserstein Distance analysis
    wasserstein_result_reactive <- eventReactive(input$generate_wasserstein, {
        req(input$wass_ref_cell_type_col, input$wass_query_cell_type_col, input$wass_pc_subset, input$n_resamples)
        
        # Convert PC subset to numeric
        pc_nums <- as.numeric(input$wass_pc_subset)
        
        # Convert n_resamples to numeric (fix for the error)
        n_resamples_num <- as.numeric(input$n_resamples)
        
        # Validate n_resamples
        if(is.na(n_resamples_num) || n_resamples_num < 1) {
            showNotification("Number of resamples must be a positive number", type = "error")
            return(NULL)
        }
        
        # Get selected cell types
        selected_cell_types <- input$wass_cell_types
        if(is.null(selected_cell_types) || length(selected_cell_types) == 0) {
            selected_cell_types <- NULL
        }
        
        # Call the Wasserstein distance function
        tryCatch({
            calculateWassersteinDistance_app(
                query_data = query_data(),
                reference_data = ref_data(),
                ref_cell_type_col = input$wass_ref_cell_type_col,
                query_cell_type_col = input$wass_query_cell_type_col,
                cell_types = selected_cell_types,
                pc_subset = pc_nums,
                n_resamples = n_resamples_num,  # Use the converted numeric value
                assay_name = input$wass_assay_name
            )
        }, error = function(e) {
            showNotification(paste("Error in Wasserstein analysis:", e$message), type = "error")
            NULL
        })
    })
    
    # Update plot cell type choices based on results
    observe({
        result <- wasserstein_result_reactive()
        if(!is.null(result)) {
            updateSelectInput(session, "wass_plot_cell_types",
                              choices = result$cell_types,
                              selected = result$cell_types)
        }
    })
    
    # Render Wasserstein Distance plot
    output$wasserstein_plot <- renderPlot({
        req(wasserstein_result_reactive())
        
        result <- wasserstein_result_reactive()
        
        # Get selected cell types for plotting
        plot_cell_types <- input$wass_plot_cell_types
        if(is.null(plot_cell_types) || length(plot_cell_types) == 0) {
            plot_cell_types <- result$cell_types
        }
        
        tryCatch({
            plotWassersteinDistance_app(
                x = result,
                cell_types = plot_cell_types
            )
        }, error = function(e) {
            showNotification(paste("Error generating Wasserstein plot:", e$message), type = "error")
            NULL
        })
    })
    
    # Wasserstein statistics output
    output$wasserstein_stats <- renderText({
        req(wasserstein_result_reactive())
        
        result <- wasserstein_result_reactive()
        prob_sup <- result$probability_superiority
        
        # Calculate summary statistics
        mean_prob <- mean(prob_sup, na.rm = TRUE)
        median_prob <- median(prob_sup, na.rm = TRUE)
        min_prob <- min(prob_sup, na.rm = TRUE)
        max_prob <- max(prob_sup, na.rm = TRUE)
        
        # Identify problematic cell types (high probability of superiority)
        high_prob_types <- names(prob_sup)[prob_sup > 0.7]
        moderate_prob_types <- names(prob_sup)[prob_sup > 0.5 & prob_sup <= 0.7]
        
        paste(
            "=== PROBABILITY OF SUPERIORITY STATISTICS ===",
            paste("Number of Cell Types Analyzed:", length(prob_sup)),
            paste("Mean Probability:", round(mean_prob, 3)),
            paste("Median Probability:", round(median_prob, 3)),
            paste("Range:", paste(round(min_prob, 3), "-", round(max_prob, 3))),
            "",
            "=== INDIVIDUAL CELL TYPE RESULTS ===",
            paste(paste0(names(prob_sup), ": ", sprintf("%.3f", prob_sup)), collapse = "\n"),
            "",
            "=== INTERPRETATION ===",
            if(length(high_prob_types) > 0) {
                paste("High Divergence (P > 0.7):", paste(high_prob_types, collapse = ", "))
            } else "",
            if(length(moderate_prob_types) > 0) {
                paste("Moderate Divergence (0.5 < P ??? 0.7):", paste(moderate_prob_types, collapse = ", "))
            } else "",
            paste("Low Divergence (P ??? 0.5):", 
                  paste(names(prob_sup)[prob_sup <= 0.5], collapse = ", ")),
            sep = "\n"
        )
    })
    
    # Wasserstein analysis summary
    output$wasserstein_analysis_summary <- renderText({
        req(wasserstein_result_reactive())
        
        result <- wasserstein_result_reactive()
        dataset_name <- switch(input$reference_dataset,
                               "reference_cells" = "Reference Marrow Myeloid Cells",
                               "reference_cells_subset" = "Reference Marrow Myeloid Cells (Promonocytes Removed)")
        
        paste(
            "=== ANALYSIS PARAMETERS ===",
            paste("Reference Dataset:", dataset_name),
            paste("Query Dataset: Query Marrow Myeloid Cells"),
            paste("Reference Cell Type Column:", input$wass_ref_cell_type_col),
            paste("Query Cell Type Column:", input$wass_query_cell_type_col),
            paste("Assay Used:", input$wass_assay_name),
            paste("Principal Components:", paste(input$wass_pc_subset, collapse = ", ")),
            paste("Number of Resamples:", as.numeric(input$n_resamples)),  # Convert here too
            "",
            "=== WASSERSTEIN DISTANCE ANALYSIS ===",
            paste("Cell Types Analyzed:", paste(result$cell_types, collapse = ", ")),
            paste("Total Distributions Computed:", length(result$cell_types) * 2),
            paste("Distance Comparisons per Cell Type:", as.numeric(input$n_resamples)^2),
            "",
            "=== METHODOLOGY ===",
            "Reference-Reference: Distances between random reference subsamples",
            "Reference-Query: Distances between reference and query samples", 
            "Probability of Superiority: P(Ref-Query distance > Ref-Ref distance)",
            "",
            "=== INTERPRETATION GUIDE ===",
            "P ??? 0.5: Query cells similar to reference distribution",
            "P > 0.7: Query cells significantly divergent from reference",
            "P < 0.3: Query cells more similar than expected (unusual)",
            sep = "\n"
        )
    })
    
    # Download handler for Wasserstein plot
    output$download_wasserstein_plot <- downloadHandler(
        filename = function() {
            paste0("Wasserstein_distance_", Sys.Date(), ".png")
        },
        content = function(file) {
            result <- wasserstein_result_reactive()
            if(!is.null(result)) {
                plot_cell_types <- input$wass_plot_cell_types
                if(is.null(plot_cell_types) || length(plot_cell_types) == 0) {
                    plot_cell_types <- result$cell_types
                }
                
                plot_obj <- plotWassersteinDistance_app(
                    x = result,
                    cell_types = plot_cell_types
                )
                if(!is.null(plot_obj)) {
                    ggsave(file, plot_obj, width = 14, height = 10, dpi = 300)
                }
            }
        }
    )
}