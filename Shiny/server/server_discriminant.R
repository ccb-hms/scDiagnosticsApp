# Discriminant analysis server functions
setup_discriminant_outputs <- function(input, output, session, ref_data, query_data) {
    
    # Reactive for Discriminant Space analysis
    discriminant_result_reactive <- eventReactive(input$generate_discriminant, {
        req(input$disc_ref_cell_type_col, input$disc_query_cell_type_col)
        
        # Get selected cell types
        selected_cell_types <- input$disc_cell_types
        if(is.null(selected_cell_types) || length(selected_cell_types) == 0) {
            selected_cell_types <- NULL
        }
        
        # Call the discriminant space function
        tryCatch({
            calculateDiscriminantSpace_app(
                reference_data = ref_data(),
                query_data = query_data(),
                ref_cell_type_col = input$disc_ref_cell_type_col,
                query_cell_type_col = input$disc_query_cell_type_col,
                cell_types = selected_cell_types,
                n_tree = input$disc_n_tree,
                n_top = input$disc_n_top,
                eigen_threshold = input$eigen_threshold,
                calculate_metrics = input$calculate_metrics,
                alpha = input$alpha,
                assay_name = input$disc_assay_name
            )
        }, error = function(e) {
            showNotification(paste("Error in discriminant analysis:", e$message), type = "error")
            NULL
        })
    })
    
    # Update discriminant vector choices based on results
    observe({
        result <- discriminant_result_reactive()
        if(!is.null(result)) {
            n_dvs <- length(result$discriminant_eigenvalues)
            updateCheckboxGroupInput(session, "dv_subset",
                                     choices = setNames(1:n_dvs, paste0("DV", 1:n_dvs)),
                                     selected = 1:min(3, n_dvs))
        }
    })
    
    # Render Discriminant Space plot
    output$discriminant_plot <- renderPlot({
        req(discriminant_result_reactive())
        
        result <- discriminant_result_reactive()
        
        tryCatch({
            plotDiscriminantSpace_app(
                x = result,
                cell_types = input$disc_cell_types,
                dv_subset = as.numeric(input$dv_subset),
                lower_facet = input$disc_lower_facet,
                diagonal_facet = input$disc_diagonal_facet,
                upper_facet = input$disc_upper_facet
            )
        }, error = function(e) {
            showNotification(paste("Error generating discriminant plot:", e$message), type = "error")
            NULL
        })
    })
    
    # Discriminant statistics output
    output$discriminant_stats <- renderText({
        req(discriminant_result_reactive())
        
        result <- discriminant_result_reactive()
        
        eigenvalues <- result$discriminant_eigenvalues
        n_dvs <- length(eigenvalues)
        
        # Calculate cumulative variance explained
        total_variance <- sum(eigenvalues)
        cumulative_variance <- cumsum(eigenvalues) / total_variance
        
        paste(
            "=== DISCRIMINANT VECTORS & EIGENVALUES ===",
            paste("Number of Discriminant Vectors:", n_dvs),
            paste("Total Discriminant Variance:", round(total_variance, 4)),
            "",
            "=== EIGENVALUES ===",
            paste(paste0("DV", 1:n_dvs, ": ", round(eigenvalues, 4)), collapse = "\n"),
            "",
            "=== VARIANCE EXPLAINED (Cumulative) ===",
            paste(paste0("DV", 1:n_dvs, ": ", sprintf("%.2f%%", cumulative_variance * 100)), collapse = "\n"),
            "",
            if(input$calculate_metrics && !is.null(result$query_cosine_similarity)) {
                paste(
                    "=== SIMILARITY METRICS ===",
                    paste("Mean Cosine Similarity:", round(mean(result$query_cosine_similarity, na.rm = TRUE), 3)),
                    paste("Mean Mahalanobis Distance:", round(mean(result$query_mahalanobis_dist, na.rm = TRUE), 3)),
                    paste("Mahalanobis Critical Value:", round(result$mahalanobis_crit, 3)),
                    paste("Cells Above Mahalanobis Threshold:", sum(result$query_mahalanobis_dist > result$mahalanobis_crit, na.rm = TRUE)),
                    sep = "\n"
                )
            } else "",
            sep = "\n"
        )
    })
    
    # Discriminant analysis summary
    output$discriminant_analysis_summary <- renderText({
        req(discriminant_result_reactive())
        
        result <- discriminant_result_reactive()
        dataset_name <- switch(input$reference_dataset,
                               "reference_cells" = "Reference Marrow Myeloid Cells",
                               "reference_cells_subset" = "Reference Marrow Myeloid Cells (Promonocytes Removed)")
        
        # Get analyzed cell types
        analyzed_cell_types <- unique(c(result$ref_proj$cell_type, result$query_proj$cell_type))
        
        paste(
            "=== ANALYSIS PARAMETERS ===",
            paste("Reference Dataset:", dataset_name),
            paste("Query Dataset: Query Marrow Myeloid Cells"),
            paste("Reference Cell Type Column:", input$disc_ref_cell_type_col),
            paste("Query Cell Type Column:", input$disc_query_cell_type_col),
            paste("Assay Used:", input$disc_assay_name),
            "",
            "=== RANDOM FOREST PARAMETERS ===",
            paste("Number of Trees:", input$disc_n_tree),
            paste("Top Variables per Comparison:", input$disc_n_top),
            "",
            "=== DISCRIMINANT ANALYSIS PARAMETERS ===",
            paste("Eigenvalue Threshold:", input$eigen_threshold),
            paste("Calculate Similarity Metrics:", input$calculate_metrics),
            if(input$calculate_metrics) paste("Significance Level (α):", input$alpha) else "",
            "",
            "=== DATA SUMMARY ===",
            paste("Cell Types Analyzed:", paste(sort(analyzed_cell_types), collapse = ", ")),
            paste("Reference Cells:", nrow(result$ref_proj)),
            paste("Query Cells:", nrow(result$query_proj)),
            paste("Discriminant Vectors Generated:", length(result$discriminant_eigenvalues)),
            "",
            "=== VISUALIZATION SETTINGS ===",
            paste("Selected Discriminant Vectors:", paste(input$dv_subset, collapse = ", ")),
            paste("Lower Facet:", input$disc_lower_facet),
            paste("Diagonal Facet:", input$disc_diagonal_facet),
            paste("Upper Facet:", input$disc_upper_facet),
            sep = "\n"
        )
    })
    
    # Download handler for discriminant plot
    output$download_discriminant_plot <- downloadHandler(
        filename = function() {
            paste0("Discriminant_space_", paste(input$dv_subset, collapse = "_"), "_", Sys.Date(), ".png")
        },
        content = function(file) {
            result <- discriminant_result_reactive()
            if(!is.null(result)) {
                plot_obj <- plotDiscriminantSpace_app(
                    x = result,
                    cell_types = input$disc_cell_types,
                    dv_subset = as.numeric(input$dv_subset),
                    lower_facet = input$disc_lower_facet,
                    diagonal_facet = input$disc_diagonal_facet,
                    upper_facet = input$disc_upper_facet
                )
                if(!is.null(plot_obj)) {
                    ggsave(file, plot_obj, width = 12, height = 10, dpi = 300)
                }
            }
        }
    )
}