# PCA server functions
setup_pca_outputs <- function(input, output, session, ref_data, query_data) {
    
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
}