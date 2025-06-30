# Overview server functions
setup_overview_outputs <- function(input, output, session, ref_data, query_data) {
    
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
}