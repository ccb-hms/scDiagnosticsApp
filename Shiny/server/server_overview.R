# Overview server functions
setup_overview_outputs <- function(input, output, session, ref_data, query_data, uploaded_data) {
    
    # Dynamic reference dataset selection UI
    output$reference_dataset_selection <- renderUI({
        # Build choices list
        choices <- list(
            "Reference Marrow Myeloid Cells" = "reference_cells",
            "Reference Marrow Myeloid Cells (Promonocytes Removed)" = "reference_cells_subset"
        )
        
        # Add uploaded reference datasets if available
        if(!is.null(uploaded_data$uploaded_ref_name)) {
            choices[[uploaded_data$uploaded_ref_name]] <- uploaded_data$uploaded_ref_key
        }
        
        selectInput("reference_dataset", "Select Reference Dataset:",
                    choices = choices,
                    selected = names(choices)[1])
    })
    
    # Dynamic query dataset selection UI
    output$query_dataset_selection <- renderUI({
        # Query choices
        choices <- list("Query Marrow Myeloid Cells" = "query_cells")
        
        # Add uploaded query datasets if available
        if(!is.null(uploaded_data$uploaded_query_name)) {
            choices[[uploaded_data$uploaded_query_name]] <- uploaded_data$uploaded_query_key
        }
        
        selectInput("query_dataset", "Select Query Dataset:",
                    choices = choices,
                    selected = names(choices)[1])
    })
    
    # Helper function to get dataset display name
    get_dataset_display_name <- function(dataset_key, uploaded_data) {
        if(is.null(dataset_key)) return("Unknown Dataset")
        
        if(startsWith(dataset_key, "uploaded_ref_")) {
            return(uploaded_data$uploaded_ref_name %||% "Uploaded Reference Dataset")
        } else if(startsWith(dataset_key, "uploaded_query_")) {
            return(uploaded_data$uploaded_query_name %||% "Uploaded Query Dataset")
        } else {
            return(switch(dataset_key,
                          "reference_cells" = "Reference Marrow Myeloid Cells",
                          "reference_cells_subset" = "Reference Marrow Myeloid Cells (Promonocytes Removed)",
                          "query_cells" = "Query Marrow Myeloid Cells",
                          dataset_key))
        }
    }
    
    # Reference dataset summary
    output$ref_summary <- renderText({
        req(ref_data(), input$reference_dataset)
        data <- ref_data()
        
        # Get dataset name
        dataset_name <- get_dataset_display_name(input$reference_dataset, uploaded_data)
        
        # Get assay names, handling empty case
        assay_names <- assayNames(data)
        assay_text <- if(length(assay_names) > 0) paste(assay_names, collapse = ", ") else "None"
        
        # Get reduced dimension names, handling empty case
        reduced_dim_names <- reducedDimNames(data)
        reduced_dim_text <- if(length(reduced_dim_names) > 0) paste(reduced_dim_names, collapse = ", ") else "None"
        
        paste(
            paste("=== REFERENCE DATASET SUMMARY ==="),
            paste("Dataset:", dataset_name),
            paste("Dimensions:", paste(dim(data), collapse = " genes x ", " cells")),
            "",
            paste("=== DATA STRUCTURE ==="),
            paste("Total Cells:", format(ncol(data), big.mark = ",")),
            paste("Total Genes:", format(nrow(data), big.mark = ",")),
            paste("Available Assays:", assay_text),
            paste("Reduced Dimensions:", reduced_dim_text),
            "",
            paste("=== METADATA COLUMNS ==="),
            paste("ColData Columns:", ncol(colData(data))),
            paste("RowData Columns:", ncol(rowData(data))),
            "",
            paste("=== MEMORY USAGE ==="),
            paste("Object Size:", format(object.size(data), units = "MB")),
            sep = "\n"
        )
    })
    
    # Query dataset summary
    output$query_summary <- renderText({
        req(query_data(), input$query_dataset)
        data <- query_data()
        
        # Get dataset name
        dataset_name <- get_dataset_display_name(input$query_dataset, uploaded_data)
        
        # Get assay names, handling empty case
        assay_names <- assayNames(data)
        assay_text <- if(length(assay_names) > 0) paste(assay_names, collapse = ", ") else "None"
        
        # Get reduced dimension names, handling empty case
        reduced_dim_names <- reducedDimNames(data)
        reduced_dim_text <- if(length(reduced_dim_names) > 0) paste(reduced_dim_names, collapse = ", ") else "None"
        
        paste(
            paste("=== QUERY DATASET SUMMARY ==="),
            paste("Dataset:", dataset_name),
            paste("Dimensions:", paste(dim(data), collapse = " genes x ", " cells")),
            "",
            paste("=== DATA STRUCTURE ==="),
            paste("Total Cells:", format(ncol(data), big.mark = ",")),
            paste("Total Genes:", format(nrow(data), big.mark = ",")),
            paste("Available Assays:", assay_text),
            paste("Reduced Dimensions:", reduced_dim_text),
            "",
            paste("=== METADATA COLUMNS ==="),
            paste("ColData Columns:", ncol(colData(data))),
            paste("RowData Columns:", ncol(rowData(data))),
            "",
            paste("=== MEMORY USAGE ==="),
            paste("Object Size:", format(object.size(data), units = "MB")),
            sep = "\n"
        )
    })
    
    # Reference cell type columns
    output$ref_cell_type_cols <- renderText({
        req(ref_data())
        cols <- colnames(colData(ref_data()))
        
        if(length(cols) == 0) {
            return("No cell metadata columns found")
        }
        
        # Try to identify likely cell type columns
        cell_type_cols <- cols[grepl("cell|type|label|cluster|annotation|class", cols, ignore.case = TRUE)]
        
        if(length(cell_type_cols) > 0) {
            paste(
                paste("=== LIKELY CELL TYPE COLUMNS ==="),
                paste(cell_type_cols, collapse = ", "),
                "",
                paste("=== ALL AVAILABLE COLUMNS ==="),
                paste(cols, collapse = ", "),
                sep = "\n"
            )
        } else {
            paste(
                paste("=== ALL AVAILABLE COLUMNS ==="),
                paste("(No obvious cell type columns detected)"),
                paste(cols, collapse = ", "),
                sep = "\n"
            )
        }
    })
    
    # Query cell type columns
    output$query_cell_type_cols <- renderText({
        req(query_data())
        cols <- colnames(colData(query_data()))
        
        if(length(cols) == 0) {
            return("No cell metadata columns found")
        }
        
        # Try to identify likely cell type columns
        cell_type_cols <- cols[grepl("cell|type|label|cluster|annotation|class", cols, ignore.case = TRUE)]
        
        if(length(cell_type_cols) > 0) {
            paste(
                paste("=== LIKELY CELL TYPE COLUMNS ==="),
                paste(cell_type_cols, collapse = ", "),
                "",
                paste("=== ALL AVAILABLE COLUMNS ==="),
                paste(cols, collapse = ", "),
                sep = "\n"
            )
        } else {
            paste(
                paste("=== ALL AVAILABLE COLUMNS ==="),
                paste("(No obvious cell type columns detected)"),
                paste(cols, collapse = ", "),
                sep = "\n"
            )
        }
    })
    
    # Dataset compatibility check
    output$compatibility_check <- renderText({
        req(ref_data(), query_data())
        
        ref_data_obj <- ref_data()
        query_data_obj <- query_data()
        
        # Check gene overlap
        ref_genes <- rownames(ref_data_obj)
        query_genes <- rownames(query_data_obj)
        common_genes <- intersect(ref_genes, query_genes)
        gene_overlap_pct <- round(100 * length(common_genes) / length(union(ref_genes, query_genes)), 1)
        
        # Check assay compatibility
        ref_assays <- assayNames(ref_data_obj)
        query_assays <- assayNames(query_data_obj)
        common_assays <- intersect(ref_assays, query_assays)
        
        # Check reduced dimension compatibility
        ref_reddims <- reducedDimNames(ref_data_obj)
        query_reddims <- reducedDimNames(query_data_obj)
        common_reddims <- intersect(ref_reddims, query_reddims)
        
        # Generate compatibility report
        compatibility_status <- if(gene_overlap_pct > 80 && length(common_assays) > 0) {
            "??? COMPATIBLE"
        } else if(gene_overlap_pct > 50) {
            "??? PARTIALLY COMPATIBLE"
        } else {
            "??? INCOMPATIBLE"
        }
        
        paste(
            paste("=== DATASET COMPATIBILITY ANALYSIS ==="),
            paste("Overall Status:", compatibility_status),
            "",
            paste("=== GENE OVERLAP ==="),
            paste("Reference Genes:", format(length(ref_genes), big.mark = ",")),
            paste("Query Genes:", format(length(query_genes), big.mark = ",")),
            paste("Common Genes:", format(length(common_genes), big.mark = ",")),
            paste("Gene Overlap:", paste0(gene_overlap_pct, "%")),
            "",
            paste("=== ASSAY COMPATIBILITY ==="),
            paste("Reference Assays:", paste(ref_assays, collapse = ", ")),
            paste("Query Assays:", paste(query_assays, collapse = ", ")),
            paste("Common Assays:", if(length(common_assays) > 0) paste(common_assays, collapse = ", ") else "None"),
            "",
            paste("=== REDUCED DIMENSIONS ==="),
            paste("Reference RedDims:", if(length(ref_reddims) > 0) paste(ref_reddims, collapse = ", ") else "None"),
            paste("Query RedDims:", if(length(query_reddims) > 0) paste(query_reddims, collapse = ", ") else "None"),
            paste("Common RedDims:", if(length(common_reddims) > 0) paste(common_reddims, collapse = ", ") else "None"),
            "",
            paste("=== RECOMMENDATIONS ==="),
            if(gene_overlap_pct < 50) {
                "??? Low gene overlap may cause issues with analyses"
            } else "",
            if(length(common_assays) == 0) {
                "??? No common assays found - check data preprocessing"
            } else "",
            if(length(common_reddims) == 0) {
                "??? No common reduced dimensions - dimensional reduction plots may not work"
            } else "",
            if(gene_overlap_pct > 80 && length(common_assays) > 0) {
                "??? Datasets appear well-matched for comparative analysis"
            } else "",
            sep = "\n"
        )
    })
}