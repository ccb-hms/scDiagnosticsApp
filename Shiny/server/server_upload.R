# Data upload server functions
setup_upload_outputs <- function(input, output, session, uploaded_data) {
    
    # Reactive values to track upload status
    values <- reactiveValues(
        reference_uploaded = FALSE,
        query_uploaded = FALSE
    )
    
    # Handle reference data upload
    observeEvent(input$upload_reference, {
        req(input$upload_reference)
        
        tryCatch({
            # Load the RDS file
            ref_data <- readRDS(input$upload_reference$datapath)
            
            # Validate it's a SingleCellExperiment
            if (!is(ref_data, "SingleCellExperiment")) {
                showNotification("Uploaded file is not a SingleCellExperiment object", type = "error", duration = 5)
                return()
            }
            
            # Store in reactive values
            dataset_name <- if(input$reference_name != "") input$reference_name else "Uploaded Reference"
            uploaded_data[[paste0("uploaded_ref_", gsub("[^A-Za-z0-9]", "_", dataset_name))]] <- ref_data
            uploaded_data$uploaded_ref_name <- dataset_name
            uploaded_data$uploaded_ref_key <- paste0("uploaded_ref_", gsub("[^A-Za-z0-9]", "_", dataset_name))
            
            values$reference_uploaded <- TRUE
            
            # Update cell type column choices
            ref_cols <- colnames(colData(ref_data))
            updateSelectInput(session, "uploaded_ref_cell_type_col", choices = ref_cols)
            
            showNotification(paste("Reference dataset", dataset_name, "uploaded successfully!"), type = "message", duration = 5)
            
        }, error = function(e) {
            showNotification(paste("Error loading reference file:", e$message), type = "error", duration = 10)
            values$reference_uploaded <- FALSE
        })
    })
    
    # Handle query data upload
    observeEvent(input$upload_query, {
        req(input$upload_query)
        
        tryCatch({
            # Load the RDS file
            query_data <- readRDS(input$upload_query$datapath)
            
            # Validate it's a SingleCellExperiment
            if (!is(query_data, "SingleCellExperiment")) {
                showNotification("Uploaded file is not a SingleCellExperiment object", type = "error", duration = 5)
                return()
            }
            
            # Store in reactive values
            dataset_name <- if(input$query_name != "") input$query_name else "Uploaded Query"
            uploaded_data[[paste0("uploaded_query_", gsub("[^A-Za-z0-9]", "_", dataset_name))]] <- query_data
            uploaded_data$uploaded_query_name <- dataset_name
            uploaded_data$uploaded_query_key <- paste0("uploaded_query_", gsub("[^A-Za-z0-9]", "_", dataset_name))
            
            values$query_uploaded <- TRUE
            
            # Update cell type column choices
            query_cols <- colnames(colData(query_data))
            updateSelectInput(session, "uploaded_query_cell_type_col", choices = query_cols)
            
            showNotification(paste("Query dataset", dataset_name, "uploaded successfully!"), type = "message", duration = 5)
            
        }, error = function(e) {
            showNotification(paste("Error loading query file:", e$message), type = "error", duration = 10)
            values$query_uploaded <- FALSE
        })
    })
    
    # Output upload status
    output$reference_uploaded <- reactive({
        values$reference_uploaded
    })
    outputOptions(output, "reference_uploaded", suspendWhenHidden = FALSE)
    
    output$query_uploaded <- reactive({
        values$query_uploaded
    })
    outputOptions(output, "query_uploaded", suspendWhenHidden = FALSE)
    
    # Reference dataset summary
    output$reference_summary_upload <- renderText({
        if(values$reference_uploaded && !is.null(uploaded_data$uploaded_ref_key)) {
            ref_data <- uploaded_data[[uploaded_data$uploaded_ref_key]]
            paste(
                paste("Dimensions:", paste(dim(ref_data), collapse = " x ")),
                paste("Assays:", paste(assayNames(ref_data), collapse = ", ")),
                paste("Reduced Dims:", paste(reducedDimNames(ref_data), collapse = ", ")),
                paste("Cells:", ncol(ref_data)),
                paste("Genes:", nrow(ref_data)),
                sep = "\n"
            )
        }
    })
    
    # Query dataset summary
    output$query_summary_upload <- renderText({
        if(values$query_uploaded && !is.null(uploaded_data$uploaded_query_key)) {
            query_data <- uploaded_data[[uploaded_data$uploaded_query_key]]
            paste(
                paste("Dimensions:", paste(dim(query_data), collapse = " x ")),
                paste("Assays:", paste(assayNames(query_data), collapse = ", ")),
                paste("Reduced Dims:", paste(reducedDimNames(query_data), collapse = ", ")),
                paste("Cells:", ncol(query_data)),
                paste("Genes:", nrow(query_data)),
                sep = "\n"
            )
        }
    })
    
    # Cell type previews
    output$ref_cell_types_preview <- renderText({
        if(values$reference_uploaded && !is.null(input$uploaded_ref_cell_type_col) && !is.null(uploaded_data$uploaded_ref_key)) {
            ref_data <- uploaded_data[[uploaded_data$uploaded_ref_key]]
            if(input$uploaded_ref_cell_type_col %in% colnames(colData(ref_data))) {
                cell_types <- unique(colData(ref_data)[[input$uploaded_ref_cell_type_col]])
                cell_types <- cell_types[!is.na(cell_types)]
                paste("Cell Types Found:", paste(head(cell_types, 10), collapse = ", "),
                      if(length(cell_types) > 10) paste("... and", length(cell_types) - 10, "more") else "")
            }
        }
    })
    
    output$query_cell_types_preview <- renderText({
        if(values$query_uploaded && !is.null(input$uploaded_query_cell_type_col) && !is.null(uploaded_data$uploaded_query_key)) {
            query_data <- uploaded_data[[uploaded_data$uploaded_query_key]]
            if(input$uploaded_query_cell_type_col %in% colnames(colData(query_data))) {
                cell_types <- unique(colData(query_data)[[input$uploaded_query_cell_type_col]])
                cell_types <- cell_types[!is.na(cell_types)]
                paste("Cell Types Found:", paste(head(cell_types, 10), collapse = ", "),
                      if(length(cell_types) > 10) paste("... and", length(cell_types) - 10, "more") else "")
            }
        }
    })
    
    # Store cell type column selections
    observe({
        if(values$reference_uploaded && !is.null(input$uploaded_ref_cell_type_col)) {
            uploaded_data$uploaded_ref_cell_type_col <- input$uploaded_ref_cell_type_col
        }
    })
    
    observe({
        if(values$query_uploaded && !is.null(input$uploaded_query_cell_type_col)) {
            uploaded_data$uploaded_query_cell_type_col <- input$uploaded_query_cell_type_col
        }
    })
}