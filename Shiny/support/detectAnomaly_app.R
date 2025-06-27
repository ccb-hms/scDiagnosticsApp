# App-adapted version of detectAnomaly function
# This version includes all functionality from the original package function

detectAnomaly_app <- function(reference_data,
                              query_data = NULL,
                              ref_cell_type_col,
                              query_cell_type_col = NULL,
                              cell_types = NULL,
                              pc_subset = 1:5,
                              n_tree = 500,
                              anomaly_threshold = 0.6,
                              assay_name = "logcounts",
                              ...) {
    
    # Load required libraries
    require(isotree)
    require(SingleCellExperiment)
    
    # Basic input validation
    if (is.null(reference_data)) {
        stop("reference_data must be provided")
    }
    
    if (!ref_cell_type_col %in% colnames(colData(reference_data))) {
        stop("ref_cell_type_col not found in reference data")
    }
    
    if (!is.null(query_data) && !is.null(query_cell_type_col)) {
        if (!query_cell_type_col %in% colnames(colData(query_data))) {
            stop("query_cell_type_col not found in query data")
        }
    }
    
    if (!assay_name %in% assayNames(reference_data)) {
        stop("Specified assay not found in reference data")
    }
    
    if (!is.null(query_data) && !assay_name %in% assayNames(query_data)) {
        stop("Specified assay not found in query data")
    }
    
    # Parameter validation
    if (!is.numeric(n_tree) || n_tree <= 0 || n_tree != as.integer(n_tree)) {
        stop("'n_tree' must be a positive integer.")
    }
    
    if (!is.numeric(anomaly_threshold) || anomaly_threshold <= 0 || anomaly_threshold >= 1) {
        stop("'anomaly_threshold' must be a positive number greater than 0 and less than 1.")
    }
    
    # Extract reference and query labels
    reference_labels <- reference_data[[ref_cell_type_col]]
    query_labels <- NULL
    if (!is.null(query_data) && !is.null(query_cell_type_col)) {
        query_labels <- query_data[[query_cell_type_col]]
    }
    
    # Get common cell types if they are not specified by user
    if (is.null(cell_types)) {
        if (is.null(query_data)) {
            cell_types <- c(as.list(na.omit(unique(reference_labels))),
                            list(na.omit(unique(reference_labels))))
        } else {
            cell_types <- c(as.list(na.omit(intersect(unique(reference_labels),
                                                      unique(query_labels)))),
                            list(na.omit(intersect(unique(reference_labels),
                                                   unique(query_labels)))))
        }
    }
    
    # Get data from reference and query datasets
    if (!is.null(pc_subset)) {
        # Use PCA space
        if (!"PCA" %in% reducedDimNames(reference_data)) {
            stop("PCA not found in reference data. Please run PCA first or set pc_subset to NULL.")
        }
        
        reference_mat <- reducedDim(reference_data, "PCA")[, pc_subset, drop = FALSE]
        
        if (!is.null(query_data)) {
            pca_output <- projectPCA(query_data = query_data,
                                     reference_data = reference_data,
                                     query_cell_type_col = query_cell_type_col,
                                     ref_cell_type_col = ref_cell_type_col,
                                     pc_subset = pc_subset,
                                     assay_name = assay_name)
            query_mat <- pca_output[pca_output[["dataset"]] == "Query",
                                    paste0("PC", pc_subset), drop = FALSE]
        }
    } else {
        # Use assay data directly
        reference_mat <- t(as.matrix(assay(reference_data, assay_name)))
        if (!is.null(query_data)) {
            query_mat <- t(as.matrix(assay(query_data, assay_name)))
        }
    }
    
    # List to store output
    output <- list()
    
    for (cell_type in cell_types) {
        
        # Filter reference and query PCA data for the current cell type
        reference_mat_subset <- na.omit(reference_mat[which(reference_labels %in% cell_type), , drop = FALSE])
        
        if (nrow(reference_mat_subset) == 0) {
            warning(paste("No reference cells found for cell type:", paste(cell_type, collapse = ", ")))
            next
        }
        
        # Build isolation forest on reference PCA data for this cell type
        isolation_forest <- isotree::isolation.forest(reference_mat_subset,
                                                      ntree = n_tree, ...)
        
        # Calculate anomaly scores for reference data
        reference_anomaly_scores <- predict(isolation_forest,
                                            newdata = reference_mat_subset,
                                            type = "score")
        
        # Calculate anomaly scores for query data if available
        query_anomaly_scores <- NULL
        query_mat_subset <- NULL
        
        if (!is.null(query_data) && !is.null(query_labels)) {
            query_mat_subset <- na.omit(query_mat[which(query_labels %in% cell_type), , drop = FALSE])
            
            if (nrow(query_mat_subset) > 0) {
                query_anomaly_scores <- predict(isolation_forest,
                                                newdata = query_mat_subset,
                                                type = "score")
            }
        }
        
        # Store cell type anomaly scores and PCA data
        list_name <- ifelse(length(cell_type) == 1, cell_type, "Combined")
        output[[list_name]] <- list()
        output[[list_name]][["reference_anomaly_scores"]] <- reference_anomaly_scores
        output[[list_name]][["reference_anomaly"]] <- reference_anomaly_scores > anomaly_threshold
        output[[list_name]][["reference_mat_subset"]] <- reference_mat_subset
        
        if (!is.null(query_mat_subset) && nrow(query_mat_subset) > 0) {
            output[[list_name]][["query_mat_subset"]] <- query_mat_subset
            output[[list_name]][["query_anomaly_scores"]] <- query_anomaly_scores
            output[[list_name]][["query_anomaly"]] <- query_anomaly_scores > anomaly_threshold
        }
        
        if (!is.null(pc_subset)) {
            output[[list_name]][["var_explained"]] <- 
                attributes(reducedDim(reference_data, "PCA"))[["percentVar"]][pc_subset]
        }
    }
    
    # Set the class of the output
    class(output) <- c(class(output), "detectAnomalyObject")
    
    # Return anomaly, PCA data and optional PCA anomaly plots for each cell type
    return(output)
}