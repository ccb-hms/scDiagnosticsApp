# App-adapted version of calculateDiscriminantSpace function
# This version includes all functionality from the original package function

calculateDiscriminantSpace_app <- function(reference_data,
                                           query_data = NULL,
                                           ref_cell_type_col,
                                           query_cell_type_col = NULL,
                                           cell_types = NULL,
                                           n_tree = 500,
                                           n_top = 20,
                                           eigen_threshold = 1e-1,
                                           calculate_metrics = FALSE,
                                           alpha = 0.01,
                                           assay_name = "logcounts") {
    
    # Load required libraries
    require(SingleCellExperiment)
    require(MASS)
    require(stats)
    
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
    
    if (!is.numeric(n_top) || n_top <= 0 || n_top != as.integer(n_top)) {
        stop("'n_top' must be a positive integer.")
    }
    
    if (!is.numeric(eigen_threshold) || eigen_threshold <= 0) {
        stop("'eigen_threshold' must be a positive number greater than 0.")
    }
    
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
        stop("'alpha' must be a positive number greater than 0 and less than 1.")
    }
    
    # Get common cell types if they are not specified by user
    if (is.null(cell_types)) {
        if (is.null(query_data)) {
            cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]])))
        } else {
            cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                           query_data[[query_cell_type_col]])))
        }
    }
    
    # Getting top variables
    var_imp <- calculateVarImpOverlap(reference_data = reference_data,
                                      ref_cell_type_col = ref_cell_type_col,
                                      n_tree = n_tree,
                                      n_top = n_top)
    
    # Get union of top n_top genes from all pairwise combinations
    all_top_genes <- unique(unlist(
        lapply(var_imp[["var_imp_ref"]], function(x) {
            x[["Gene"]][seq_len(n_top)]
        })))
    
    # Create a single discriminant model using all_top_genes
    # Extract reference matrix using top genes
    ref_mat <- t(as.matrix(assay(reference_data, assay_name)))[, all_top_genes]
    
    # Compute within-class and between-class scatter matrix
    sw <- sb <- matrix(0, length(all_top_genes), length(all_top_genes))
    overall_mean <- colMeans(ref_mat)
    
    for (cell_type in cell_types) {
        # Extract matrix for cell type
        ref_mat_subset <- ref_mat[which(reference_data[[ref_cell_type_col]] == cell_type), ]
        
        if (nrow(ref_mat_subset) == 0) {
            warning(paste("No reference cells found for cell type:", cell_type))
            next
        }
        
        # Ledoit-Wolf estimation for the current class
        lw_cov <- ledoitWolf_app(ref_mat_subset)
        # Update within-class scatter matrix
        sw <- sw + (nrow(ref_mat_subset) - 1) * lw_cov
        
        # Update between-class scatter matrix
        class_mean <- colMeans(ref_mat_subset)
        sb <- sb + nrow(ref_mat_subset) * (class_mean - overall_mean) %*%
            t(class_mean - overall_mean)
    }
    
    # Ensure symmetry of sw
    sw <- (sw + t(sw)) / 2
    
    # Check if sw is invertible
    if (det(sw) == 0 || any(is.na(sw))) {
        # Add regularization if sw is not invertible
        sw <- sw + diag(1e-6, nrow(sw))
    }
    
    # Solve generalized eigenvalue problem
    eig <- eigen(solve(sw) %*% sb)
    
    # Sort eigenvectors by eigenvalues
    discriminant_eigenvalues <- Re(eig[["values"]][which(Re(eig[["values"]]) > eigen_threshold), drop = FALSE])
    discriminant_eigenvectors <- Re(eig[["vectors"]][, which(Re(eig[["values"]]) > eigen_threshold), drop = FALSE])
    
    if (length(discriminant_eigenvalues) == 0) {
        stop("No discriminant vectors found above the eigen_threshold. Consider lowering the threshold.")
    }
    
    rownames(discriminant_eigenvectors) <- all_top_genes
    colnames(discriminant_eigenvectors) <- paste0("DV", seq_len(ncol(discriminant_eigenvectors)))
    
    # Compute projected data for reference
    ref_proj <- data.frame(
        ref_mat[which(reference_data[[ref_cell_type_col]] %in% cell_types), ] %*%
            discriminant_eigenvectors,
        reference_data[[ref_cell_type_col]][reference_data[[ref_cell_type_col]] %in% cell_types])
    colnames(ref_proj) <- c(paste0("DV", seq_len(ncol(discriminant_eigenvectors))), "cell_type")
    
    # Create discriminant_output
    discriminant_output <- list()
    discriminant_output[["discriminant_eigenvalues"]] <- discriminant_eigenvalues
    discriminant_output[["discriminant_eigenvectors"]] <- discriminant_eigenvectors
    discriminant_output[["ref_proj"]] <- ref_proj
    
    # Computations for query data
    if (!is.null(query_data)) {
        
        # Projection on discriminant space
        query_mat <- t(as.matrix(assay(query_data, assay_name)))[, all_top_genes]
        query_proj <- data.frame(
            query_mat[which(query_data[[query_cell_type_col]] %in% cell_types), ] %*% discriminant_eigenvectors,
            query_data[[query_cell_type_col]][query_data[[query_cell_type_col]] %in% cell_types])
        colnames(query_proj) <- c(paste0("DV", seq_len(ncol(discriminant_eigenvectors))), "cell_type")
        discriminant_output[["query_proj"]] <- query_proj
        
        if (calculate_metrics) {
            
            # Cosine similarity and Mahalanobis distance calculations
            cosine_similarity <- mahalanobis_dist <- numeric(nrow(query_proj))
            mahalanobis_crit <- numeric(length(cell_types))
            
            for (type_idx in seq_along(cell_types)) {
                type <- cell_types[type_idx]
                query_cells_of_type <- query_proj[query_proj[, "cell_type"] == type,
                                                  paste0("DV", seq_len(length(discriminant_eigenvalues))), drop = FALSE]
                ref_cells_of_type <- ref_proj[ref_proj[, "cell_type"] == type,
                                              paste0("DV", seq_len(length(discriminant_eigenvalues))), drop = FALSE]
                
                # Skip if we have no cells of this type
                if (nrow(query_cells_of_type) == 0 || nrow(ref_cells_of_type) == 0) {
                    next
                }
                
                # Calculate mean of reference projection for this cell type
                ref_mean <- colMeans(ref_cells_of_type)
                
                # Calculate covariance of reference projection for this cell type
                ref_cov <- cov(ref_cells_of_type)
                
                # Check if covariance matrix is invertible
                if (any(is.na(ref_cov)) || determinant(ref_cov)[["modulus"]][1] <= 0) {
                    # If not invertible, use a regularized version
                    ref_cov <- ledoitWolf_app(ref_cells_of_type)
                }
                
                # Calculate Mahalanobis distance
                mahalanobis_dist[query_proj[, "cell_type"] == type] <- mahalanobis(
                    query_cells_of_type, ref_mean, ref_cov)
                
                # Calculate cosine similarity
                cosine_similarity[query_proj[, "cell_type"] == type] <-
                    apply(query_cells_of_type, 1,
                          function(x, y) return(sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))),
                          y = ref_mean)
            }
            discriminant_output[["query_mahalanobis_dist"]] <- mahalanobis_dist
            discriminant_output[["mahalanobis_crit"]] <- qchisq(1 - alpha, df = length(cell_types))
            discriminant_output[["query_cosine_similarity"]] <- cosine_similarity
        }
    }
    
    # Set class and return
    class(discriminant_output) <- c(class(discriminant_output), "calculateDiscriminantSpaceObject")
    return(discriminant_output)
}

# Ledoit-Wolf covariance matrix estimation function
ledoitWolf_app <- function(class_data) {
    
    # Sample covariance matrix
    sample_cov <- cov(class_data)
    
    # Check for zero column means and replace them with a small constant if necessary
    col_means <- colMeans(class_data)
    col_means[col_means == 0] <- 1e-10
    
    # Calculate the shrinkage target (identity matrix scaled by the average variance)
    mean_variance <- mean(diag(sample_cov))
    shrinkage_target <- diag(mean_variance, ncol(class_data), ncol(class_data))
    
    # Calculate the shrinkage intensity
    phi_hat <- sum((class_data - col_means)^2) / (nrow(class_data) - 1)
    shrinkage_intensity <- (1 / nrow(class_data)) * min(phi_hat, mean_variance^2)
    
    # Ledoit-Wolf estimated covariance matrix
    lw_cov <- (1 - shrinkage_intensity) * sample_cov + shrinkage_intensity * shrinkage_target
    
    # Ensure symmetry
    lw_cov <- (lw_cov + t(lw_cov)) / 2
    
    # Return Ledoit-Wolf estimated covariance matrix
    return(lw_cov)
}