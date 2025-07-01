# App version of plotMarkerExpression function (UPDATE the subtitle line)
plotMarkerExpression_app <- function(reference_data,
                                     query_data,
                                     ref_cell_type_col,
                                     query_cell_type_col,
                                     cell_type,
                                     gene_name,
                                     assay_name = "logcounts",
                                     normalization = c("z_score", "min_max", "rank", "none")) {
    
    # Match normalization argument
    normalization <- match.arg(normalization)
    
    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_type,
                  assay_name = assay_name)
    
    # Get common cell types if they are not specified by user
    if(is.null(cell_type)){
        cell_type <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                      query_data[[query_cell_type_col]])))
    }
    
    # Check if gene_name is present in both query_data and reference_data
    if (!(gene_name %in% rownames(assay(query_data)) &&
          gene_name %in% rownames(assay(reference_data)))) {
        stop("gene_name: \'",
             gene_name,
             "\' is not present in the row names of both \'query_data\' and \'reference_data\'.")
    }
    
    # Get expression of the specified gene for reference and query datasets
    ref_gene_expression <- assay(reference_data, assay_name)[gene_name, ]
    query_gene_expression <- assay(query_data, assay_name)[gene_name, ]
    ref_gene_expression_specific <- assay(
        reference_data, assay_name)[gene_name,
                                    which(reference_data[[ref_cell_type_col]] %in%
                                              cell_type)]
    query_gene_expression_specific <- assay(
        query_data, assay_name)[gene_name,
                                which(query_data[[query_cell_type_col]] %in%
                                          cell_type)]
    
    # Transformation functions
    .quantileTransformation <- function(x) {
        return(100 * rank(x, ties.method = "average") / length(x))
    }
    
    .zScoreTransformation <- function(x) {
        if(sd(x) == 0) return(rep(0, length(x)))
        return((x - mean(x)) / sd(x))
    }
    
    .minMaxTransformation <- function(x) {
        if(max(x) == min(x)) return(rep(0.5, length(x)))
        return((x - min(x)) / (max(x) - min(x)))
    }
    
    # Apply selected normalization
    if (normalization == "rank") {
        ref_gene_expression_norm <-
            .quantileTransformation(ref_gene_expression)
        query_gene_expression_norm <-
            .quantileTransformation(query_gene_expression)
        ref_gene_expression_specific_norm <-
            .quantileTransformation(ref_gene_expression_specific)
        query_gene_expression_specific_norm <-
            .quantileTransformation(query_gene_expression_specific)
        x_label <-
            paste("Quantile Rank Normalized Gene Expression:", gene_name)
        
    } else if (normalization == "z_score") {
        ref_gene_expression_norm <-
            .zScoreTransformation(ref_gene_expression)
        query_gene_expression_norm <-
            .zScoreTransformation(query_gene_expression)
        ref_gene_expression_specific_norm <-
            .zScoreTransformation(ref_gene_expression_specific)
        query_gene_expression_specific_norm <-
            .zScoreTransformation(query_gene_expression_specific)
        x_label <-
            paste("Z-Score Normalized Gene Expression:", gene_name)
        
    } else if (normalization == "min_max") {
        ref_gene_expression_norm <-
            .minMaxTransformation(ref_gene_expression)
        query_gene_expression_norm <-
            .minMaxTransformation(query_gene_expression)
        ref_gene_expression_specific_norm <-
            .minMaxTransformation(ref_gene_expression_specific)
        query_gene_expression_specific_norm <-
            .minMaxTransformation(query_gene_expression_specific)
        x_label <- paste("Min-Max Normalized Gene Expression:", gene_name)
    } else if (normalization == "none") {
        ref_gene_expression_norm <-
            ref_gene_expression
        query_gene_expression_norm <-
            query_gene_expression
        ref_gene_expression_specific_norm <-
            ref_gene_expression_specific
        query_gene_expression_specific_norm <-
            query_gene_expression_specific
        x_label <-
            paste("Log-Normalized Gene Expression:", gene_name)
        
    }
    
    # Create a combined vector of gene expression values
    combined_gene_expression <- c(
        ref_gene_expression_norm,
        query_gene_expression_norm,
        ref_gene_expression_specific_norm,
        query_gene_expression_specific_norm)
    
    # Create a grouping vector for dataset types
    dataset_types <- rep(c("Reference", "Query", "Reference", "Query"),
                         times = c(
                             length(ref_gene_expression),
                             length(query_gene_expression),
                             length(ref_gene_expression_specific),
                             length(query_gene_expression_specific)))
    
    # Combine the gene expression values and dataset types
    data <- data.frame(
        GeneExpression = combined_gene_expression,
        Dataset = dataset_types,
        plot_type = rep(c("Overall Distribution", "Cell Type-Specific Distribution"),
                        times = c(length(ref_gene_expression) +
                                      length(query_gene_expression),
                                  length(ref_gene_expression_specific) +
                                      length(query_gene_expression_specific))))
    
    # Create a stacked density plot
    plot_obj <- ggplot2::ggplot(
        data,
        ggplot2::aes(x = .data[["GeneExpression"]],
                     fill = .data[["Dataset"]])) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::facet_wrap(~ .data[["plot_type"]], scales = "free") +
        ggplot2::labs(title = paste("Gene Expression Distribution Analysis:", gene_name),
                      subtitle = paste("Cell Type:", cell_type),  # Updated for single cell type
                      x = x_label,
                      y = "Density") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            strip.background = ggplot2::element_rect(fill = "white",
                                                     color = "black"),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10),
            plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5))
    return(plot_obj)
}