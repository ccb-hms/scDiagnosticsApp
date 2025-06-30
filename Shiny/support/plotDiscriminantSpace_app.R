# App-adapted version of plot.calculateDiscriminantSpaceObject function
# This version includes all functionality from the original package function

plotDiscriminantSpace_app <- function(x,
                                      cell_types = NULL,
                                      dv_subset = NULL,
                                      lower_facet = "scatter",
                                      diagonal_facet = "ridge",
                                      upper_facet = "blank",
                                      ...) {
    
    # Load required libraries
    require(ggplot2)
    require(GGally)
    require(ggridges)
    require(MASS)
    require(rlang)
    
    # Check if query data is available in the object
    if (!("query_proj" %in% names(x))) {
        stop("There is no query data to plot.")
    }
    
    # Extract all cell types if not specified
    if (is.null(cell_types)) {
        cell_types <- unique(c(x[["ref_proj"]][["cell_type"]],
                               x[["query_proj"]][["cell_type"]]))
    }
    
    # Filter data to include only requested cell types
    ref_data <- x[["ref_proj"]][x[["ref_proj"]][["cell_type"]] %in% cell_types, ]
    query_data <- x[["query_proj"]][x[["query_proj"]][["cell_type"]] %in% cell_types, ]
    
    # Get discriminant vector columns
    dv_cols <- grep("^DV", colnames(ref_data), value = TRUE)
    total_dvs <- length(dv_cols)
    
    # Handle dv_subset parameter
    if (is.null(dv_subset)) {
        # Use all available discriminant vectors if dv_subset is NULL
        dv_subset <- seq_len(total_dvs)
    } else {
        # Check if specified dv_subset is valid
        if (max(dv_subset) > total_dvs) {
            stop(sprintf("Invalid dv_subset. Maximum available discriminant vectors: %d.", total_dvs))
        }
        if (min(dv_subset) < 1) {
            stop("Invalid dv_subset. Indices must be positive integers.")
        }
    }
    
    # Filter to the specified subset of discriminant vectors
    dv_cols <- dv_cols[dv_subset]
    
    # Check that we have at least one discriminant vector
    if (length(dv_cols) == 0) {
        stop("No valid discriminant vectors specified.")
    }
    
    # Combine reference and query data
    plot_data <- rbind(
        data.frame(ref_data, dataset = "Reference"),
        data.frame(query_data, dataset = "Query")
    )
    
    # Get the discriminant vector names for the plot
    plot_names <- paste0("DV", dv_subset)
    
    # Create a new data frame with selected DVs
    dv_df <- data.frame(matrix(0, nrow = nrow(plot_data), ncol = length(dv_subset)))
    colnames(dv_df) <- plot_names
    
    for (i in 1:length(dv_subset)) {
        dv_df[, i] <- plot_data[, paste0("DV", dv_subset[i])]
    }
    
    # Create a cell type dataset column for coloring
    cell_type_dataset <- paste(plot_data[["dataset"]], plot_data[["cell_type"]], sep = " ")
    
    # Define the order of cell type and dataset combinations
    order_combinations <- paste(
        rep(c("Reference", "Query"), length(cell_types)),
        rep(sort(cell_types), each = 2))
    
    cell_type_dataset <- factor(cell_type_dataset, levels = order_combinations)
    
    # Generate colors for cell types and datasets
    cell_type_colors <- generateColors(order_combinations, paired = TRUE)
    
    # Add the cell type dataset column to the DV data frame
    dv_df[["cell_type_dataset"]] <- cell_type_dataset
    
    # Create a simple plot to extract the legend using GGally::grab_legend
    legend_plot <- ggplot2::ggplot(dv_df, ggplot2::aes(x = dv_df[, 1], y = dv_df[, 2])) +
        ggplot2::geom_point(ggplot2::aes(color = cell_type_dataset)) +
        ggplot2::scale_color_manual(values = cell_type_colors, name = "Cell Type") +
        ggplot2::theme(
            legend.position = "right",
            legend.box = "vertical",
            legend.key = ggplot2::element_rect(fill = "white"),
            legend.background = ggplot2::element_blank()
        ) +
        ggplot2::guides(color = ggplot2::guide_legend(ncol = 1))
    
    # Scatterplot facet function
    .scatterFunc <- function(data, mapping, ...) {
        ggplot2::ggplot(data = data, mapping = mapping) +
            ggplot2::geom_point(alpha = 0.5, size = 1,
                                ggplot2::aes(color = cell_type_dataset)) +
            ggplot2::scale_color_manual(values = cell_type_colors) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                legend.position = "none"
            )
    }
    
    # Contour facet function
    .smoothContourFunc <- function(data, mapping, ...) {
        x_name <- rlang::as_name(mapping[["x"]])
        y_name <- rlang::as_name(mapping[["y"]])
        
        # Start with empty plot
        p <- ggplot2::ggplot() + ggplot2::theme_minimal()
        
        # Process each cell type separately
        for (ct in unique(data[["cell_type_dataset"]])) {
            subset_data <- data[data[["cell_type_dataset"]] == ct, ]
            
            # Skip if too few points
            if (nrow(subset_data) < 10) next
            
            # Use smoother contours with fewer levels
            p <- p + ggplot2::stat_density_2d(
                data = subset_data,
                mapping = ggplot2::aes(x = !!rlang::sym(x_name), y = !!rlang::sym(y_name)),
                contour = TRUE,
                adjust = 1.5,
                bins = 5,
                color = cell_type_colors[which(levels(data[["cell_type_dataset"]]) == ct)],
                linewidth = 0.5,
                na.rm = TRUE
            )
        }
        
        p + ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                legend.position = "none",
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                axis.title = ggplot2::element_blank()
            )
    }
    
    # Ellipse facet function
    .robustEllipseFunc <- function(data, mapping, ...) {
        # Function to calculate robust ellipses
        createEllipse <- function(d) {
            if (nrow(d) < 10) return(NULL)
            
            x_var <- rlang::as_name(mapping[["x"]])
            y_var <- rlang::as_name(mapping[["y"]])
            
            # Extract variables
            x <- d[[x_var]]
            y <- d[[y_var]]
            
            # Calculate robust center and covariance
            center <- c(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE))
            
            # Calculate ellipse points
            theta <- seq(0, 2 * pi, length.out = 100)
            ellipse <- MASS::cov.trob(cbind(x, y))
            
            # Get ellipse coordinates
            ev <- eigen(ellipse[["cov"]])
            a <- sqrt(ev[["values"]][1]) * 2.45
            b <- sqrt(ev[["values"]][2]) * 2.45
            
            # Create ellipse coordinates
            angle <- atan2(ev[["vectors"]][2, 1], ev[["vectors"]][1, 1])
            ellipse_x <- center[1] + a * cos(theta) * cos(angle) - b * sin(theta) * sin(angle)
            ellipse_y <- center[2] + a * cos(theta) * sin(angle) + b * sin(theta) * cos(angle)
            
            return(data.frame(x = ellipse_x, y = ellipse_y))
        }
        
        p <- ggplot2::ggplot(data, mapping)
        
        # Split by cell type and create robust ellipses
        for (ct in unique(data[["cell_type_dataset"]])) {
            subset_data <- data[data[["cell_type_dataset"]] == ct, ]
            ellipse_data <- createEllipse(subset_data)
            
            if (!is.null(ellipse_data)) {
                p <- p + ggplot2::geom_path(
                    data = ellipse_data,
                    ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
                    color = cell_type_colors[which(levels(data[["cell_type_dataset"]]) == ct)],
                    linewidth = 0.7
                )
            }
        }
        
        p + ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                legend.position = "none",
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                axis.title = ggplot2::element_blank()
            )
    }
    
    # Blank facet function
    .blankFunc <- function(data, mapping, ...) {
        ggplot2::ggplot() +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                legend.position = "none",
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                axis.title = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank()
            )
    }
    
    # Ridge diagonal facet
    .ridgeFunc <- function(data, mapping, ...) {
        x_var_name <- rlang::as_name(mapping[["x"]])
        
        # Create a long-format data frame for ridge plot
        plot_data <- data.frame(
            value = data[[x_var_name]],
            group = data[["cell_type_dataset"]]
        )
        
        # Create ridge plot with ggridges
        suppressMessages({
            p <- ggplot2::ggplot(plot_data,
                                 ggplot2::aes(x = .data[["value"]],
                                              y = .data[["group"]],
                                              fill = .data[["group"]])) +
                ggridges::geom_density_ridges(
                    alpha = 0.7,
                    scale = 2,
                    rel_min_height = 0.01,
                    quantile_lines = FALSE
                ) +
                ggplot2::scale_fill_manual(values = cell_type_colors) +
                ggplot2::scale_y_discrete(limits = rev(levels(cell_type_dataset))) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                    axis.title = ggplot2::element_blank(),
                    axis.text.y = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank(),
                    legend.position = "none"
                )
        })
        
        return(p)
    }
    
    # Density diagonal facet
    .densityFunc <- function(data, mapping, ...) {
        ggplot2::ggplot(data = data, mapping = mapping) +
            ggplot2::geom_density(ggplot2::aes(fill = cell_type_dataset,
                                               color = cell_type_dataset),
                                  alpha = 0.5) +
            ggplot2::scale_fill_manual(values = cell_type_colors) +
            ggplot2::scale_color_manual(values = cell_type_colors) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                legend.position = "none"
            )
    }
    
    # Boxplot diagonal facet
    .boxplotFunc <- function(data, mapping, ...) {
        x_name <- rlang::as_name(mapping[["x"]])
        
        # Create a long-format data frame for horizontal boxplot
        plot_data <- data.frame(
            value = data[[x_name]],
            group = data[["cell_type_dataset"]]
        )
        
        # Create horizontal boxplot
        ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["value"]],
                                                y = .data[["group"]],
                                                fill = .data[["group"]])) +
            ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.6) +
            ggplot2::scale_fill_manual(values = cell_type_colors) +
            ggplot2::scale_y_discrete(limits = rev(levels(cell_type_dataset))) +
            ggplot2::labs(x = "", y = "") +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                axis.title.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
                legend.position = "none",
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank()
            )
    }
    
    # Select the plot types based on user input
    if (lower_facet == "scatter") {
        lower_plot <- .scatterFunc
    } else if (lower_facet == "contour") {
        lower_plot <- .smoothContourFunc
    } else if (lower_facet == "ellipse") {
        lower_plot <- .robustEllipseFunc
    } else if (lower_facet == "blank") {
        lower_plot <- .blankFunc
    }
    
    if (diagonal_facet == "density") {
        diag_plot <- .densityFunc
    } else if (diagonal_facet == "boxplot") {
        diag_plot <- .boxplotFunc
    } else if (diagonal_facet == "ridge") {
        diag_plot <- .ridgeFunc
    } else if (diagonal_facet == "blank") {
        diag_plot <- .blankFunc
    }
    
    if (upper_facet == "blank") {
        upper_plot <- .blankFunc
    } else if (upper_facet == "scatter") {
        upper_plot <- .scatterFunc
    } else if (upper_facet == "contour") {
        upper_plot <- .smoothContourFunc
    } else if (upper_facet == "ellipse") {
        upper_plot <- .robustEllipseFunc
    }
    
    # Create pairs plot using GGally
    plot_obj <- suppressMessages(
        GGally::ggpairs(
            dv_df,
            columns = seq_len(length(dv_subset)),
            mapping = ggplot2::aes(color = cell_type_dataset),
            lower = list(continuous = lower_plot),
            upper = list(continuous = upper_plot),
            diag = list(continuous = diag_plot),
            progress = FALSE,
            legend = GGally::grab_legend(legend_plot)
        )
    )
    
    # Add black frame around facet titles
    plot_obj <- plot_obj +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.5),
            strip.text = ggplot2::element_text(color = "black"),
            plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::ggtitle("Discriminant Space Projection")
    
    return(plot_obj)
}