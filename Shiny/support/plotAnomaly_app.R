# App-adapted version of plot.detectAnomalyObject function
# This version includes all functionality from the original package function

plotAnomaly_app <- function(x,
                            cell_type = NULL,
                            pc_subset = NULL,
                            data_type = "query",
                            n_tree = 500,
                            upper_facet = "blank",
                            diagonal_facet = "density",
                            ...) {
    
    # Load required libraries
    require(ggplot2)
    require(GGally)
    require(ggridges)
    require(MASS)
    require(isotree)
    require(rlang)
    require(scales)
    
    # Check if PCA was used for computations
    if (!("var_explained" %in% names(x[[names(x)[1]]]))) {
        stop("The plot function can only be used if 'pc_subset' is not NULL.")
    }
    
    # Check input for cell type
    if (is.null(cell_type)) {
        cell_type <- "Combined"
    } else {
        if (!(cell_type %in% names(x))) {
            stop("'cell_type' is not available in 'x'.")
        }
    }
    
    # Check input for pc_subset
    if (!is.null(pc_subset)) {
        if (!all(pc_subset %in% seq_len(ncol(x[[cell_type]][["reference_mat_subset"]])))) {
            stop("'pc_subset' is out of range.")
        }
    } else {
        pc_subset <- seq_len(ncol(x[[cell_type]][["reference_mat_subset"]]))
    }
    
    # Filter and prepare data based on data type
    if (is.null(x[[cell_type]][["query_mat_subset"]]) && data_type == "query") {
        stop("There is no query data available in the 'detectAnomaly' object.")
    } else {
        if (data_type == "query") {
            data_subset <- x[[cell_type]][["query_mat_subset"]][, pc_subset, drop = FALSE]
            anomaly <- x[[cell_type]][["query_anomaly"]]
        } else if (data_type == "reference") {
            data_subset <- x[[cell_type]][["reference_mat_subset"]][, pc_subset, drop = FALSE]
            anomaly <- x[[cell_type]][["reference_anomaly"]]
        }
    }
    
    # Add variance explained to column names
    colnames(data_subset) <- paste0(
        "PC", pc_subset,
        " (", sprintf("%.1f%%", x[[cell_type]][["var_explained"]][pc_subset]), ")")
    
    # Create a data frame with PC values and anomaly info
    pc_df <- data.frame(data_subset)
    colnames(pc_df) <- colnames(data_subset)
    pc_df[["anomaly"]] <- factor(anomaly)
    
    # Colors for anomalous and non-anomalous data
    anomaly_colors <- c("TRUE" = "red", "FALSE" = "black")
    anomaly_fill_colors <- c("TRUE" = "#FFB6B6", "FALSE" = "#DDDDDD")
    
    # Train isolation forest for each PC combination
    isolation_forests <- list()
    for (i in seq_along(pc_subset)) {
        for (j in seq_along(pc_subset)) {
            if (i < j) {
                pc_i <- paste0("PC", pc_subset[i])
                pc_j <- paste0("PC", pc_subset[j])
                
                if (data_type == "reference") {
                    train_data <- x[[cell_type]][["reference_mat_subset"]][, c(pc_i, pc_j)]
                } else {
                    train_data <- x[[cell_type]][["reference_mat_subset"]][, c(pc_i, pc_j)]
                }
                
                isolation_forests[[paste(pc_i, pc_j, sep = "-")]] <-
                    isotree::isolation.forest(train_data, ntree = n_tree, ...)
            }
        }
    }
    
    # Function to create scatter plot with isolation forest background
    .anomalyScatterFunc <- function(data, mapping, ...) {
        x_name <- rlang::as_name(mapping[["x"]])
        y_name <- rlang::as_name(mapping[["y"]])
        
        # Get PC numbers from column names
        pc_x <- as.numeric(sub("PC([0-9]+).*", "\\1", x_name))
        pc_y <- as.numeric(sub("PC([0-9]+).*", "\\1", y_name))
        
        # Get isolation forest for this PC combination
        if_key <- paste0("PC", min(pc_x, pc_y), "-PC", max(pc_x, pc_y))
        if_model <- isolation_forests[[if_key]]
        
        # Get ranges for grid
        x_range <- range(data[[x_name]])
        y_range <- range(data[[y_name]])
        x_buffer <- 0.1 * diff(x_range)
        y_buffer <- 0.1 * diff(y_range)
        
        # Create grid for background
        x_seq <- seq(x_range[1] - x_buffer, x_range[2] + x_buffer, length.out = 150)
        y_seq <- seq(y_range[1] - y_buffer, y_range[2] + y_buffer, length.out = 150)
        grid <- expand.grid(x = x_seq, y = y_seq)
        
        # Calculate step sizes for better tile width/height
        x_step <- diff(x_seq)[1]
        y_step <- diff(y_seq)[1]
        
        # Get anomaly scores for grid
        grid_data <- data.frame(grid)
        colnames(grid_data) <- c(paste0("PC", pc_x), paste0("PC", pc_y))
        
        scores <- predict(if_model, newdata = grid_data, type = "score")
        
        # Create background data for plotting
        background_data <- data.frame(
            x_value = grid[["x"]],
            y_value = grid[["y"]],
            anomaly_score = scores,
            width = x_step,
            height = y_step
        )
        
        # Create gradient colors
        gradient_colors <- c("#2CA25F", "#FFFFE5", "#FF8585")
        stops <- c(0, 0.55, 1)
        
        # Create plot with background and points
        p <- ggplot2::ggplot() +
            ggplot2::geom_tile(
                data = background_data,
                ggplot2::aes(
                    x = .data[["x_value"]],
                    y = .data[["y_value"]],
                    fill = .data[["anomaly_score"]],
                    width = .data[["width"]],
                    height = .data[["height"]]
                )
            ) +
            ggplot2::scale_fill_gradientn(
                colors = gradient_colors,
                values = scales::rescale(stops),
                limits = c(0, 1),
                guide = "none"
            ) +
            ggplot2::geom_point(
                data = data,
                ggplot2::aes(
                    x = !!rlang::sym(x_name),
                    y = !!rlang::sym(y_name),
                    color = anomaly
                ),
                shape = 16,
                size = 1.5,
                alpha = 0.45
            ) +
            ggplot2::scale_color_manual(
                values = anomaly_colors,
                name = "Anomalous"
            ) +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::scale_y_continuous(expand = c(0, 0)) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(),
                legend.position = "none"
            )
        
        return(p)
    }
    
    # Function for diagonal density plots with separate anomaly curves
    .densityDiagFunc <- function(data, mapping, ...) {
        x_name <- rlang::as_name(mapping[["x"]])
        
        ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x_name))) +
            ggplot2::geom_density(
                ggplot2::aes(fill = anomaly, color = anomaly),
                alpha = 0.3,
                linewidth = 0.8
            ) +
            ggplot2::scale_color_manual(values = anomaly_colors) +
            ggplot2::scale_fill_manual(values = anomaly_fill_colors) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                legend.position = "none"
            )
    }
    
    # Function for diagonal ridge plots with separate anomaly distributions
    .ridgeDiagFunc <- function(data, mapping, ...) {
        x_name <- rlang::as_name(mapping[["x"]])
        
        plot_data <- data.frame(
            value = data[[x_name]],
            anomaly = data[["anomaly"]]
        )
        
        suppressMessages({
            p <- ggplot2::ggplot(plot_data,
                                 ggplot2::aes(x = .data[["value"]],
                                              y = .data[["anomaly"]],
                                              fill = .data[["anomaly"]],
                                              color = .data[["anomaly"]])) +
                ggridges::geom_density_ridges(
                    alpha = 0.6,
                    scale = 1.8,
                    rel_min_height = 0.01,
                    quantile_lines = FALSE
                ) +
                ggplot2::scale_fill_manual(values = anomaly_fill_colors) +
                ggplot2::scale_color_manual(values = anomaly_colors) +
                ggplot2::scale_y_discrete(limits = rev(levels(data[["anomaly"]]))) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                    axis.title = ggplot2::element_blank(),
                    axis.text.y = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank(),
                    legend.position = "none",
                    panel.grid = ggplot2::element_blank()
                )
        })
        
        return(p)
    }
    
    # Function for diagonal boxplot with separate anomaly distributions
    .boxplotDiagFunc <- function(data, mapping, ...) {
        x_name <- rlang::as_name(mapping[["x"]])
        
        plot_data <- data.frame(
            value = data[[x_name]],
            anomaly = data[["anomaly"]]
        )
        
        ggplot2::ggplot(plot_data,
                        ggplot2::aes(x = .data[["value"]],
                                     y = .data[["anomaly"]],
                                     fill = .data[["anomaly"]])) +
            ggplot2::geom_boxplot(
                alpha = 0.7,
                outlier.size = 0.5,
                width = 0.6,
                color = "black"
            ) +
            ggplot2::scale_fill_manual(values = anomaly_fill_colors) +
            ggplot2::scale_y_discrete(limits = rev(levels(data[["anomaly"]]))) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(),
                axis.title = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                legend.position = "none"
            )
    }
    
    # Function for blank panels
    .blankFunc <- function(data, mapping, ...) {
        ggplot2::ggplot() +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(),
                axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank()
            )
    }
    
    # Function for contour plots separated by anomaly status
    .contourFunc <- function(data, mapping, ...) {
        x_name <- rlang::as_name(mapping[["x"]])
        y_name <- rlang::as_name(mapping[["y"]])
        
        p <- ggplot2::ggplot(data, ggplot2::aes(
            x = !!rlang::sym(x_name),
            y = !!rlang::sym(y_name))) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(),
                axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                legend.position = "none"
            )
        
        # Add contours for non-anomalous data
        data_normal <- data[data[["anomaly"]] == "FALSE",]
        if (nrow(data_normal) >= 10) {
            p <- p + ggplot2::stat_density_2d(
                data = data_normal,
                ggplot2::aes(
                    x = !!rlang::sym(x_name),
                    y = !!rlang::sym(y_name)
                ),
                color = anomaly_colors[["FALSE"]],
                linewidth = 0.5,
                contour = TRUE,
                bins = 5
            )
        }
        
        # Add contours for anomalous data
        data_anomaly <- data[data[["anomaly"]] == "TRUE",]
        if (nrow(data_anomaly) >= 10) {
            p <- p + ggplot2::stat_density_2d(
                data = data_anomaly,
                ggplot2::aes(
                    x = !!rlang::sym(x_name),
                    y = !!rlang::sym(y_name)
                ),
                color = anomaly_colors[["TRUE"]],
                linewidth = 0.5,
                contour = TRUE,
                bins = 5
            )
        }
        
        return(p)
    }
    
    # Function for ellipse plots separated by anomaly status
    .ellipseFunc <- function(data, mapping, ...) {
        x_name <- rlang::as_name(mapping[["x"]])
        y_name <- rlang::as_name(mapping[["y"]])
        
        # Function to calculate and create ellipses
        createEllipse <- function(d) {
            if (nrow(d) < 10) return(NULL)
            
            x <- d[[x_name]]
            y <- d[[y_name]]
            
            cov_mat <- cov(cbind(x, y), use = "pairwise.complete.obs")
            center <- c(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE))
            
            ellipse <- MASS::cov.trob(cbind(x, y))
            
            ev <- eigen(ellipse[["cov"]])
            a <- sqrt(ev[["values"]][1]) * 2.45
            b <- sqrt(ev[["values"]][2]) * 2.45
            
            angle <- atan2(ev[["vectors"]][2, 1], ev[["vectors"]][1, 1])
            
            theta <- seq(0, 2 * pi, length.out = 100)
            ellipse_x <- center[1] + a * cos(theta) * cos(angle) - b * sin(theta) * sin(angle)
            ellipse_y <- center[2] + a * cos(theta) * sin(angle) + b * sin(theta) * cos(angle)
            
            return(data.frame(x = ellipse_x, y = ellipse_y))
        }
        
        p <- ggplot2::ggplot() +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(),
                axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                legend.position = "none"
            )
        
        # Create ellipses for normal data
        data_normal <- data[data[["anomaly"]] == "FALSE",]
        ellipse_normal <- createEllipse(data_normal)
        if (!is.null(ellipse_normal)) {
            p <- p + ggplot2::geom_path(
                data = ellipse_normal,
                ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
                color = anomaly_colors[["FALSE"]],
                linewidth = 0.7
            )
        }
        
        # Create ellipses for anomalous data
        data_anomaly <- data[data[["anomaly"]] == "TRUE",]
        ellipse_anomaly <- createEllipse(data_anomaly)
        if (!is.null(ellipse_anomaly)) {
            p <- p + ggplot2::geom_path(
                data = ellipse_anomaly,
                ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
                color = anomaly_colors[["TRUE"]],
                linewidth = 0.7
            )
        }
        
        return(p)
    }
    
    # Select diagonal facet type based on input
    if (diagonal_facet == "density") {
        diag_func <- .densityDiagFunc
    } else if (diagonal_facet == "ridge") {
        diag_func <- .ridgeDiagFunc
    } else if (diagonal_facet == "boxplot") {
        diag_func <- .boxplotDiagFunc
    } else if (diagonal_facet == "blank") {
        diag_func <- .blankFunc
    }
    
    # Select upper facet type based on input
    if (upper_facet == "blank") {
        upper_func <- .blankFunc
    } else if (upper_facet == "contour") {
        upper_func <- .contourFunc
    } else if (upper_facet == "ellipse") {
        upper_func <- .ellipseFunc
    }
    
    # Create a simple plot to extract the legend
    legend_plot <- ggplot2::ggplot(pc_df, ggplot2::aes(x = pc_df[, 1], y = pc_df[, 2])) +
        ggplot2::geom_point(ggplot2::aes(color = anomaly)) +
        ggplot2::scale_color_manual(values = anomaly_colors, name = "Anomalous") +
        ggplot2::theme(
            legend.position = "right",
            legend.box = "vertical",
            legend.key = ggplot2::element_rect(fill = "white")
        )
    
    # Create pairs plot using GGally
    plot_obj <- suppressMessages(
        GGally::ggpairs(
            pc_df,
            columns = seq_len(length(pc_subset)),
            mapping = ggplot2::aes(color = anomaly),
            lower = list(continuous = .anomalyScatterFunc),
            diag = list(continuous = diag_func),
            upper = list(continuous = upper_func),
            progress = FALSE,
            legend = GGally::grab_legend(legend_plot)
        )
    )
    
    # Add black frame around facet titles and center the title
    plot_obj <- plot_obj +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.5),
            strip.text = ggplot2::element_text(color = "black"),
            plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::ggtitle(paste0("Isolation Forest Anomaly Plot: ", cell_type))
    
    return(plot_obj)
}