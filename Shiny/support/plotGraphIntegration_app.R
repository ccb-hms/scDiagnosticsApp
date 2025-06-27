# App-adapted version of plot.calculateGraphIntegrationObject function
# This version includes all functionality from the original package function

plotGraphIntegration_app <- function(x,
                                     plot_type = "community_network",
                                     color_by = "cell_type",
                                     max_nodes = 2000,
                                     point_size = 0.8,
                                     exclude_reference_only = FALSE) {
    
    # Load required libraries
    require(ggplot2)
    require(igraph)
    
    # Define colors
    colors <- list(
        high_query_prop = "#DC143C",
        cross_mixing = "#FF8C00",
        local_inconsistent = "#9932CC",
        well_integrated = "#2E8B57",
        reference_only = "#4169E1",
        small_community = "gray70",
        excellent = "#228B22",
        concerning = "#FF6347",
        problematic = "#DC143C"
    )
    
    # Helper function for pluralization
    pluralize <- function(count, singular, plural = paste0(singular, "s")) {
        if (count == 1) return(paste(count, singular))
        else return(paste(count, plural))
    }
    
    if (plot_type == "community_network") {
        community_composition <- x[["community_composition"]]
        
        if (nrow(community_composition) == 0) {
            p <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5,
                                  label = "No communities to display",
                                  size = 8, color = "gray50") +
                ggplot2::labs(title = "Community Network") +
                ggplot2::theme_void() +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
            return(p)
        }
        
        # Filter out reference-only communities if requested
        if (exclude_reference_only) {
            reference_only_mask <- community_composition[["query_proportion"]] < 0.1
            community_composition <- community_composition[!reference_only_mask, ]
        }
        
        if (nrow(community_composition) == 0) {
            p <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5,
                                  label = "No communities to display after filtering",
                                  size = 8, color = "gray50") +
                ggplot2::labs(title = "Community Network") +
                ggplot2::theme_void() +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
            return(p)
        }
        
        # Calculate community centroids
        community_centroids <- data.frame(
            community = integer(0),
            x = numeric(0),
            y = numeric(0),
            stringsAsFactors = FALSE
        )
        
        for (comm in community_composition[["community"]]) {
            comm_cells <- x[["cell_info"]][["community"]] == comm
            if (sum(comm_cells) > 0) {
                centroid_x <- mean(x[["graph_info"]][["layout"]][comm_cells, 1])
                centroid_y <- mean(x[["graph_info"]][["layout"]][comm_cells, 2])
                community_centroids <- rbind(community_centroids, data.frame(
                    community = comm,
                    x = centroid_x,
                    y = centroid_y,
                    stringsAsFactors = FALSE
                ))
            }
        }
        
        # ========== PROPER EDGE CALCULATION ==========
        # Calculate adaptive threshold based on overall graph connectivity
        original_edges <- x[["graph_info"]][["edges"]]
        cell_communities <- x[["cell_info"]][["community"]]
        total_cells <- length(cell_communities)
        total_edges <- nrow(original_edges)
        
        # Base connectivity rate in the graph
        base_connectivity <- total_edges / choose(total_cells, 2)
        
        # Adaptive threshold: 2x the base rate, with bounds
        adaptive_threshold <- max(0.005, min(0.025, 1.5 * base_connectivity))
        
        # Create edge lookup for fast checking
        edge_lookup <- paste(original_edges[, 1], original_edges[, 2], sep = "_")
        edge_lookup <- c(edge_lookup, paste(original_edges[, 2], original_edges[, 1], sep = "_"))
        edge_set <- as.environment(as.list(setNames(rep(TRUE, length(edge_lookup)), edge_lookup)))
        
        edge_data <- data.frame(x = numeric(0), y = numeric(0),
                                xend = numeric(0), yend = numeric(0),
                                connection_strength = numeric(0))
        
        communities <- community_composition[["community"]]
        n_communities <- length(communities)
        
        # Check all pairs of communities for connections
        if (n_communities > 1) {
            for (i in seq_len(n_communities - 1)) {
                for (j in seq(i + 1, n_communities)) {
                    comm1 <- communities[i]
                    comm2 <- communities[j]
                    
                    cells_comm1 <- which(cell_communities == comm1)
                    cells_comm2 <- which(cell_communities == comm2)
                    
                    if (length(cells_comm1) > 0 && length(cells_comm2) > 0) {
                        # Adaptive sampling based on community sizes
                        max_samples <- min(250, length(cells_comm1) * length(cells_comm2))
                        n_samples <- max(100, max_samples)
                        
                        sampled_comm1 <- sample(cells_comm1, n_samples, replace = TRUE)
                        sampled_comm2 <- sample(cells_comm2, n_samples, replace = TRUE)
                        
                        # Count connections
                        connections <- 0
                        for (k in seq_len(n_samples)) {
                            pair_key <- paste(sampled_comm1[k], sampled_comm2[k], sep = "_")
                            if (exists(pair_key, envir = edge_set)) {
                                connections <- connections + 1
                            }
                        }
                        
                        connection_rate <- connections / n_samples
                        
                        # Use adaptive threshold
                        if (connection_rate >= adaptive_threshold) {
                            comm1_idx <- which(community_centroids$community == comm1)
                            comm2_idx <- which(community_centroids$community == comm2)
                            
                            if (length(comm1_idx) > 0 && length(comm2_idx) > 0) {
                                edge_data <- rbind(edge_data, data.frame(
                                    x = community_centroids$x[comm1_idx],
                                    y = community_centroids$y[comm1_idx],
                                    xend = community_centroids$x[comm2_idx],
                                    yend = community_centroids$y[comm2_idx],
                                    connection_strength = connection_rate,
                                    stringsAsFactors = FALSE
                                ))
                            }
                        }
                    }
                }
            }
        }
        # ========== END EDGE CALCULATION ==========
        
        # Prepare node data
        node_data <- merge(community_centroids, community_composition, by = "community")
        
        if (color_by == "cell_type") {
            # Get all cell types and generate paired colors
            all_cell_types <- sort(unique(x[["cell_info"]][["cell_type"]]))
            cell_types_cases <- c()
            for(cell_type_case in all_cell_types){
                cell_types_cases <- c(cell_types_cases,
                                      paste(cell_type_case, "(Mixed)"),
                                      paste(cell_type_case, "(Pure)"))
            }
            paired_colors <- generateColors(cell_types_cases, paired = TRUE)
            
            # Determine dominant cell type and whether community is mixed
            node_data[["dominant_cell_type"]] <- sapply(node_data[["community"]], function(comm) {
                comm_cells <- x[["cell_info"]][x[["cell_info"]][["community"]] == comm, ]
                cell_type_counts <- table(comm_cells[["cell_type"]])
                names(cell_type_counts)[which.max(cell_type_counts)]
            })
            node_data[["is_mixed"]] <- node_data[["n_cell_types"]] > 1
            node_data[["cell_type_case"]] <- paste(node_data[["dominant_cell_type"]],
                                                   ifelse(node_data[["is_mixed"]], "(Mixed)", "(Pure)"))
            
            # Classify communities for shapes
            node_data[["shape_type"]] <- "Well Integrated"
            node_data[["shape_type"]][node_data[["community"]] %in% x[["cross_type_mixing"]][["community"]]] <- "Cross-Type Mixing"
            node_data[["shape_type"]][node_data[["community"]] %in% x[["high_query_prop_analysis"]][["community"]]] <- "High Query Proportion"
            node_data[["shape_type"]][node_data[["query_proportion"]] < 0.1 & node_data[["shape_type"]] == "Well Integrated"] <- "High Reference Proportion"
            
            p <- ggplot2::ggplot() +
                {if (nrow(edge_data) > 0) {
                    ggplot2::geom_segment(data = edge_data,
                                          ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                                       xend = .data[["xend"]], yend = .data[["yend"]],
                                                       alpha = .data[["connection_strength"]]),
                                          color = "gray60", linewidth = 0.8)
                }} +
                ggplot2::geom_point(data = node_data,
                                    ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                                 color = .data[["cell_type_case"]],
                                                 shape = .data[["shape_type"]],
                                                 size = .data[["total_cells"]]),
                                    alpha = 0.8) +
                ggplot2::scale_color_manual(values = paired_colors, name = "Dominant Cell Type") +
                ggplot2::scale_shape_manual(
                    values = c("Well Integrated" = 16, "Cross-Type Mixing" = 18,
                               "High Query Proportion" = 17, "High Reference Proportion" = 15),
                    name = "Community Classification"
                ) +
                ggplot2::scale_size_continuous(range = c(3, 10), trans = "sqrt", guide = "none") +
                ggplot2::scale_alpha_continuous(range = c(0.3, 1.0), guide = "none") +
                ggplot2::geom_text(data = node_data,
                                   ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                                label = .data[["community"]]),
                                   color = "black", size = 3, fontface = "bold") +
                ggplot2::labs(
                    title = "Community Network: Cell Type Composition",
                    subtitle = paste0("Showing ", nrow(node_data), " communities with ",
                                      nrow(edge_data), " connections"),
                    caption = "Node size proportional to # cells; Shape indicates community classification; Edge opacity proportional to connection strength"
                ) +
                ggplot2::theme_void() +
                ggplot2::theme(
                    plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                    plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
                    legend.position = "right"
                )
            
        } else if (color_by == "community_type") {
            # Original community type classification
            node_data[["issue_type"]] <- "Well Integrated"
            node_data[["issue_type"]][node_data[["community"]] %in% x[["cross_type_mixing"]][["community"]]] <- "Cross-Type Mixing"
            node_data[["issue_type"]][node_data[["community"]] %in% x[["high_query_prop_analysis"]][["community"]]] <- "High Query Proportion"
            node_data[["issue_type"]][node_data[["query_proportion"]] < 0.1 & node_data[["issue_type"]] == "Well Integrated"] <- "High Reference Proportion"
            
            p <- ggplot2::ggplot() +
                {if (nrow(edge_data) > 0) {
                    ggplot2::geom_segment(data = edge_data,
                                          ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                                       xend = .data[["xend"]], yend = .data[["yend"]],
                                                       alpha = .data[["connection_strength"]]),
                                          color = "gray60", linewidth = 0.8)
                }} +
                ggplot2::geom_point(data = node_data,
                                    ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                                 color = .data[["issue_type"]],
                                                 size = .data[["total_cells"]]),
                                    alpha = 0.8) +
                ggplot2::scale_color_manual(
                    values = c("High Query Proportion" = colors[["high_query_prop"]],
                               "Cross-Type Mixing" = colors[["cross_mixing"]],
                               "Well Integrated" = colors[["well_integrated"]],
                               "High Reference Proportion" = colors[["reference_only"]]),
                    name = "Community Type"
                ) +
                ggplot2::scale_size_continuous(range = c(3, 10), trans = "sqrt", guide = "none") +
                ggplot2::scale_alpha_continuous(range = c(0.3, 1.0), guide = "none") +
                ggplot2::geom_text(data = node_data,
                                   ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                                label = .data[["community"]]),
                                   color = "black", size = 3, fontface = "bold") +
                ggplot2::labs(
                    title = "Community Network: Inter-Community Connection Strength",
                    subtitle = paste0("Showing ", nrow(node_data), " communities with ",
                                      nrow(edge_data), " connections"),
                    caption = "Node size proportional to # cells in community; Edge opacity proportional to connection strength"
                ) +
                ggplot2::theme_void() +
                ggplot2::theme(
                    plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                    plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
                    legend.position = "right"
                )
        }
        
    } else if (plot_type == "summary") {
        summary_data <- x[["annotation_consistency"]]
        summary_data[["total_issues"]] <- summary_data[["high_query_prop_cells"]] + summary_data[["true_cross_mixing_cells"]]
        
        local_summary <- x[["local_inconsistency_summary"]]
        summary_data <- merge(summary_data, local_summary[, c("cell_type", "locally_inconsistent_cells")], by = "cell_type")
        summary_data[["total_issues"]] <- summary_data[["total_issues"]] + summary_data[["locally_inconsistent_cells"]]
        summary_data[["total_issue_rate"]] <- summary_data[["total_issues"]] / summary_data[["total_query_cells"]]
        
        summary_data[["issue_category"]] <- cut(summary_data[["total_issue_rate"]],
                                                breaks = c(-Inf, 0.05, 0.15, Inf),
                                                labels = c("Excellent", "Concerning", "Problematic"))
        
        p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = reorder(.data[["cell_type"]], .data[["total_issue_rate"]]))) +
            ggplot2::geom_col(ggplot2::aes(y = .data[["total_issue_rate"]], fill = .data[["issue_category"]]),
                              alpha = 0.8, color = "black", linewidth = 0.3) +
            ggplot2::scale_fill_manual(
                values = c("Excellent" = colors[["excellent"]],
                           "Concerning" = colors[["concerning"]],
                           "Problematic" = colors[["problematic"]]),
                name = "Issue Level"
            ) +
            ggplot2::geom_text(ggplot2::aes(y = .data[["total_issue_rate"]],
                                            label = paste0(.data[["total_issues"]], "/", .data[["total_query_cells"]])),
                               vjust = -0.5, size = 3) +
            ggplot2::labs(
                title = "Total Annotation Issues by Cell Type",
                subtitle = "Includes query-only, cross-mixing, and local inconsistencies",
                x = "Cell Type",
                y = "Proportion of Query Cells with Issues",
                caption = paste0("Mean issue rate: ", round(mean(summary_data[["total_issue_rate"]]), 3))
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                legend.position = "bottom"
            ) +
            ggplot2::coord_cartesian(ylim = c(0, max(1, max(summary_data[["total_issue_rate"]]) * 1.1)))
        
    } else if (plot_type == "community_data") {
        community_composition <- x[["community_composition"]]
        
        if (nrow(community_composition) == 0) {
            p <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No communities to display",
                                  size = 8, color = "gray50") +
                ggplot2::labs(title = "Community Overview") +
                ggplot2::theme_void() +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
            return(p)
        }
        
        # Filter out reference-only communities if requested
        if (exclude_reference_only) {
            reference_only_mask <- community_composition[["query_proportion"]] < 0.1
            community_composition <- community_composition[!reference_only_mask, ]
        }
        
        if (nrow(community_composition) == 0) {
            p <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No communities to display after filtering",
                                  size = 8, color = "gray50") +
                ggplot2::labs(title = "Community Overview") +
                ggplot2::theme_void() +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
            return(p)
        }
        
        # Classify communities
        community_composition[["issue_type"]] <- "Well Integrated"
        community_composition[["issue_type"]][community_composition[["community"]] %in% x[["cross_type_mixing"]][["community"]]] <- "Cross-Type Mixing"
        community_composition[["issue_type"]][community_composition[["community"]] %in% x[["high_query_prop_analysis"]][["community"]]] <- "High Query Proportion"
        community_composition[["issue_type"]][community_composition[["query_proportion"]] < 0.1 & community_composition[["issue_type"]] == "Well Integrated"] <- "High Reference Proportion"
        
        p <- ggplot2::ggplot(community_composition, ggplot2::aes(x = .data[["query_proportion"]], y = .data[["total_cells"]])) +
            ggplot2::geom_point(ggplot2::aes(color = .data[["issue_type"]], size = .data[["n_cell_types"]]), alpha = 0.7) +
            ggplot2::scale_color_manual(
                values = c("High Query Proportion" = colors[["high_query_prop"]],
                           "Cross-Type Mixing" = colors[["cross_mixing"]],
                           "Well Integrated" = colors[["well_integrated"]],
                           "High Reference Proportion" = colors[["reference_only"]]),
                name = "Community Type"
            ) +
            ggplot2::scale_size_continuous(range = c(2, 8), name = "# Cell Types") +
            ggplot2::geom_text(ggplot2::aes(label = .data[["community"]], vjust = -0.75), size = 3) +
            ggplot2::labs(
                title = "Community Data Overview",
                x = "Query Proportion",
                y = "Total Cells in Community"
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                legend.position = "right"
            )
    }
    
    return(p)
}