# Marker expression server functions
setup_marker_outputs <- function(input, output, session, ref_data, query_data) {
    
    # Get available genes (intersection of reference and query) - with proper error handling
    available_genes <- reactive({
        # Check if both datasets are available and inputs exist
        if(is.null(input$reference_dataset) || is.null(input$query_dataset)) {
            return(character(0))
        }
        
        tryCatch({
            ref_data_obj <- ref_data()
            query_data_obj <- query_data()
            
            if(is.null(ref_data_obj) || is.null(query_data_obj)) {
                return(character(0))
            }
            
            ref_genes <- rownames(ref_data_obj)
            query_genes <- rownames(query_data_obj)
            common_genes <- intersect(ref_genes, query_genes)
            sort(common_genes)  # Sort alphabetically
        }, error = function(e) {
            return(character(0))
        })
    })
    
    # Update gene choices when datasets change - with proper dependency
    observe({
        req(input$reference_dataset, input$query_dataset)
        
        genes <- available_genes()
        if(length(genes) > 0) {
            updateSelectizeInput(session, "gene_name",
                                 choices = setNames(genes, genes),
                                 server = TRUE)
        }
    })
    
    # Reactive for Marker Expression plot
    marker_plot_reactive <- eventReactive(input$generate_marker, {
        req(input$marker_ref_cell_type_col, input$marker_query_cell_type_col, 
            input$gene_name, input$marker_cell_type,
            input$reference_dataset, input$query_dataset)
        
        # Validate gene name
        if(input$gene_name == "" || is.null(input$gene_name)) {
            showNotification("Please select a gene name", type = "warning")
            return(NULL)
        }
        
        # Get selected cell type (single selection)
        selected_cell_type <- input$marker_cell_type
        if(is.null(selected_cell_type) || selected_cell_type == "") {
            showNotification("Please select a cell type", type = "warning")
            return(NULL)
        }
        
        # Call the marker expression function
        tryCatch({
            plotMarkerExpression_app(
                reference_data = ref_data(),
                query_data = query_data(),
                ref_cell_type_col = input$marker_ref_cell_type_col,
                query_cell_type_col = input$marker_query_cell_type_col,
                cell_type = selected_cell_type,  # Single cell type
                gene_name = input$gene_name,
                assay_name = input$marker_assay_name,
                normalization = input$normalization
            )
        }, error = function(e) {
            showNotification(paste("Error generating marker plot:", e$message), type = "error")
            NULL
        })
    })
    
    # Render Marker Expression plot
    output$marker_plot <- renderPlot({
        marker_plot_reactive()
    })
    
    # Calculate expression statistics
    marker_stats_reactive <- reactive({
        req(input$generate_marker > 0, input$gene_name, input$marker_cell_type,
            input$reference_dataset, input$query_dataset)
        
        tryCatch({
            # Get gene expression data
            ref_expr <- assay(ref_data(), input$marker_assay_name)[input$gene_name, ]
            query_expr <- assay(query_data(), input$marker_assay_name)[input$gene_name, ]
            
            # Get cell type specific expression
            ref_ct_indices <- which(colData(ref_data())[[input$marker_ref_cell_type_col]] == input$marker_cell_type)
            query_ct_indices <- which(colData(query_data())[[input$marker_query_cell_type_col]] == input$marker_cell_type)
            
            ref_ct_expr <- ref_expr[ref_ct_indices]
            query_ct_expr <- query_expr[query_ct_indices]
            
            list(
                ref_overall = ref_expr,
                query_overall = query_expr,
                ref_celltype = ref_ct_expr,
                query_celltype = query_ct_expr
            )
        }, error = function(e) {
            NULL
        })
    })
    
    # Marker statistics output
    output$marker_stats <- renderText({
        stats_data <- marker_stats_reactive()
        if(is.null(stats_data)) return("No statistics available")
        
        # Calculate summary statistics
        ref_overall_stats <- c(
            mean = mean(stats_data$ref_overall, na.rm = TRUE),
            median = median(stats_data$ref_overall, na.rm = TRUE),
            sd = sd(stats_data$ref_overall, na.rm = TRUE),
            n = length(stats_data$ref_overall)
        )
        
        query_overall_stats <- c(
            mean = mean(stats_data$query_overall, na.rm = TRUE),
            median = median(stats_data$query_overall, na.rm = TRUE),
            sd = sd(stats_data$query_overall, na.rm = TRUE),
            n = length(stats_data$query_overall)
        )
        
        ref_ct_stats <- c(
            mean = mean(stats_data$ref_celltype, na.rm = TRUE),
            median = median(stats_data$ref_celltype, na.rm = TRUE),
            sd = sd(stats_data$ref_celltype, na.rm = TRUE),
            n = length(stats_data$ref_celltype)
        )
        
        query_ct_stats <- c(
            mean = mean(stats_data$query_celltype, na.rm = TRUE),
            median = median(stats_data$query_celltype, na.rm = TRUE),
            sd = sd(stats_data$query_celltype, na.rm = TRUE),
            n = length(stats_data$query_celltype)
        )
        
        paste(
            paste("=== OVERALL DISTRIBUTION STATISTICS ==="),
            paste("Gene:", input$gene_name),
            paste("Assay:", input$marker_assay_name),
            "",
            paste("Reference Dataset (All Cells):"),
            paste("  Mean:", round(ref_overall_stats[["mean"]], 3)),
            paste("  Median:", round(ref_overall_stats[["median"]], 3)),
            paste("  SD:", round(ref_overall_stats[["sd"]], 3)),
            paste("  N:", ref_overall_stats[["n"]]),
            "",
            paste("Query Dataset (All Cells):"),
            paste("  Mean:", round(query_overall_stats[["mean"]], 3)),
            paste("  Median:", round(query_overall_stats[["median"]], 3)),
            paste("  SD:", round(query_overall_stats[["sd"]], 3)),
            paste("  N:", query_overall_stats[["n"]]),
            "",
            paste("=== CELL TYPE-SPECIFIC STATISTICS ==="),
            paste("Cell Type:", input$marker_cell_type),
            "",
            paste("Reference Dataset (Selected Cell Type):"),
            paste("  Mean:", round(ref_ct_stats[["mean"]], 3)),
            paste("  Median:", round(ref_ct_stats[["median"]], 3)),  
            paste("  SD:", round(ref_ct_stats[["sd"]], 3)),
            paste("  N:", ref_ct_stats[["n"]]),
            "",
            paste("Query Dataset (Selected Cell Type):"),
            paste("  Mean:", round(query_ct_stats[["mean"]], 3)),
            paste("  Median:", round(query_ct_stats[["median"]], 3)),
            paste("  SD:", round(query_ct_stats[["sd"]], 3)),
            paste("  N:", query_ct_stats[["n"]]),
            sep = "\n"
        )
    })
    
    # Marker analysis summary
    output$marker_analysis_summary <- renderText({
        req(input$generate_marker > 0, input$reference_dataset, input$query_dataset)
        
        dataset_name <- switch(input$reference_dataset,
                               "reference_cells" = "Reference Marrow Myeloid Cells",
                               "reference_cells_subset" = "Reference Marrow Myeloid Cells (Promonocytes Removed)",
                               input$reference_dataset)
        
        query_name <- switch(input$query_dataset,
                             "query_cells" = "Query Marrow Myeloid Cells",
                             input$query_dataset)
        
        normalization_desc <- switch(input$normalization,
                                     "z_score" = "Z-Score: Centers and scales expression within each dataset",
                                     "min_max" = "Min-Max: Scales expression to 0-1 range within each dataset",
                                     "rank" = "Quantile Rank: Converts to percentile ranks (0-100 scale)",
                                     "none" = "None: Preserves original expression values")
        
        paste(
            "=== ANALYSIS PARAMETERS ===",
            paste("Reference Dataset:", dataset_name),
            paste("Query Dataset:", query_name),
            paste("Reference Cell Type Column:", input$marker_ref_cell_type_col),
            paste("Query Cell Type Column:", input$marker_query_cell_type_col),
            paste("Target Gene:", input$gene_name),
            paste("Assay Used:", input$marker_assay_name),
            paste("Selected Cell Type:", input$marker_cell_type),
            "",
            "=== NORMALIZATION METHOD ===",
            normalization_desc,
            "",
            "=== PLOT INTERPRETATION ===",
            "Overall Distribution: Shows gene expression across all cells",
            "Cell Type-Specific: Shows expression within selected cell type only",
            "",
            "Similar distributions suggest good dataset alignment",
            "Different distributions may indicate:",
            "  ??? Technical batch effects",
            "  ??? Biological differences between datasets", 
            "  ??? Different cell populations or states",
            "  ??? Processing or normalization differences",
            "",
            "=== RECOMMENDATIONS ===",
            "??? Compare multiple marker genes for comprehensive assessment",
            "??? Try different normalization methods if distributions differ",
            "??? Consider batch correction if systematic differences observed",
            sep = "\n"
        )
    })
    
    # Download handler for marker plot
    output$download_marker_plot <- downloadHandler(
        filename = function() {
            paste0("Marker_expression_", gsub("[^A-Za-z0-9]", "_", input$gene_name), "_", 
                   gsub("[^A-Za-z0-9]", "_", input$marker_cell_type), "_", Sys.Date(), ".png")
        },
        content = function(file) {
            plot_obj <- marker_plot_reactive()
            if(!is.null(plot_obj)) {
                ggsave(file, plot_obj, width = 12, height = 8, dpi = 300)
            }
        }
    )
}