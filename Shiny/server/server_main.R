# Source individual server components
source("server/server_overview.R", local = TRUE)
source("server/server_pca.R", local = TRUE)
source("server/server_discriminant.R", local = TRUE)
source("server/server_anomaly.R", local = TRUE)
source("server/server_graph.R", local = TRUE)
source("server/server_wasserstein.R", local = TRUE)
source("server/server_marker.R", local = TRUE)
source("server/server_dimred.R", local = TRUE)
source("server/server_upload.R", local = TRUE)

# Define main server
server <- function(input, output, session) {
    
    # Reactive values for uploaded data
    uploaded_data <- reactiveValues()
    
    # Initialize default dataset selections to prevent errors
    observe({
        if(is.null(input$reference_dataset)) {
            updateSelectInput(session, "reference_dataset", selected = "reference_cells")
        }
        if(is.null(input$query_dataset)) {
            updateSelectInput(session, "query_dataset", selected = "query_cells")
        }
    })
    
    # Reactive values for datasets - these will be used throughout the app
    ref_data <- reactive({
        req(input$reference_dataset)
        if(startsWith(input$reference_dataset, "uploaded_ref_")) {
            uploaded_data[[input$reference_dataset]]
        } else {
            switch(input$reference_dataset,
                   "reference_cells" = reference_cells,
                   "reference_cells_subset" = reference_cells_subset)
        }
    })
    
    query_data <- reactive({
        req(input$query_dataset)
        if(startsWith(input$query_dataset, "uploaded_query_")) {
            uploaded_data[[input$query_dataset]]
        } else {
            query_cells
        }
    })
    
    # Helper function to get cell type column for dataset
    get_cell_type_col <- function(dataset_key) {
        if(is.null(dataset_key)) return(NULL)
        
        if(startsWith(dataset_key, "uploaded_ref_")) {
            uploaded_data$uploaded_ref_cell_type_col
        } else if(startsWith(dataset_key, "uploaded_query_")) {
            uploaded_data$uploaded_query_cell_type_col
        } else {
            NULL  # Will use default column selection logic
        }
    }
    
    # Update cell type column choices based on selected datasets (for PCA)
    observe({
        req(input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_cols <- colnames(colData(ref_data()))
            query_cols <- colnames(colData(query_data()))
            
            # Use uploaded dataset's configured cell type column if available
            ref_default <- get_cell_type_col(input$reference_dataset)
            query_default <- get_cell_type_col(input$query_dataset)
            
            updateSelectInput(session, "ref_cell_type_col",
                              choices = ref_cols,
                              selected = if(!is.null(ref_default) && ref_default %in% ref_cols) ref_default 
                              else if(length(ref_cols) > 0) ref_cols[1] else NULL)
            
            updateSelectInput(session, "query_cell_type_col", 
                              choices = query_cols,
                              selected = if(!is.null(query_default) && query_default %in% query_cols) query_default
                              else if(length(query_cols) > 0) query_cols[1] else NULL)
        }
    })
    
    # Update cell type column choices for graph analysis
    observe({
        req(input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_cols <- colnames(colData(ref_data()))
            query_cols <- colnames(colData(query_data()))
            
            ref_default <- get_cell_type_col(input$reference_dataset)
            query_default <- get_cell_type_col(input$query_dataset)
            
            updateSelectInput(session, "graph_ref_cell_type_col",
                              choices = ref_cols,
                              selected = if(!is.null(ref_default) && ref_default %in% ref_cols) ref_default
                              else if(length(ref_cols) > 0) ref_cols[1] else NULL)
            
            updateSelectInput(session, "graph_query_cell_type_col", 
                              choices = query_cols,
                              selected = if(!is.null(query_default) && query_default %in% query_cols) query_default
                              else if(length(query_cols) > 0) query_cols[1] else NULL)
        }
    })
    
    # Update cell type column choices for anomaly detection
    observe({
        req(input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_cols <- colnames(colData(ref_data()))
            query_cols <- colnames(colData(query_data()))
            
            ref_default <- get_cell_type_col(input$reference_dataset)
            query_default <- get_cell_type_col(input$query_dataset)
            
            updateSelectInput(session, "anomaly_ref_cell_type_col",
                              choices = ref_cols,
                              selected = if(!is.null(ref_default) && ref_default %in% ref_cols) ref_default
                              else if(length(ref_cols) > 0) ref_cols[1] else NULL)
            
            updateSelectInput(session, "anomaly_query_cell_type_col", 
                              choices = query_cols,
                              selected = if(!is.null(query_default) && query_default %in% query_cols) query_default
                              else if(length(query_cols) > 0) query_cols[1] else NULL)
        }
    })
    
    # Update cell type column choices for discriminant analysis
    observe({
        req(input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_cols <- colnames(colData(ref_data()))
            query_cols <- colnames(colData(query_data()))
            
            ref_default <- get_cell_type_col(input$reference_dataset)
            query_default <- get_cell_type_col(input$query_dataset)
            
            updateSelectInput(session, "disc_ref_cell_type_col",
                              choices = ref_cols,
                              selected = if(!is.null(ref_default) && ref_default %in% ref_cols) ref_default
                              else if(length(ref_cols) > 0) ref_cols[1] else NULL)
            
            updateSelectInput(session, "disc_query_cell_type_col", 
                              choices = query_cols,
                              selected = if(!is.null(query_default) && query_default %in% query_cols) query_default
                              else if(length(query_cols) > 0) query_cols[1] else NULL)
        }
    })
    
    # Update cell type column choices for Wasserstein analysis
    observe({
        req(input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_cols <- colnames(colData(ref_data()))
            query_cols <- colnames(colData(query_data()))
            
            ref_default <- get_cell_type_col(input$reference_dataset)
            query_default <- get_cell_type_col(input$query_dataset)
            
            updateSelectInput(session, "wass_ref_cell_type_col",
                              choices = ref_cols,
                              selected = if(!is.null(ref_default) && ref_default %in% ref_cols) ref_default
                              else if(length(ref_cols) > 0) ref_cols[1] else NULL)
            
            updateSelectInput(session, "wass_query_cell_type_col", 
                              choices = query_cols,
                              selected = if(!is.null(query_default) && query_default %in% query_cols) query_default
                              else if(length(query_cols) > 0) query_cols[1] else NULL)
        }
    })
    
    # Update cell type column choices for Marker analysis
    observe({
        req(input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_cols <- colnames(colData(ref_data()))
            query_cols <- colnames(colData(query_data()))
            
            ref_default <- get_cell_type_col(input$reference_dataset)
            query_default <- get_cell_type_col(input$query_dataset)
            
            updateSelectInput(session, "marker_ref_cell_type_col",
                              choices = ref_cols,
                              selected = if(!is.null(ref_default) && ref_default %in% ref_cols) ref_default
                              else if(length(ref_cols) > 0) ref_cols[1] else NULL)
            
            updateSelectInput(session, "marker_query_cell_type_col", 
                              choices = query_cols,
                              selected = if(!is.null(query_default) && query_default %in% query_cols) query_default
                              else if(length(query_cols) > 0) query_cols[1] else NULL)
        }
    })
    
    # Update available cell types (for PCA)
    observe({
        req(input$ref_cell_type_col, input$query_cell_type_col, input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_types <- unique(colData(ref_data())[[input$ref_cell_type_col]])
            query_types <- unique(colData(query_data())[[input$query_cell_type_col]])
            all_types <- sort(unique(c(ref_types, query_types)))
            all_types <- all_types[!is.na(all_types)]
            
            updateSelectInput(session, "cell_types", choices = all_types)
        }
    })
    
    # Update available cell types for discriminant analysis
    observe({
        req(input$disc_ref_cell_type_col, input$disc_query_cell_type_col, input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_types <- unique(colData(ref_data())[[input$disc_ref_cell_type_col]])
            query_types <- unique(colData(query_data())[[input$disc_query_cell_type_col]])
            all_types <- sort(unique(c(ref_types, query_types)))
            all_types <- all_types[!is.na(all_types)]
            
            updateSelectInput(session, "disc_cell_types", choices = all_types)
        }
    })
    
    # Update available cell types for anomaly detection
    observe({
        req(input$anomaly_ref_cell_type_col, input$anomaly_query_cell_type_col, input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_types <- unique(colData(ref_data())[[input$anomaly_ref_cell_type_col]])
            query_types <- unique(colData(query_data())[[input$anomaly_query_cell_type_col]])
            all_types <- sort(unique(c(ref_types, query_types)))
            all_types <- all_types[!is.na(all_types)]
            
            updateSelectInput(session, "anomaly_cell_types", choices = all_types)
            updateSelectInput(session, "anomaly_cell_type_plot", choices = c("Combined", all_types), selected = "Combined")
        }
    })
    
    # Update available cell types for graph analysis
    observe({
        req(input$graph_ref_cell_type_col, input$graph_query_cell_type_col, input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_types <- unique(colData(ref_data())[[input$graph_ref_cell_type_col]])
            query_types <- unique(colData(query_data())[[input$graph_query_cell_type_col]])
            all_types <- sort(unique(c(ref_types, query_types)))
            all_types <- all_types[!is.na(all_types)]
            
            updateSelectInput(session, "graph_cell_types", choices = all_types)
        }
    })
    
    # Update available cell types for Wasserstein analysis
    observe({
        req(input$wass_ref_cell_type_col, input$wass_query_cell_type_col, input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_types <- unique(colData(ref_data())[[input$wass_ref_cell_type_col]])
            query_types <- unique(colData(query_data())[[input$wass_query_cell_type_col]])
            all_types <- sort(unique(c(ref_types, query_types)))
            all_types <- all_types[!is.na(all_types)]
            
            updateSelectInput(session, "wass_cell_types", choices = all_types)
        }
    })
    
    # Update available cell types for Marker analysis
    observe({
        req(input$marker_ref_cell_type_col, input$marker_query_cell_type_col, input$reference_dataset, input$query_dataset)
        
        if(!is.null(ref_data()) && !is.null(query_data())) {
            ref_types <- unique(colData(ref_data())[[input$marker_ref_cell_type_col]])
            query_types <- unique(colData(query_data())[[input$marker_query_cell_type_col]])
            all_types <- sort(unique(c(ref_types, query_types)))
            all_types <- all_types[!is.na(all_types)]
            
            updateSelectInput(session, "marker_cell_type", choices = all_types)
        }
    })
    
    # Call modular server functions
    setup_overview_outputs(input, output, session, ref_data, query_data, uploaded_data)
    setup_upload_outputs(input, output, session, uploaded_data)
    setup_pca_outputs(input, output, session, ref_data, query_data)
    setup_discriminant_outputs(input, output, session, ref_data, query_data)
    setup_anomaly_outputs(input, output, session, ref_data, query_data)
    setup_graph_outputs(input, output, session, ref_data, query_data)
    setup_wasserstein_outputs(input, output, session, ref_data, query_data)
    setup_marker_outputs(input, output, session, ref_data, query_data)
    setup_dimred_outputs(input, output, session, ref_data, query_data)
}