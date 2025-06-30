# Source individual server components
source("server/server_overview.R", local = TRUE)
source("server/server_pca.R", local = TRUE)
source("server/server_discriminant.R", local = TRUE)
source("server/server_graph.R", local = TRUE)
source("server/server_anomaly.R", local = TRUE)

# Define main server
server <- function(input, output, session) {
    
    # Reactive values for datasets
    ref_data <- reactive({
        switch(input$reference_dataset,
               "reference_cells" = reference_cells,
               "reference_cells_subset" = reference_cells_subset)
    })
    
    query_data <- reactive({
        query_cells
    })
    
    # Update cell type column choices based on selected datasets (for PCA)
    observe({
        ref_cols <- colnames(colData(ref_data()))
        query_cols <- colnames(colData(query_data()))
        
        # Filter for likely cell type columns
        ref_cell_type_cols <- ref_cols
        query_cell_type_cols <- query_cols
        
        updateSelectInput(session, "ref_cell_type_col",
                          choices = ref_cell_type_cols,
                          selected = if(length(ref_cell_type_cols) > 0) ref_cell_type_cols[1] else NULL)
        
        updateSelectInput(session, "query_cell_type_col", 
                          choices = query_cell_type_cols,
                          selected = if(length(query_cell_type_cols) > 0) query_cell_type_cols[1] else NULL)
    })
    
    # Update cell type column choices for graph analysis
    observe({
        ref_cols <- colnames(colData(ref_data()))
        query_cols <- colnames(colData(query_data()))
        
        ref_cell_type_cols <- ref_cols
        query_cell_type_cols <- query_cols
        
        updateSelectInput(session, "graph_ref_cell_type_col",
                          choices = ref_cell_type_cols,
                          selected = if(length(ref_cell_type_cols) > 0) ref_cell_type_cols[1] else NULL)
        
        updateSelectInput(session, "graph_query_cell_type_col", 
                          choices = query_cell_type_cols,
                          selected = if(length(query_cell_type_cols) > 0) query_cell_type_cols[1] else NULL)
    })
    
    # Update cell type column choices for anomaly detection
    observe({
        ref_cols <- colnames(colData(ref_data()))
        query_cols <- colnames(colData(query_data()))
        
        ref_cell_type_cols <- ref_cols
        query_cell_type_cols <- query_cols
        
        updateSelectInput(session, "anomaly_ref_cell_type_col",
                          choices = ref_cell_type_cols,
                          selected = if(length(ref_cell_type_cols) > 0) ref_cell_type_cols[1] else NULL)
        
        updateSelectInput(session, "anomaly_query_cell_type_col", 
                          choices = query_cell_type_cols,
                          selected = if(length(query_cell_type_cols) > 0) query_cell_type_cols[1] else NULL)
    })
    
    # Update cell type column choices for discriminant analysis
    observe({
        ref_cols <- colnames(colData(ref_data()))
        query_cols <- colnames(colData(query_data()))
        
        ref_cell_type_cols <- ref_cols
        query_cell_type_cols <- query_cols
        
        updateSelectInput(session, "disc_ref_cell_type_col",
                          choices = ref_cell_type_cols,
                          selected = if(length(ref_cell_type_cols) > 0) ref_cell_type_cols[1] else NULL)
        
        updateSelectInput(session, "disc_query_cell_type_col", 
                          choices = query_cell_type_cols,
                          selected = if(length(query_cell_type_cols) > 0) query_cell_type_cols[1] else NULL)
    })
    
    # Update available cell types (for PCA)
    observe({
        req(input$ref_cell_type_col, input$query_cell_type_col)
        
        ref_types <- unique(colData(ref_data())[[input$ref_cell_type_col]])
        query_types <- unique(colData(query_data())[[input$query_cell_type_col]])
        all_types <- sort(unique(c(ref_types, query_types)))
        all_types <- all_types[!is.na(all_types)]
        
        updateSelectInput(session, "cell_types", choices = all_types)
    })
    
    # Update available cell types for graph analysis
    observe({
        req(input$graph_ref_cell_type_col, input$graph_query_cell_type_col)
        
        ref_types <- unique(colData(ref_data())[[input$graph_ref_cell_type_col]])
        query_types <- unique(colData(query_data())[[input$graph_query_cell_type_col]])
        all_types <- sort(unique(c(ref_types, query_types)))
        all_types <- all_types[!is.na(all_types)]
        
        updateSelectInput(session, "graph_cell_types", choices = all_types)
    })
    
    # Update available cell types for anomaly detection
    observe({
        req(input$anomaly_ref_cell_type_col, input$anomaly_query_cell_type_col)
        
        ref_types <- unique(colData(ref_data())[[input$anomaly_ref_cell_type_col]])
        query_types <- unique(colData(query_data())[[input$anomaly_query_cell_type_col]])
        all_types <- sort(unique(c(ref_types, query_types)))
        all_types <- all_types[!is.na(all_types)]
        
        updateSelectInput(session, "anomaly_cell_types", choices = all_types)
        updateSelectInput(session, "anomaly_cell_type_plot", choices = c("Combined", all_types), selected = "Combined")
    })
    
    # Update available cell types for discriminant analysis
    observe({
        req(input$disc_ref_cell_type_col, input$disc_query_cell_type_col)
        
        ref_types <- unique(colData(ref_data())[[input$disc_ref_cell_type_col]])
        query_types <- unique(colData(query_data())[[input$disc_query_cell_type_col]])
        all_types <- sort(unique(c(ref_types, query_types)))
        all_types <- all_types[!is.na(all_types)]
        
        updateSelectInput(session, "disc_cell_types", choices = all_types)
    })
    
    # Call modular server functions
    setup_overview_outputs(input, output, session, ref_data, query_data)
    setup_pca_outputs(input, output, session, ref_data, query_data)
    setup_discriminant_outputs(input, output, session, ref_data, query_data)
    setup_anomaly_outputs(input, output, session, ref_data, query_data)
    setup_graph_outputs(input, output, session, ref_data, query_data)
}