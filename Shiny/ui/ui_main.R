# Load required libraries for UI
source("ui/ui_overview.R", local = TRUE)
source("ui/ui_pca.R", local = TRUE)
source("ui/ui_discriminant.R", local = TRUE)
source("ui/ui_graph.R", local = TRUE)
source("ui/ui_anomaly.R", local = TRUE)
source("ui/ui_wasserstein.R", local = TRUE)
source("ui/ui_marker.R", local = TRUE)
source("ui/ui_dimred.R", local = TRUE)
source("ui/ui_upload.R", local = TRUE)

# Define UI
ui <- dashboardPage(
    dashboardHeader(title = "scDiagnostics Interactive App"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Data Overview", tabName = "overview", icon = icon("table")),
            menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
            menuItem("PCA Projection", tabName = "pca", icon = icon("chart-line")),
            menuItem("Discriminant Space", tabName = "discriminant", icon = icon("search-plus")),
            menuItem("Graph Integration", tabName = "graph", icon = icon("project-diagram")),
            menuItem("Anomaly Detection", tabName = "anomaly", icon = icon("exclamation-triangle")),
            menuItem("Wasserstein Distance", tabName = "wasserstein", icon = icon("wave-square")),
            menuItem("Marker Expression", tabName = "marker", icon = icon("dna")),
            menuItem("Gene Expression DimRed", tabName = "dimred", icon = icon("project-diagram"))
        )
        # REMOVED THE DATASET SELECTION FROM SIDEBAR
    ),
    
    dashboardBody(
        # Add custom CSS for better styling
        tags$head(
            tags$style(HTML("
                .content-wrapper, .right-side {
                    background-color: #f4f4f4;
                }
                .box {
                    border-radius: 5px;
                }
                .nav-tabs-custom > .nav-tabs > li.active {
                    border-top-color: #3c8dbc;
                }
            "))
        ),
        
        tabItems(
            overview_tab,
            upload_tab,
            pca_tab,
            discriminant_tab,
            graph_tab,
            anomaly_tab,
            wasserstein_tab,
            marker_tab,
            dimred_tab
        )
    )
)