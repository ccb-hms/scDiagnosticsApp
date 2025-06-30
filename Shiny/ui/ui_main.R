# Source individual UI components
source("ui/ui_overview.R", local = TRUE)
source("ui/ui_pca.R", local = TRUE)
source("ui/ui_discriminant.R", local = TRUE)
source("ui/ui_graph.R", local = TRUE)
source("ui/ui_anomaly.R", local = TRUE)
source("ui/ui_placeholder.R", local = TRUE)

# Define main UI
ui <- dashboardPage(
    dashboardHeader(title = "scDiagnostics Interactive App"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Data Overview", tabName = "overview", icon = icon("table")),
            menuItem("PCA Projection", tabName = "pca", icon = icon("project-diagram")),
            menuItem("Discriminant Projection", tabName = "discriminant", icon = icon("sitemap")),
            menuItem("Anomaly Detection", tabName = "anomaly", icon = icon("exclamation-triangle")),
            menuItem("Graph Integration", tabName = "graph", icon = icon("network-wired")),
            menuItem("Cell Type Distribution", tabName = "distribution", icon = icon("chart-pie")),
            menuItem("Marker Gene Analysis", tabName = "markers", icon = icon("dna")),
            menuItem("Quality Metrics", tabName = "quality", icon = icon("chart-line"))
        )
    ),
    
    dashboardBody(
        tabItems(
            overview_tab,
            pca_tab,
            discriminant_tab,
            graph_tab,
            anomaly_tab,
            distribution_tab,
            markers_tab,
            quality_tab
        )
    )
)