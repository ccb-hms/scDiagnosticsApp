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
            menuItem("Anomaly Detection", tabName = "anomaly", icon = icon("exclamation-triangle")),
            menuItem("Graph Integration", tabName = "graph", icon = icon("project-diagram")),
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
            anomaly_tab,
            graph_tab,
            wasserstein_tab,
            marker_tab,
            dimred_tab
        ),

        # Citation footer
        tags$footer(
            style = "margin-top: 20px; padding: 10px; text-align: center; font-size: 12px; color: #666; border-top: 1px solid #d2d6de;",
            HTML(paste0(
                "If you use this app in published research, please cite: Christidis A, Ghazi A, Chawla S, ",
                "Turaga N, Gentleman R, Geistlinger L (2026). scDiagnostics: systematic assessment of cell type ",
                "annotation in single-cell transcriptomics data. <em>Briefings in Bioinformatics</em>, 27(5), bbag496. ",
                "doi: <a href='https://doi.org/10.1093/bib/bbag496' target='_blank'>10.1093/bib/bbag496</a>."
            ))
        )
    )
)