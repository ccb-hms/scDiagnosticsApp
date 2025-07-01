# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(ggplot2)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater)
library(scran)
library(igraph)
library(cluster)
library(GGally)
library(ggridges)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(shinycssloaders)
library(randomForest)
library(corpcor)
library(isotree)
library(transport)
library(Matrix)
library(methods)

# Main application file
source("global.R")
source("ui/ui_main.R")
source("server/server_main.R")

# Run the application
shinyApp(ui = ui, server = server)

