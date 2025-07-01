# Set maximum file upload size to 100MB
options(shiny.maxRequestSize = 100*1024^2)

# Load example datasets
query_cells <- readRDS("data/query_marrow_myeloid.rds")
reference_cells <- readRDS("data/reference_marrow_myeloid.rds")
reference_cells_subset <- readRDS("data/reference_subset_marrow_myeloid.rds")

# Source support functions
source("support/argumentCheck.R")
source("support/generateColors.R")
source("support/projectPCA.R")
source("support/calculateVarImpOverlap.R")
source("support/plotCellTypePCA_app.R")
source("support/calculateGraphIntegration_app.R")
source("support/plotGraphIntegration_app.R")
source("support/detectAnomaly_app.R")
source("support/plotAnomaly_app.R")
source("support/calculateDiscriminantSpace_app.R")
source("support/plotDiscriminantSpace_app.R")
source("support/calculateWassersteinDistance_app.R")
source("support/plotWassersteinDistance_app.R")
source("support/plotMarkerExpression_app.R")
source("support/plotGeneExpressionDimred_app.R")