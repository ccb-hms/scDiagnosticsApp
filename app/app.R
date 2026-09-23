# Temporary placeholder.
#
# This exists only to give the Connect content one successful deployment, so
# that environment variables can be set on it. Connect refuses to accept
# variables on content that has never published, and the real app cannot
# publish until GITHUB_PAT is fixed. Deploy this, set the variable, then
# revert the commit that introduced it.

library(shiny)

ui <- fluidPage(
    titlePanel("scDiagnostics"),
    p("This application is being set up. Please check back shortly.")
)

server <- function(input, output, session) {}

shinyApp(ui, server)
