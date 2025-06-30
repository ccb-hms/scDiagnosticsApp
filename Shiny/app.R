# Main application file
source("global.R")
source("ui/ui_main.R")
source("server/server_main.R")

# Run the application
shinyApp(ui = ui, server = server)
