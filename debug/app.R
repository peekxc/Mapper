library(grapher)
# Define server logic for random distribution app ----
server <- function(input, output) {
    
    output$mapper_graph <- renderUI({ grapher::grapherOutput("grapher") })
    
    output$grapher <- renderGrapher({ G })
    # renderUI({ textOutput("test") })
    # output$grapher <- renderGrapher({ G })
}


# ui <- pageWithSidebar(
#     headerPanel("Simple shiny app with Grapher"),
#     sidebarPanel(),
#     mainPanel( grapherOutput(outputId = "grapher", height = "100vh") )
# )


ui <- htmlTemplate("~/mapper/debug/www/index.html")
#jquery_depend <- file.path(system.file(file.path("dashboard", "www", "js"), package = "Mapper"))

# jq <- htmltools::htmlDependency("jquery", version = "3.3.1", src = jquery_depend, script = "jquery-3.3.1.min.js")

shinyApp(ui = ui, server)
