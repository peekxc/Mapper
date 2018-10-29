library(shiny)


ui <- fluidPage(
    titlePanel("Multiscale Mapper demo"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("number_intervals", "Resolution:", min = 1, max = 50, value = 1), 
            sliderInput("percent_overlap", "Gain:", min = 0, max = 1, value = 0)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           grapherOutput("grapher")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$grapher <- renderGrapher({
        m$plot_interactive()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
