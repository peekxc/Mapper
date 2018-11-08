## endpoints.R
## Attaches the default endpoints to the Mapper dashboard app.
## Author: Matt Piekenbrock

## Number of intervals
output$num_intervals <- renderUI({
  numericInput(inputId = "num_intervals", label = "Number of Intervals", value = head(M$cover$number_intervals, 1), min = 1, step = 1, width = '100%')
})

## Output parameter slider
output$overlap <- renderUI({
  sliderInput(inputId = "overlap",
              label = "Overlap Percentage",
              min = 0, max = 100, value = head(M$cover$percent_overlap, 1), round = -2, width = '100%')
})

## Generate node min input
output$node_min <- renderUI({
  numericInput(inputId = "node_min", label = "Node min", value = dash_config$node_min, min = 1, width = '100%', step = 1L)
})

## Generate node max input
output$node_max <- renderUI({
  numericInput(inputId = "node_max", label = "Node max", value = dash_config$node_max, min = 1, width = '100%', step = 1L)
})

## Recenter network
output$recenter <- renderUI({
  shiny::actionButton(inputId = "recenter", label = "Center graph", width = "100%")
})

## Color functions
output$colorFunc <- renderUI({
  color_f_choices <- ls(color_funcs)
  shiny::selectInput(inputId = "colorFunc", label = "Coloring Function",
                     choices = color_f_choices, width = "100%")
})

## Node scaling function
output$scale <- renderUI({
  shiny::selectInput(inputId = "scale", label = "Scaling Function", selected = "Logarithmic",
                     choices = c("Linear", "Logarithmic"), width = "100%")
})

output$selection_mode <- renderUI({
  shiny::radioButtons(inputId = "selection_mode", label = "Filter on Selection", choices = c("TRUE", "FALSE"), selected = "TRUE", width = "100%")
})


## Mapper reactive network panel
# output$mapper_network_ui <- renderUI({ grapher::grapherOutput(outputId = "grapher") })
# output$grapher <- grapher::renderGrapher({ G })

## Panel 1 - Interactive Table
output$panel1 <- renderUI({ DT::DTOutput(outputId = "mdatatable") })
output$mdatatable <- DT::renderDT({
  DT::datatable(X, fillContainer = TRUE,
                extensions = c("Scroller", "Buttons", "ColReorder", "FixedColumns", "Scroller"),
                options = list(dom = 'Bfrtip', buttons = I('colvis'), scrollX = TRUE, colReorder = TRUE, fixedColumns = TRUE, deferRender = TRUE, scrollY = 200, scroller = TRUE))
})
## Panel 2 - Filter plot
# updateFilterPlot <- function(M){
#   output$panel2 <- renderPlot({ M$cover$plotFilterSpace(show_ls_bounds = TRUE, show_lsfi = TRUE) },
#                               height = function() { session$clientData[["output_panel2_height"]] },
#                               width = function(){ session$clientData[["output_panel2_width"]] }
#   )
# }
# updateFilterPlot(M)

## Network force options

# Call this function with an input (such as `textInput("text", NULL, "Search")`) if you
# want to add an input to the navbar
dropDownItems <- function(inputs) {
  for(i in 1:length(inputs)){
    input_el <- inputs[[i]]
    input_el$attribs$class <- paste(input_el$attribs$class, "dropdown_item")
    inputs[[i]] <- input_el
  }
  inputs
}

output$force_opts <- renderUI({
  dropDownItems(list(
    shiny::numericInput(inputId = "charge", label = "Charge", value = -50, step = 10),
    shiny::numericInput(inputId = "gravity", label = "Gravity", value = 0, step = 0.001),
    shiny::numericInput(inputId = "link_strength", label = "Link Strength", value = 0.2, step = 0.1),
    shiny::numericInput(inputId = "link_distance", label = "Link Distance", value = 60, step = 25)
  ))
})
