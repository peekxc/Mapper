## global.R
## Attaches the default sources to the Mapper dashboard app.
## Author: Matt Piekenbrock

## Normalize between 0-1, unless all the same
normalize <- function(x) { 
  if (all(x == x[1])){ return(rep(1, length(x))) }
  else {  (x - min(x))/(max(x) - min(x)) }
}

mapper_graph <- reactive({
  am <- M$.simplicial_complex$as_adjacency_matrix()
  igraph::graph_from_adjacency_matrix(am, mode = "undirected", add.colnames = NA) 
})

## Makes the default configuration depend on M
json_config <- reactive({
  MG <- mapper_graph()
  grapher::getDefaultJsonConfig(network=MG)
})

local({
  ## The current environment passed in
  c_env <- environment()
  
  ## Computes the mean filter value for each node
  c_env[["Mean filter value"]] <- function(M, ...){
    agg_pt_fv <- sapply(M$vertices, function(n_idx){ apply(as.matrix(M$cover$filter_values[n_idx,]), 1, mean)})
    agg_node_val <- sapply(agg_pt_fv, mean)
    return(agg_node_val)
  }
  
  ## Computes the mean data value for each node
  c_env[["Mean data value"]] <- function(M, ...){
    agg_pt_x <- sapply(M$vertices, function(n_idx){ apply(as.matrix(M$X[n_idx,]), 1, mean) })
    agg_node_x <- sapply(agg_pt_x, mean)
    return(agg_node_x)
  }
  
  ## Density of nodes
  c_env[["Density"]] <- function(M, ...){
    return(sapply(M$vertices, length))
  }
  rm(c_env)
})

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

## Utility function to make a fancy datatables object with all the cool extensions
fancytable <- function(dt){
  DT::datatable(dt, fillContainer = TRUE,
                extensions = c("Scroller", "Buttons", "ColReorder", "FixedColumns", "Scroller"),
                options = list(dom = 'Bfrtip', buttons = I('colvis'), scrollX = TRUE, colReorder = TRUE, fixedColumns = TRUE, deferRender = TRUE, scrollY = 200, scroller = TRUE)
  )
}

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



