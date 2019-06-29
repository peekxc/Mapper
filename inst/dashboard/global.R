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
  requireNamespace("grapher")
  grapher::getDefaultJsonConfig(network=MG)
})

## Utility function to make a fancy datatables object with all the cool extensions
fancytable <- function(dt){
  DT::datatable(dt, fillContainer = TRUE, class = "compact cell-border stripe",
                extensions = c("Scroller", "Buttons", "ColReorder", "FixedColumns", "Scroller"),
                options = list(dom = 'Bfrtip', buttons = I('colvis'), scrollX = TRUE, colReorder = TRUE, fixedColumns = TRUE, deferRender = TRUE, scrollY = 200, scroller = TRUE)
  )
}

# Call this function with an input (such as `textInput("text", NULL, "Search")`) to add a list of dropdown inputs to the navbar
dropDownItems <- function(inputs) {
  for(i in 1:length(inputs)){
    input_el <- inputs[[i]]
    input_el$attribs$class <- paste(input_el$attribs$class, "dropdown_item")
    inputs[[i]] <- input_el
  }
  return(inputs)
}



