## Re-useable shiny modules for interacting with a Mapper visualization.  

#' @title mapperVis - shiny UI module 
#' @description Module for configuring the visualization options for Mapper.
#' @export
mapperVisOutput <- function(id, output_type=c("grapher")){
  ns <- NS(id)
  if (output_type == "grapher"){
    # browser()
    grapher::grapherOutput(outputId = ns("widget_out"), width = "100%", height = "100%")
  }
}

#' @title mapperVis - shiny server module
#' @description Returns a set of reactiveValues for configuring the mapper visualization
#' @export
mapperVis <- function(input, output, session, mapper, X, output_type){
  # browser()
  rv <- reactiveValues()
  rv$output_type <- isolate(output_type())
  rv$mapper <- isolate(mapper())
  rv$X <- isolate(X())
  rv$vertex_sizes <- NULL
  rv$selected_points <- NULL
    
  ## Render the mapper
  # browser()
  if (isolate(rv$output_type) == "grapher"){
    rv$output_id <- session$ns("widget_out")
    output[["widget_out"]] <- grapher::renderGrapher({
      rv$mapper$as_grapher(construct_widget = TRUE)
    })
  }
  return(rv) ## Return the reactive values
}

#' @title Node sizing input/UI module
#' @export
nodeSizingInput <- function(id, node_min_rng = c(1, 25), node_max_rng = c(1, 50), node_scaling = c("linear", "logarithmic")) {
  ns <- NS(id)
  tagList(
    numericInput(inputId = ns("node_min"), label = "Node min size", value = diff(node_min_rng)/2, min = node_min_rng[[1]], max = node_min_rng[[2]], width = '100%', step = 1L),
    numericInput(inputId = ns("node_max"), label = "Node max size", value = diff(node_max_rng)/2, min = node_max_rng[[1]], max = node_max_rng[[2]], width = '100%', step = 1L),
    selectInput(inputId = ns("scale"), label = "Node scaling", selected = "Logarithmic", choices = c("Linear", "Logarithmic"), width = "100%")
  )
}

#' @title Node sizing server function 
#' @export
nodeSizing <- function(input, output, session, mapper_vis, v_f){
  ## Get reactive node sizes 
  normalize <- function(x) { (x - min(x))/(diff(range(x))) }
  vertex_radii <- reactive({
    scale_f <- switch (toupper(input$scale), "LINEAR" = function(x){ x },  "LOGARITHMIC" = function(x){ log(x) } )
    if (length(mapper_vis$mapper$vertices) > 0){
      vertex_values <- unname(v_f()(mapper_vis$mapper))
      return((input$node_max - input$node_min)*normalize(scale_f(vertex_values)) + input$node_min)
    } else { return(NULL) }
  })

  ## Register observers
  observeEvent({ input$node_min; input$node_max; input$scale }, {
    req(input$scale)
    if (mapper_vis$output_type == "grapher"){
      # print(mapper_vis$output_id)
      grapher::setNodeSize(mapper_vis$output_id, vertex_radii())
    }
  })
  
  ## Return node sizes for other uses
  return(vertex_radii)
}

#' @title Node coloring input/UI module
#' @export
nodeColorInput <- function(id, color_f){
  stopifnot(is.list(color_f))
  stopifnot(all(sapply(color_f, is.function)))
  ns <- NS(id)
  f_names <- names(color_f)
  selectInput(ns("node_color_f"), label = "node_color_f", choices = f_names, selected = f_names[[1]], width = "100%")
}

#' @title Node coloring associated server module
#' @export
nodeColor <- function(input, output, session, mapper_vis, color_f){
  # browser()
  ## Bin the node colors 
  v_colors <- reactive({ bin_color(color_f()[[input$node_color_f]](mapper_vis$mapper, mapper_vis$X)) })
  
  ## Register observers
  observeEvent(input$node_color_f, {
    req(input$node_color_f)
    if (mapper_vis$output_type == "grapher"){
      grapher::setNodeColor(id = mapper_vis$output_id, color = v_colors())
    }
  })
  return(v_colors)
}

# nodeHighlightInput <- function(id, ){
#   ns <- NS(id)
#   ns("highlighted")
# }

#' @title nodeHighlight shiny server module
#' @description Server-only function to 'highlight' nodes. This does the equivalent of changing the nodes colors, 
#' but adjusts low values to have lower opacity, and registers a reactive expression to to reset when 
#' no data points are selected. Takes as input a boolean vector indicated which _points_ should be highlighted. 
#' @export
nodeHighlight <- function(mapper_vis, reset){
  ## Register observers
  observeEvent(mapper_vis$selected_points, { ## Enact the highlighting
    cdata <- mapper_vis$selected_points ## data.frame w/ boolean vector selected_ as column
    v_col <- switch(as.character(all(cdata[["selected_"]] == FALSE)), 
      "TRUE" = reset(), 
      "FALSE" = {
        selected_idx <- which(cdata[["selected_"]])
        v_weight <- sapply(mapper_vis$mapper$vertices, function(v) sum(v %in% selected_idx))
        bin_color(v_weight, col_pal = viridis::inferno(100, begin = 0.15, end = 0.8, alpha = seq(0.65, 1, length.out = 100)))
    })
    
    if (mapper_vis$output_type == "grapher"){
      grapher::setNodeColor(id = mapper_vis$output_id, color = v_col) 
    }
  })
}


#' @title Mapper input options - shiny UI module
#' @export
mapperOptionsInput <- function(id, mapper){
  # shiny::numericInput()
}

#' @title Mapper input options - shiny server module
#' @export
mapperOptions <- function(input, output, session, mapper){
  # Mapper cover parameter observer
  observeEvent({ input$num_intervals; input$overlap }, {
    # browser()
    req(input$num_intervals, input$overlap)
    M$cover$number_intervals <- input$num_intervals
    M$cover$percent_overlap <- input$overlap
    M$cover$construct_cover()
    M$compute_k_skeleton(1L)
    rv$grapher_widget <- M$as_grapher(construct_widget = TRUE)
    output$grapher <- grapher::renderGrapher({ rv$grapher_widget })
  })
}

#' @title brushScatterOutput - shiny UI module
#' @export
brushScatterOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns("plot"), brush = ns("brush"))
}

#' @title brushScatter - shiny server module
#' @export
brushScatter <- function(input, output, session, mapper_vis, X, xvar, yvar, reset) {
  
  ## Yields the data frame with an additional column "selected_" that indicates whether that observation is brushed
  observeEvent(input$brush, {
    req(input$brush)
    mapper_vis$selected_points <- brushedPoints(X(), input$brush, xvar = xvar(), yvar = yvar(), allRows = TRUE)
  })
  
  ## Show which points are highlighted
  output$plot <- renderPlot({
    cdata <- mapper_vis$selected_points ## data.frame w/ boolean vector selected_ as column
    if (is.null(cdata) || is.data.frame(cdata) && all(cdata[["selected_"]] == FALSE)){
      pt_color <- adjustcolor("black", alpha.f = 0.85)
    } else {
      pt_color <- ifelse(cdata[["selected_"]], adjustcolor("black", alpha.f = 0.85), adjustcolor("black", alpha.f = 0.15))
    }
    plot(X()[, c(xvar(), yvar())], pch = 20, col = pt_color, cex = 0.8)
  })
  
  ## Attach brush observer
  # browser()
  # observeEvent(input$brush, {
  #   selected <- which(brushedPoints(X(), input$brush, xvar = xvar(), yvar = yvar(), allRows = TRUE)[["selected_"]])
  #   v_col <- switch(as.character(all(selected == FALSE)), 
  #     "TRUE" = reset(), 
  #     "FALSE" = {
  #       v_weight <<- sapply(mapper_vis$mapper$vertices, function(v) sum(v %in% selected))
  #       highlight_pal <<- viridis::inferno(100, begin = 0.15, end = 0.8, alpha = seq(0.65, 1, length.out = 100))
  #       bin_color(v_weight, col_pal = highlight_pal)
  #     })
  #   if (mapper_vis$output_type == "grapher"){
  #     grapher::setNodeColor(id = mapper_vis$output_id, color = v_col) 
  #   }
  # })
  # return(dataWithSelection)
}

#' @title linkedTableOutput - shiny UI module
#' @export
linkedTableOutput <- function(id){
  ns <- NS(id)
  DT::DTOutput(ns("data_table")) 
}

#' @title linkedTableOutput - shiny server module
#' @description Utility function to make a datatables object with nice default extensions.
#' @export
linkedTable <- function(input, output, session, dt, config=NULL, options=NULL){
  
  ## Make default config 
  def_cfg <- list(data = dt, fillContainer = TRUE, class = "compact cell-border stripe",
                  extensions = c("Scroller", "Buttons", "ColReorder", "FixedColumns", "Scroller"))
  def_opt <- list(dom = 'Bfrtip', buttons = I('colvis'), scrollX = TRUE, colReorder = TRUE,
                  fixedColumns = TRUE, deferRender = TRUE, scrollY = 200, scroller = TRUE)
  if (!missing(config) && !is.null(config)){ for (vn in names(config)){ def_cfg[[vn]] <- config[[vn]] } }
  if (!missing(options) && !is.null(options)){ for (vn in names(options)){ def_opt[[vn]] <- options[[vn]] } }
  def_cfg[["options"]] <- def_opt
  
  ## Make the reactive datatable expression
  datatable <- reactive({ do.call(DT::datatable, def_cfg) })
  
  ## Render the data table
  output$data_table <- DT::renderDataTable({ datatable() })
  
  ## Attach observers
  base_id <- "data_table"
  observeEvent(input[[paste0(base_id, "_rows_selected")]], {
    ## do something
    print("hello")
  })
  return(datatable)
}


# api <- list(highlight=NULL, set_node_color=NULL, set_edge_color=NULL, set_node_size=NULL)
# 
# highlight <- function(){
#   grapher::setNodeSize(id = grapher_id, r = as.integer(vertex_radii))
# }
# if (output_type == "grapher"){
#   if (!missing(highlight_nodes)){
#     observeEvent(highlight_nodes()){
#       highlighted <- highlight_nodes()
#       if (all(highlighted == FALSE)){ print("reset") }
#       v_weight <- sapply(mapper$vertices, function(v) sum(v %in% highlighted))
#       v_col <- bin_color(v_weight, col_pal = viridis::inferno(100, begin = 0.15, end = 0.8, alpha = seq(0.65, 1, length.out = 100)))
#       grapher::setNodeColor(id = "grapher", color = v_col)
#     }
#   } 
#   if (!missing(set_node_size)){
#     observeEvent(set_node_size()){
#       node_sizes <- set_node_size()
#     }
#   }
#   if (!missing(set_node_color)){
#     observeEvent(set_node_color()){
#       
#     }
#   }
#   if (!missing(set_edge_color)){
#     
#   }
# 
# }
# grapher::center()
# data_xy <- reactive({
#   d <- dim(data())[[2]]
#   if (d == 1){ 
#     x <- data()
#     new_data <- structure(as.data.frame(cbind(x, runif(length(x)))), names = c("V1", "V2")) 
#     return(list(data=new_data, x = "V1", y = "V2"))
#   }
#   else if (d > 1){ return(list(data=data(), x = x(), y = y()))}
# })  d <- dim(M$cover$filter_values)[[2]]
# if (d == 1){
#   fv.brush <- brushedPoints(fv, input$panel2_brush, xvar = cn[[1]], yvar = cn[[2]], allRows = TRUE)
#   plot(fv, ylim = c(0, 1), pch = 20, col=ifelse(fv.brush$selected_, adjustcolor("black", alpha.f = 0.90), adjustcolor("black", alpha.f = 0.45)), 
#        ylab = "random jitter", xlab = "Filter value")
# } else if (d == 2){
#   plot(M$cover$filter_values, pch = 20)
