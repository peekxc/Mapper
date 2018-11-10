## server.R
## Sets up the shiny server function. Assumes upon sourcing that a valid MapperRef object is defined in the calling
## environment associated with the symbol 'M'
if (!exists("M")){ stop("Error creating server; there was no 'MapperRef' object associated with the name 'M'.") }
if (!exists("G")){ stop("Error creating server; there was no 'grapher' object associated with the name 'G'.") }
if (!exists("X")){ stop("Error creating server; there was no 'data.frame' object associated with the name 'X'.") }
if (!exists("color_funcs")){ stop("Error creating server; there was no 'color_funcs' object.") }

# Define server logic
server <- function(input, output, session) {

  rv <- reactiveValues(grapher_widget = G)
  output$mapper_network_ui <- renderUI({ grapher::grapherOutput(outputId = "grapher") })
  output$grapher <- grapher::renderGrapher({ rv$grapher_widget })
  
  ## ========== Sources ==========
  

  ## ========= Endpoints =========

  ## ========= Observers =========

  # Observer to center the mapper
  observeEvent(input$recenter, {
    print("Centering grapher")
    grapher::center(id = "grapher")
  })

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

  # Node min/max size observer
  observeEvent({ input$node_min; input$node_max; input$scale }, {
    req(input$scale)
    # browser()
    if (toupper(input$scale) == "LINEAR"){
      scale_f <- function(x){ return(x) }
    } else if (toupper(input$scale) == "LOGARITHMIC"){
      scale_f <- function(x){ return(log(x)) }
    }
    print(sprintf("Adjusting node scale: min = %d, max = %d", input$node_min, input$node_max))
    if (length(M$vertices) > 0){
      vertex_sizes <- sapply(M$vertices, length)
      vertex_radii <- (input$node_max - input$node_min)*normalize(scale_f(vertex_sizes)) + input$node_min
      setNodeSize(id = "grapher", r = as.integer(vertex_radii))
    }
  })

  # Node color observer
  observeEvent(input$colorFunc, {
    #browser()
    f <- color_funcs[[input$colorFunc]]
    res <- f(M, X)
    binned_color <- bin_color(res, output_format = 'hex9')
    setNodeColor(id = "grapher", color = binned_color)
  })

  # Force input change
  observeEvent({ input$charge; input$gravity; input$link_strength; input$link_distance }, {
    grapher::enableForce(id = "grapher")
    # G$updateForce(id = "grapher_canvas",
    #               charge = list(strength=input$charge),
    #               gravity = list(strength=input$gravity),
    #               link = list(strength=input$link_strength, distance=input$link_distance))
  })

  # Simple row selecter for when nodes are selected
  proxy <- DT::dataTableProxy('mytable')
  observeEvent({ input$node_selected; input$selection_mode }, {
    req({ input$node_selected; input$selection_mode })
    print(sprintf("Node selected: %d (mode = %s)", input$node_selected, input$selection_mode))
    if (!is.null(input$node_selected) && input$selection_mode == "TRUE"){
      if (input$node_selected != -1 && input$node_selected < length(M$vertices)){
        # DT::selectRows(proxy, M$G$nodes[[input$node_selected]])
        DT::replaceData(proxy, X[M$vertices[[input$node_selected+1]],])
      } else {
        DT::replaceData(proxy, X)
        # DT::reloadData(proxy, clearSelection = "row")
      }
    }
  })
}
