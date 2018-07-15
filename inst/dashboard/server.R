## server.R
## Sets up the shiny server function. Assumes upon sourcing that a valid MapperRef object is defined in the calling
## environment associated with the symbol 'M'
if (!exists("M")){ stop("Error creating server; there was no 'MapperRef' object associated with the name 'M'.") }
if (!exists("G")){ stop("Error creating server; there was no 'grapher' object associated with the name 'G'.") }
if (!exists("X")){ stop("Error creating server; there was no 'data.frame' object associated with the name 'X'.") }
if (!exists("color_funcs")){ stop("Error creating server; there was no 'color_funcs' object.") }
if (!all.equal(G$mapper_obj, M)){ stop("grapher instance not synchronized") }

# Define server logic
server <- function(input, output, session) {

  ## Reactive to get a updated mapper; only needed when conductive/reactive-based dependencies need to be explicitly constructed
  # getMapper <- reactive({ M })

  ## ========== Sources ==========
  source(system.file(file.path("dashboard", "reactivity", "sources.R"), package = "Mapper"), local = TRUE)

  ## ========= Endpoints =========
  source(system.file(file.path("dashboard", "reactivity", "endpoints.R"), package = "Mapper"), local = TRUE)

  ## ========= Observers =========

  # Observer to center the mapper
  observeEvent(input$recenter, {
    print("Centering grapher")
    G$center(id = "grapher_canvas")
  })

  # Mapper cover parameter observer
  observeEvent({ input$num_intervals; input$overlap }, {
    M$cover$setResolution(input$num_intervals)
    M$cover$setOverlap(input$overlap)
    M$update()
    G$setDefaultJsonConfig() # update json configuration with new mapper
    cat(sprintf("Num edges: %d\n", sum(M$G$adjacency == 1)/2))
    G$updateNetwork("grapher_canvas", net = G$exportJSON())
    updateFilterPlot(M)
  })

  # Node min/max size observer
  observeEvent({ input$node_min; input$node_max }, {
    req(input$scale, input$colorFunc)
    print(sprintf("Adjusting node scale: min = %d, max = %d", input$node_min, input$node_max))
    G$updateNodeSize(id = "grapher_canvas", node_min = input$node_min, node_max = input$node_max)
  })

  # Node color observer
  observeEvent(input$colorFunc, {
    f <- color_funcs[[input$colorFunc]]
    color_res <- f(M, X)
    rbw_pal <- rev(rainbow(100, start = 0, end = 4/6))
    binned_idx <- cut(color_res, breaks = 100, labels = F)
    binned_color <- substr(rbw_pal[binned_idx], start = 0, stop = 7)
    G$updateNodeColor(id = "grapher_canvas", color = binned_color)
  })

  # Node scaling observer
  observeEvent(input$scale, {
    if (toupper(input$scale) == "LINEAR"){
      G$json_config$node_scale <- function(x){ return(x) }
    } else if (toupper(input$scale) == "LOGARITHMIC"){
      G$json_config$node_scale <- function(x){ return(log(x)) }
    }
    G$updateNodeSize(id = "grapher_canvas", node_min = input$node_min, node_max = input$node_max)
  })

  # Force input change
  observeEvent({ input$charge; input$gravity; input$link_strength; input$link_distance }, {
    G$updateForce(id = "grapher_canvas",
                  charge = list(strength=input$charge),
                  gravity = list(strength=input$gravity),
                  link = list(strength=input$link_strength, distance=input$link_distance))
  })

  # Simple row selecter for when nodes are selected
  proxy <- DT::dataTableProxy('mytable')
  observeEvent({ input$node_selected; input$selection_mode }, {
    req({ input$node_selected; input$selection_mode })
    print(sprintf("Node selected: %d (mode = %s)", input$node_selected, input$selection_mode))
    if (!is.null(input$node_selected) && input$selection_mode == "TRUE"){
      if (input$node_selected != -1 && input$node_selected < length(M$G$nodes)){
        # DT::selectRows(proxy, M$G$nodes[[input$node_selected]])
        DT::replaceData(proxy, X[M$G$nodes[[input$node_selected+1]],])
      } else {
        DT::replaceData(proxy, X)
        # DT::reloadData(proxy, clearSelection = "row")
      }
    }
  })
}
