## server.R
## Sets up the shiny server function. Assumes upon sourcing that a valid MapperRef object is defined in the calling
## environment associated with the symbol 'M'
if (!exists("M", envir = .GlobalEnv)){ stop("Error creating server; there was no 'MapperRef' object associated with the name 'M'.") }
if (!exists("G", envir = .GlobalEnv)){ stop("Error creating server; there was no 'grapher' object associated with the name 'G'.") }
if (!exists("X", envir = .GlobalEnv)){ stop("Error creating server; there was no 'data.frame' object associated with the name 'X'.") }
if (!exists("color_funcs", envir = .GlobalEnv)){ stop("Error creating server; there was no 'color_funcs' object.") }

# load_env <- function(from, to){
#   symbols <- grep("\\..*", ls(from), ignore.case = TRUE, value = TRUE)
#   for (sym in symbols) {
#     new_sym <- strsplit(sym, split = "\\.")[[1]][2]
#     assign(x = new_sym, value = from[[sym]], envir = to)
#   }
# }

# Define server logic
server <- function(input, output, session) {

  ## ========== Sources ==========
  rv <- reactiveValues(grapher_widget = G)
  d <- dim(M$cover$filter_values)[[2]]
  if (d == 1){ ## d == 1, 
    fv <- as.data.frame(cbind(M$cover$filter_values, runif(nrow(M$cover$filter_values))))
    fv_cn <- colnames(fv)
  } else if (d == 2){
    fv <- as.data.frame(cbind(M$cover$filter_values, runif(nrow(M$cover$filter_values))))
    fv_cn <- colnames(fv)
  }

  ## ========= Endpoints =========
  output$num_intervals <- renderUI({ ## Number of intervals
    numericInput(inputId = "num_intervals", label = "Number of Intervals", value = head(M$cover$number_intervals, 1), min = 1, step = 1, width = '100%')
  })
  output$overlap <- renderUI({ ## Output parameter slider
    sliderInput(inputId = "overlap",
                label = "Overlap Percentage",
                min = 0, max = 100, value = head(M$cover$percent_overlap, 1), round = -2, width = '100%')
  })
  output$node_min <- renderUI({ ## Generate node min input
    numericInput(inputId = "node_min", label = "Node min", value = dash_config$node_min, min = 1, width = '100%', step = 1L)
  })
  output$node_max <- renderUI({ ## Generate node max input
    numericInput(inputId = "node_max", label = "Node max", value = dash_config$node_max, min = 1, width = '100%', step = 1L)
  })
  output$recenter <- renderUI({  ## Recenter network
    shiny::actionButton(inputId = "recenter", label = "Center graph", width = "100%")
  })
  output$colorFunc <- renderUI({ ## Color functions
    color_f_choices <- ls(color_funcs)
    shiny::selectInput(inputId = "colorFunc", label = "Coloring Function",
                       choices = color_f_choices, width = "100%")
  })
  output$scale <- renderUI({ ## Node scaling function
    shiny::selectInput(inputId = "scale", label = "Scaling Function", selected = "Logarithmic",
                       choices = c("Linear", "Logarithmic"), width = "100%")
  })
  output$selection_mode <- renderUI({
    shiny::radioButtons(inputId = "selection_mode", label = "Filter on Selection", choices = c("TRUE", "FALSE"), selected = "TRUE", width = "100%")
  })
  output$force_opts <- renderUI({ ## force options
    dropDownItems(list(
      shiny::numericInput(inputId = "charge", label = "Charge", value = -50, step = 10),
      shiny::numericInput(inputId = "gravity", label = "Gravity", value = 0, step = 0.001),
      shiny::numericInput(inputId = "link_strength", label = "Link Strength", value = 0.2, step = 0.1),
      shiny::numericInput(inputId = "link_distance", label = "Link Distance", value = 60, step = 25)
    ))
  })
  
  ## Top panel - Mapper reactive network panel
  output$mapper_network_ui <- renderUI({ grapher::grapherOutput(outputId = "grapher") })
  output$grapher <- grapher::renderGrapher({ rv$grapher_widget })
  
  ## Bottom panel 1 - Interactive Table
  output$panel1 <- renderUI({ DT::DTOutput(outputId = "mdatatable") })
  output$mdatatable <- DT::renderDT({ fancytable(X) })
  
  ## Bottom panel 2 - Multiple plot outputs
  output$panel2 <- renderUI({ 
    tabsetPanel(
      tabPanel("Filter Values", plotOutput("filter_plot", brush = brushOpts(id = "panel2_brush", fill = "#ccc", direction = "x"))),
      tabPanel("PP", plotOutput("pp_plot"))
    )
  })
  output$filter_plot <- renderPlot({
    d <- dim(M$cover$filter_values)[[2]]
    if (d == 1){
      fv.brush <- brushedPoints(fv, input$panel2_brush, xvar = cn[[1]], yvar = cn[[2]], allRows = TRUE)
      plot(fv, ylim = c(0, 1), pch = 20, col=ifelse(fv.brush$selected_, adjustcolor("black", alpha.f = 0.90), adjustcolor("black", alpha.f = 0.45)), 
           ylab = "random jitter", xlab = "Filter value")
    } else if (d == 2){
      plot(M$cover$filter_values, pch = 20)
    }
    
  })
  
  
  
  ## ========= Observers =========
  observeEvent(input$panel2_brush, {
    req(input$panel2_brush)
    # browser()
    bp <- brushedPoints(fv, input$panel2_brush, xvar = cn[[1]], yvar = cn[[2]], allRows = TRUE)
    highlighted <- which(bp$selected_)
    v_weight <- sapply(M$vertices, function(v) sum(v %in% highlighted))
    v_col <- bin_color(v_weight, col_pal = viridis::inferno(100, begin = 0.15, end = 0.8, alpha = seq(0.65, 1, length.out = 100)))
    grapher::setNodeColor(id = "grapher", color = v_col)
    # print("brushing clicked")
  })
  
  observeEvent(input$recenter, { grapher::center(id = "grapher") })  # Observer to center the mapper

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
      grapher::setNodeSize(id = "grapher", r = as.integer(vertex_radii))
    }
  })

  # Node color observer
  observeEvent(input$colorFunc, {
    #browser()
    f <- color_funcs[[input$colorFunc]]
    res <- f(M, X)
    binned_color <- bin_color(res, output_format = 'hex9')
    grapher::setNodeColor(id = "grapher", color = binned_color)
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
  proxy <- DT::dataTableProxy('mdatatable')
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
