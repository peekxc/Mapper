library("shiny")

print(ls(Mapper:::.dash_env))

## Ensure required variables are in the global environment
if (!exists("M", envir = .GlobalEnv)){ stop("Error creating server; there was no 'MapperRef' object associated with the name 'M'.") }
if (!exists("X", envir = .GlobalEnv)){ stop("Error creating server; there was no 'data.frame' object associated with the name 'X'.") }
if (!exists("color_funcs", envir = .GlobalEnv)){ stop("Error creating server; there was no 'color_funcs' object.") }
if (!exists("dash_config", envir = .GlobalEnv)){ stop("Error creating server; there was no 'dash_config' object.") }

## Make the UI
ui <- navbarPage(title = "Mapper dashboard",
  tabPanel("Visualization", 
   sidebarLayout(
     sidebarPanel(
       nodeSizingInput("node_sizing"),
       nodeColorInput("node_color", color_f = as.list(color_funcs)),
       textOutput("summary")
     ), 
     mainPanel(
       fluidRow(column(12, mapperVisOutput("mapper_vis"), style="height: 60vh !important;")), 
       fluidRow(
         column(6, brushScatterOutput("scatter_out")),
         column(6, linkedTableOutput("data_table"))
       )
     )
   )
  ), 
  navbarMenu("Force Options") 
    # actionButton("force_button", label = "Enable force")
  #)
)

## Make the server
server <- function(input, output, session) {

  ## Variables
  cnames <- colnames(X)
  
  ## Choose output type to handling interaction with visualization
  mapper_vis <- callModule(mapperVis, id="mapper_vis", mapper=reactive(M), X=reactive(X), output_type=reactive("grapher"))
  
  ## Configure interactive observers specific to mapper
  node_sizing_f <- reactive(function(mapper){ sapply(mapper$vertices, length) })
  node_sizes <- callModule(nodeSizing, id="node_sizing", mapper_vis=mapper_vis, v_f=node_sizing_f)
  node_colors <- callModule(nodeColor, id="node_color", mapper_vis=mapper_vis, color_f=reactive(as.list(color_funcs)))
  
  ## Compose generic reactive input modules/sources for interactivity
  callModule(brushScatter, id="scatter_out", mapper_vis=mapper_vis, X = reactive(X), xvar = reactive(cnames[[1]]), yvar = reactive(cnames[[2]]), reset = node_colors)

  ##
  nodeHighlight(mapper_vis, reset=node_colors)
  
  browser()
  # output$summary <- renderText({
  #   sprintf("%d observation(s) selected", sum(selected_data()[["selected_"]]))
  # })

  ## Register observers
  # observeEvent(selected_data, { ## Enact the highlighting
  #   selected <- selected_data()[["selected_"]]
  #   selected_idx <- which(selected)
  #   v_col <- switch(as.character(all(selected == FALSE)), 
  #     "TRUE" = node_colors(), 
  #     "FALSE" = {
  #       v_weight <- sapply(mapper_vis$mapper$vertices, function(v) sum(v %in% selected_idx))
  #       bin_color(v_weight, col_pal = viridis::inferno(100, begin = 0.15, end = 0.8, alpha = seq(0.65, 1, length.out = 100)))
  #   })
  #   print("highlight color: ")
  #   print(unique(v_col))
  #   if (mapper_vis$output_type == "grapher"){
  #     grapher::setNodeColor(id = mapper_vis$output_id, color = v_col) 
  #   }
  # })

  ## Outputs
  callModule(linkedTable, id="data_table", dt=X) ## construct the data table

  session$onSessionEnded(stopApp)
}

shinyApp(ui=ui, server=server)

