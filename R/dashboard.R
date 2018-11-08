#' Mapper dashboard
#' @description Starts an interactive dashboard for analyzing Mapper constructions using Shiny.
#' @param M MapperRef object.
#' @param X data.frame to display in an interactive table. See details. 
#' @param node_color_f List of color functions. See details.
#' @param dash_config Dashboard configuration. See details. 
#'
#' @details The data.frame supplied by the \code{X} argument is displayed in an interactive table in the dashboard.
#' Supplying the original data as a data.frame with named columns may be helpful in situations where the actual \code{X} mapper was constructed with required 
#' preprocessing (such as feature scaling), whereas the original, unprocessed data.frame is better suited for direct analysis.
#'
#' If supplied, \code{node_color_f} must be a named list of functions which accepts as input a MapperRef object as the \code{M} argument and returns as output a vector of numerical values, one per vertex.
#' These values are automatically converted to color hex codes based on the current palette. If \code{node_color_f} isn't supplied, a few default color functions are supplied.
#' 
#' @import shiny
#' @export
dashboard <- function(M, X, node_color_f = "default", dash_config=list(node_min=5L, node_max=10L)){
  library("shiny")
  
  ## Make sure shorthand 'M' is defined
  if (!"MapperRef" %in% class(M)){ stop("'dashboard' must take as input a MapperRef instance.") }
  G <- (M$as_grapher(construct_widget = TRUE) %>% grapher::enableForce())

  ## Get the UI components
  ui_file <- system.file(file.path("dashboard", "ui.R"), package = "Mapper")
  source(file = ui_file, local = TRUE)

  ## Node color functions
  if (missing(node_color_f) || node_color_f == "default"){
    if (is(X, "dist")){ stop("Default color functions cannot be applied when 'X' is given as a dist object.") }
    if (is.null(colnames(X))) { colnames(X) <- paste0("dim", 1:dim(X)[2]) }
    color_funcs <- new.env(parent = .BaseNamespaceEnv)
    color_file <- system.file(file.path("dashboard", "components", "default_color_functions.R"), package = "Mapper")
    source(color_file, local = color_funcs)
    sapply(colnames(X), function(dim_name){ color_funcs[[dim_name]] <- make_Dim_f(dim_name) })
  } else {
    color_funcs <- node_color_f
  }

  ## Get the server components
  server_file <- system.file(file.path("dashboard", "server.R"), package = "Mapper")
  source(file = server_file, local = TRUE)

  ## Return the shiny app
  # shiny::shinyApp(ui = ui, server = server)
  shiny::runApp(list(ui=ui, server=server))
}

## Auxillary function to make a closure for each dimension of X
make_Dim_f <- function(dim_name){
  function(M, X, ...){
    agg_dim_val <- sapply(M$vertices, function(n_idx){
      tmp <- as.matrix(X[n_idx, which(dim_name == colnames(X))])
      apply(tmp, 1, mean)
    })
    agg_node_val <- sapply(agg_dim_val, mean)
    return(agg_node_val)
  }
}
