#' Mapper dashboard
#' @description Starts an interactive dashboard for analyzing Mapper constructions using Shiny.
#' @param mapper_obj MapperRef object.
#' @param X data.frame to display in lower table and subset on node selection.
#' @param node_color_f List of color functions. See details.
#'
#' @details The data.frame supplied by the \code{X} argument is displayed in an interactive table in the dashboard.
#' By default, if \code{X} is not supplied it is assigned the original data the Mapper object was constructed with.
#' Supplying a data.frame with named columns may be helpful in situations where the point cloud data the mapper was constructed with required 
#' preprocessing (such as feature scaling), whereas the original, unprocessed data.frame is better suited for direct analysis.
#'
#' If supplied, \code{node_color_f} must be a named list of functions which accepts as input a MapperRef object and returns as output a vector of numerical values, one per node.
#' These values are automatically converted to color hex codes based on the current palette. 
#' 
#' If this isn't supplied, a few default color functions are supplied.
#'
#' @import shiny htmlwidgets
#' @importFrom DT DTOutput renderDT
#' @export
dashboard <- function(mapper_obj, X = mapper_obj$X, node_color_f = "default"){

  ## Make sure shorthand 'M' is defined
  if (!"MapperRef" %in% class(mapper_obj)){ stop("'dashboard' must take as input a mapper reference object.") }
  # M <- mapper_obj
  # G <- Mapper::grapher(M)
  # if (is(X, "dist")){ X <- data.frame(index=1:attr(M$X, "Size")) }
  # else { X <- as.data.frame(X) }
  # 
  # ## Get the UI components
  # ui_file <- system.file(file.path("dashboard", "ui.R"), package = "Mapper")
  # source(file = ui_file, local = TRUE)
  # 
  # ## Node color functions
  # if (missing(node_color_f) || node_color_f == "default"){
  #   if (is(X, "dist")){ stop("Default color functions cannot be applied when 'X' is given as a dist object.") }
  #   if (is.null(colnames(X))) { colnames(X) <- paste0("dim", 1:dim(X)[2]) }
  #   color_funcs <- new.env(parent = .BaseNamespaceEnv)
  #   color_file <- system.file(file.path("dashboard", "components", "default_color_functions.R"), package = "Mapper")
  #   source(color_file, local = color_funcs)
  #   sapply(colnames(X), function(dim_name){ color_funcs[[dim_name]] <- make_Dim_f(dim_name) })
  # } else {
  #   color_funcs <- node_color_f
  # }
  # 
  # ## Get the server components
  # server_file <- system.file(file.path("dashboard", "server.R"), package = "Mapper")
  # source(file = server_file, local = TRUE)
  # 
  # ## Return the shiny app
  # shiny::shinyApp(ui = ui, server = server)
}

## Auxillary function to make a closure for each dimension of X
make_Dim_f <- function(dim_name){
  function(M, X, ...){
    agg_dim_val <- sapply(M$G$nodes, function(n_idx){
      tmp <- as.matrix(X[n_idx, which(dim_name == colnames(X))])
      apply(tmp, 1, mean)
    })
    agg_node_val <- sapply(agg_dim_val, mean)
    return(agg_node_val)
  }
}
