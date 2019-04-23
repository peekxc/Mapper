#' Mapper dashboard
#' @description Starts an interactive dashboard for analyzing Mapper constructions using Shiny.
#' @param M MapperRef object.
#' @param X data.frame to display in an interactive table. See details. 
#' @param node_color_f Named list of color functions. See details.
#' @param dash_config Dashboard configuration. See details. 
#'
#' @details 
#' \strong{WARNING: THIS APP IS EXPERIMENTAL} \cr
#' \cr  
#' The data.frame supplied by the \code{X} argument is displayed in an interactive table in the dashboard.
#' Supplying the original data as a data.frame with named columns may be helpful in situations where the actual \code{X} mapper was constructed with required 
#' preprocessing (such as feature scaling), whereas the original, unprocessed data.frame is better suited for direct analysis.
#'
#' If supplied, \code{node_color_f} must be a named list of functions which accepts as input a MapperRef object as the \code{M} argument and returns as output a vector of numerical values, one per vertex.
#' These values are automatically converted to color hex codes based on the current palette. If \code{node_color_f} isn't supplied, a few default color functions are supplied.
dashboard <- function(M, X, node_color_f = "default", dash_config=list(node_min=5L, node_max=10L)){
  requireNamespace("shiny", quietly = TRUE)
  requireNamespace("grapher", quietly = TRUE)
  
  ## Make sure shorthand 'M' is defined
  # stopifnot(nrow(M$X) == nrow(X))
  if (!"MapperRef" %in% class(M)){ stop("'dashboard' must take as input a MapperRef instance.") }

  ## Node color functions
  if (missing(node_color_f) || node_color_f == "default"){
    if (is(X, "dist")){ stop("Default color functions cannot be applied when 'X' is given as a dist object.") }
    if (is.null(colnames(X))) { colnames(X) <- paste0("dim", 1:dim(X)[2]) }
    color_funcs <- make_default_color_f(X)
  } else {
    color_funcs <- as.environment(node_color_f)
  }
  
  ## Get the app directory
  appDir <- system.file("dashboard_module", package = "Mapper")
  if (appDir == "") { stop("Could not find dashboard directory. Try re-installing `Mapper`.", call. = FALSE) }
  
  .dash_env[["test"]] <- 1L
  
  ## Make sure these are available to the app
  # mapply(assign, x=c("M", "X", "color_funcs", "dash_config"), value=list(M, X, color_funcs, dash_config), envir = .GlobalEnv)
  .GlobalEnv$M <- M
  .GlobalEnv$X <- X
  .GlobalEnv$color_funcs <- color_funcs
  .GlobalEnv$dash_config <- dash_config
  on.exit(rm(list=c("M", "X", "color_funcs", "dash_config"), envir=.GlobalEnv))
  
  ## Run the shiny app
  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}

## Generates a set of default functions for coloring the mapper
make_default_color_f <- function(X){
  ## The current environment passed in
  c_env <- new.env(parent = emptyenv())
  
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
  
  ## Auxillary function to make a closure for each averaging each dimension of a given data.frame X
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
  
  ## Make one function per column
  for (dim_name in colnames(X)){
    c_env[[dim_name]] <- make_Dim_f(dim_name)
  }
  return(c_env)
}
