## shiny_modules.R
## Author: Matt Piekenbrock

#' Data input module UI function
#' @description Creates the UI to select a data set for use with \code{ShinyMapper}.
#' @param id character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#' @param datasets optional named list of data sets. See details. 
#' @details if \code{datasets} is supplied, it must be a named list, where each name corresponds
#' to the data set name and each element is either a matrix or a function that accepts as input 
#' a set of indices and returns either a data matrix or a dist object. See \code{\link{MapperRef}} for details. 
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements.
#' @export
MapperDataInput <- function(id, datasets=NULL){
  ns <- NS(id)
  data_supplied <- !(missing(datasets) || is.null(datasets))
  if (!data_supplied){
    extdata <- file.path(system.file(package="Mapper"), "extdata")
    datasets <- list(
      "circle" = function(){ readRDS(file.path(extdata, "noisy_circle.rds")) }, 
      "WVS" = function() { readRDS(file.path(extdata, "wvs_us_wave6_ex.rds")) }
    )
    selectInput(ns("data_name"), label = "Data", choices = datasets, selected = head(names(datasets), 1))
  } else {
    functor <- function(symbol){ function(){ eval(parse(text=symbol), envir = .GlobalEnv) } }
    datasets <- structure(lapply(ls(.GlobalEnv), functor), names = ls(.GlobalEnv))
    selectizeInput(ns("data_name"), label = "Data", choices = datasets, multiple = FALSE, selectize = TRUE)
  }
}

#' Shiny data 
#' @param input, output, session standard \code{shiny} boilerplate
#' @return reactive expression that returns a function containing the data set.
MapperData <- function(input, output, session){
  return(reactive({ input$data_name }))
}

#' Mapper filter shiny UI module
#' @description Creates the UI to select a data set for use with \code{ShinyMapper}.
#' @param id character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#' @param datasets optional named list of data sets. See details. 
#' @details if \code{datasets} is supplied, it must be a named list, where each name corresponds
#' to the data set name and each element is either a matrix or a function that accepts as input 
#' a set of indices and returns either a data matrix or a dist object. See \code{\link{MapperRef}} for details. 
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements.
#' @export
MapperFilterInput <- function(){
  selectInput(ns("filter_name"), label = "Filter", choices = c("PCA", "ICA", "Eccentricity", "MDS", "LE", "tSNE", "UMAP", "Custom"))
}

#' Shiny Mapper module UI function
#' @description Creates UI elements that parameterize the Mapper construction for Shiny apps.
#' @param id character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements
#' @import shiny
ShinyMapperInput <- function(id){
  # stopifnot(is.list(X), !is.null(names(X)))
  ns <- NS(id)
  tagList(
    # selectInput(ns("data_name"), label = "Input data", choices = names(X), selected = head(names(X), 1)),
    actionButton(ns("construct_mapper"), label = "Construct Mapper", width = "100%"),
    selectInput(ns("filter_name"), label = "Filter", choices = c("PCA", "ICA", "Eccentricity", "MDS", "LE", "tSNE", "UMAP", "Custom")),
    uiOutput(ns("filter_options")),
    selectInput(ns("cover_type"), label="Cover", 
                choices = c("fixed interval", "restrained interval", "adaptive", "ball"), 
                selected = "fixed interval"),
    uiOutput(ns("cover_options")),
    selectInput(ns("metric"), label = "Metric", choices=proxy::pr_DB$get_entry_names(), 
                selected = "Euclidean"),
    selectInput(ns("clustering"), label="Clustering algorithm (or linkage criterion)", 
                choices=c("ward.D", "ward.D2", "single", "complete", 
                          "average", "mcquitty", "median", "centroid", "custom"), 
                selected = "single"),
    selectInput(ns("cutoff_heuristic"), label = "Cutoff Heuristic", 
                choices = c("histogram", "continuous"), 
                selected="continuous")
  )
}
# M$use_cover(filter_values = .GlobalEnv[[input$filter_name]], typename = input$cover_type)
ShinyCoverInput <- function(id){
  ns <- NS(id)
  tagList(
    conditionalPanel(condition = "(input.cover_type == 'fixed interval') || (input.cover_type == 'restrained interval')", 
      tagList(
        sliderInput(ns("number_intervals"), label = "Number intervals", min = 1, max = 50, value = 10, step = 1),
        sliderInput(ns("percent_overlap"), label = "Percent overlap", min = 0, max = 100, value = 50, step = 1) 
      )
    ), 
    conditionalPanel(condition = "input.cover_type == 'adaptive interval'", 
      tagList(
        sliderInput(ns("number_intervals"), label = "Number intervals", min = 1, max = 50, value = 10, step = 1),
        sliderInput(ns("percent_overlap"), label = "Percent overlap", min = 0, max = 100, value = 50, step = 1), 
        numericInput(ns("quantile_method"), label = "Quantile method", min = 1, max = 9, step = 1L, value = 7) 
      )
    ), 
    conditionalPanel(condition = "input.cover_type == 'ball'", 
      tagList(
       sliderInput(ns("epsilon"), label = "Radius (normalized)", min = 0, max = 1, step = 0.001)
      )
    )
  )
  
}


#' Shiny Mapper module server function
#' @description module server function for shiny mapper.
#' @param input, output, session standard \code{shiny} boilerplate
#' @param dataset reactive expression that returns the data set. 
#' @param cover reactive expression that returns the cover. 
#' @param M optional MapperRef instance (non-reactive). See details. 
#' @details If \code{M} is supplied, the parameters of the current Mapper are copied to the inputs, and that instance 
#' is then updated using the controls. Otherwise, a new \code{\link{MapperRef}} instance is constructed and returned from 
#' this method, parameterized by the associated inputs.
#' @return rective expression that constructs the Mapper.
ShinyMapper <- function(input, output, session, dataset, cover, M = NULL){
  mapper_supplied <- !missing(M)
  if (mapper_supplied){ stopifnot(is(M, "MapperRef")) }
  M <- if (mapper_supplied) { M } else { new(MapperRef) }
  MapperObj <- reactive({
    M <- MapperRef$new(dataset)
    M$cover <- cover()
    M$use_distance_measure(input$metric)
    M$use_clustering_algorithm(cl = input$clustering, cutoff_method = input$cutoff_heuristic)
    return(M$construct_k_skeleton(k=1L))
  })
  return(MapperObj)
}
