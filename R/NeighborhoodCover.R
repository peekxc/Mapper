#' Ball Cover
#'
#' @docType class
#' @description This class provides a cover whose open sets are formed by \deqn{\epsilon}-balls centered
#' about each point. Using this class requires the \code{RANN} package to be installed, and thus explicitly assumes
#' the filter space endowed with the euclidean metric.
#'
#' This differs from the BallCover in that it does NOT union intersecting cover sets.
#'
#' @field epsilon := radius of the ball to form around each point
#' @author Cory Brunsion, Yara Skaf
#' @export

library(proxy)

# Seed methods: SPEC (specify index), RAND (random index), ECC (seed with highest eccentricity data point)
# Default: specified index using first data point (seed_method = "SPEC", seed_index = 1)
#' @export
NeighborhoodCover <- R6::R6Class(
  classname = "NeighborhoodCover",
  inherit = CoverRef,
  public = list(k=NULL, seed_index=1, seed_method="SPEC")
)

## initialize ------
#' @export
NeighborhoodCover$set("public", "initialize", function(...){
  super$initialize(typename="neighborhood")
  params <- list(...)
  if ("k" %in% names(params)){ self$k <- params[["k"]] }
  if ("seed_index" %in% names(params)){ self$seed_index <- params[["seed_index"]] }
  if ("seed_method" %in% names(params)){ self$seed_method <- params[["seed_method"]] }
})

## validate ------
NeighborhoodCover$set("public", "validate", function(filter){
  ## Get filter values
  fv <- filter()
  f_size <- nrow(fv)

  ## validate parameters
  stopifnot(!is.null(self$k)) # require nieghborhood size
  stopifnot(self$k >= 2 && self$k <= f_size) # size must be at least 2 points
  stopifnot(self$seed_index <= f_size && self$seed_index > 0) # seed index must be within the range of the data indices
  stopifnot(all(self$seed_method == "RAND") || all(self$seed_method == "SPEC") || all(self$seed_method == "ECC")) # must use one of available seed methods
})

## format ----
NeighborhoodCover$set("public", "format", function(...){
  titlecase <- function(x){
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }
  sprintf("%s Cover: (k = %s, seed_index = %s)", titlecase(private$.typename), self$k, self$seed_index)
})

## construct_cover ------
NeighborhoodCover$set("public", "construct_cover", function(filter, index=NULL){
  if (!requireNamespace("RANN", quietly = TRUE)){
    stop("Package \"RANN\" is needed for to use this cover.", call. = FALSE)
  }
  self$validate(filter)

  ## Get filter values
  fv <- filter()
  f_size <- nrow(fv)

  ## Set the seed index if necessary
  if(is.null(index)){
    if(all(self$seed_method == "RAND")) {
      self$seed_index = sample(1:f_size, 1)
    }
    if(all(self$seed_method == "ECC")) {
      ecc = eccentricity(from=fv, x=fv)
      self$seed_index = which.max(ecc)
    }
  }

  eps_lm <- landmarks(x=fv, k=self$k, seed_index=self$seed_index) # compute landmark set

  ## Construct the index set
  self$index_set <- as.character(attr(eps_lm,"names"))
  self$level_sets <- eps_lm

  if (!missing(index)){ return(self$level_sets[[index]]) }

  ## Always return self
  invisible(self)
})
