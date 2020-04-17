#' Landmark Ball Cover
#'
#' @docType class
#' @description This class provides a cover whose open sets are formed by balls centered about each
#' point in a landmark set. Given a radius \deqn{\epsilon}, choose a set of landmark points via the
#' algorithm presented in Dłotko to produce a cover by balls of radius \deqn{\epsilon}. Alternatively,
#' given a number of cover sets \code{n}, choose \code{n} landmarks via maxmin algorithm. If no seed
#' or seed_method is specified, default behavior uses the first data point as the seed.
#'
#' This differs from BallCover.R in that it does NOT union intersecting cover sets.
#'
#' Using this class requires the \code{RANN} package to be installed, and thus explicitly assumes
#' the filter space endowed with the euclidean metric.
#'
#' @field epsilon := radius of the ball to form around each landmark point
#' @field num_sets := desired number of balls/cover sets
#' @field seed_index := index of data point to use as the seed for landmark set calculation
#' @field seed_method := method to select a seed ("SPEC" : user specified index | "RAND" : random index
#'  | "ECC" : point with highest eccentricity in the filter space)
#' @author Yara Skaf, Cory Brunsion
#' @family cover
#' @references Dłotko, Paweł. "Ball Mapper: A Shape Summary for Topological Data Analysis." (2019). Web.

library(proxy)

#' @export
LandmarkBallCover <- R6::R6Class(
  classname = "LandmarkBallCover",
  inherit = CoverRef,
  public = list(epsilon=NULL, num_sets=NULL, seed_index=1, seed_method="SPEC")
)

## initialize ------
#' @export
LandmarkBallCover$set("public", "initialize", function(...){
  super$initialize(typename="landmark_ball")
  params <- list(...)
  if ("epsilon" %in% names(params)){ self$epsilon <- params[["epsilon"]] }
  if ("num_sets" %in% names(params)){ self$num_sets <- params[["num_sets"]] }
  if ("seed_index" %in% names(params)){ self$seed_index <- params[["seed_index"]] }
  if ("seed_method" %in% names(params)){ self$seed_method <- params[["seed_method"]] }
})

## validate ------
LandmarkBallCover$set("public", "validate", function(filter){
  ## Get filter values
  fv <- filter()
  f_size <- nrow(fv)

  ## validate parameters
  stopifnot(!is.null(self$epsilon) || !is.null(self$num_sets)) # require either radius or # balls
  stopifnot(is.null(self$num_sets) || (self$num_sets <= f_size && self$num_sets > 0)) # cannot have more cover sets than data points
  stopifnot(is.null(self$epsilon) || self$epsilon >= 0) # radius must be positive
  stopifnot(self$seed_index <= f_size && self$seed_index > 0) # seed index must be within the range of the data indices
  stopifnot(all(self$seed_method == "RAND") || all(self$seed_method == "SPEC") || all(self$seed_method == "ECC")) # must use one of available seed methods
})

## format ----
LandmarkBallCover$set("public", "format", function(...){
  titlecase <- function(x){
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }
  if (!is.null(self$num_sets)) {
    sprintf("%s Cover: (number of sets = %s, seed index = %s)", titlecase(private$.typename), self$num_sets, self$seed_index)
  }else if (!is.null(self$epsilon)){
    sprintf("%s Cover: (epsilon = %.2f, seed index = %s)", titlecase(private$.typename), self$epsilon, self$seed_index)
  }
})

## construct_cover ------
LandmarkBallCover$set("public", "construct_cover", function(filter, index=NULL){
  if (!requireNamespace("RANN", quietly = TRUE)){
    stop("Package \"RANN\" is needed for to use this cover.", call. = FALSE)
  }
  self$validate(filter)

  ## Get filter values
  fv <- filter()
  f_size <- nrow(fv)

  ## Construct the balls
  if(is.null(index)){
    ## Set the seed index if necessary
    if(all(self$seed_method == "RAND")) { self$seed_index = sample(1:f_size, 1) }
    if(all(self$seed_method == "ECC")) {  self$seed_index = which.max(eccentricity(from=fv, x=fv)) }

    ## Compute the landmark set
    if (!is.null(self$num_sets)) { eps_lm <- unique(landmarks(x=fv, n=self$num_sets, seed_index=self$seed_index))
    } else if (!is.null(self$epsilon)) { eps_lm <- landmarks(x=fv, eps=self$epsilon, seed_index=self$seed_index) }

    ## Construct the index set
    self$index_set <- as.character(eps_lm)

    ## Get distance from each point to landmarks
    dist_to_lm <- proxy::dist(fv, fv[eps_lm,,drop=FALSE])
    pts_within_eps <- function(lm_dist){ which(lm_dist <= self$epsilon) }

    ## Calculate an epsilon if one was not given
    if (!is.null(self$num_sets)) {
      sortedDists = matrix(apply(dist_to_lm,1,sort),nrow=f_size,byrow=TRUE)
      max = which.max(sortedDists[,1])
      self$epsilon = sortedDists[max,1] # radius should be distance of the farthest pt from the landmark set so that all pts are in at least one ball
    }

    x = apply(dist_to_lm, 2, pts_within_eps)

    # if all level sets contain the same number of points, apply returns a matrix -> need to split columns into list elements
    if(is.matrix(x)){ self$level_sets <- structure(split(x, rep(1:ncol(x), each = nrow(x))), names=self$index_set)
    } else { self$level_sets <- structure(as.list(x), names=self$index_set) }

    print(self$level_sets)
  }
  if (!missing(index)){ return(self$level_sets[[index]]) }

  ## Always return self
  invisible(self)
})
