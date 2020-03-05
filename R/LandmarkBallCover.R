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
LandmarkBallCover <- R6::R6Class(
  classname = "LandmarkBallCover",
  inherit = CoverRef,
  public = list(epsilon=NULL, num_sets=NULL, seed_index=1, seed_method="SPEC")
)

## initialize ------
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
  stopifnot(!is.null(self$epsilon) || !is.null(self$num_sets))
})

## format ----
LandmarkBallCover$set("public", "format", function(...){
  titlecase <- function(x){
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }
  if (!is.null(self$num_sets)) {
    sprintf("%s Cover: (num_sets = %s, seed_index = %s)", titlecase(private$.typename), self$num_sets, self$seed_index)
  }else{
    if (!is.null(self$epsilon)) {
      sprintf("%s Cover: (epsilon = %.2f, seed_index = %s)", titlecase(private$.typename), self$epsilon, self$seed_index)
    }
  }
})

## construct_cover ------
LandmarkBallCover$set("public", "construct_cover", function(filter, index=NULL){
  if (!requireNamespace("RANN", quietly = TRUE)){
    stop("Package \"RANN\" is needed for to use this cover.", call. = FALSE)
  }
  self$validate()

  ## Get filter values
  fv <- filter()
  f_size <- nrow(fv)

  ## Set the seed index if necessary
  if(is.null(index)){
    if(all(self$seed_method == "RAND")) {
      self$seed_index = sample(1:f_size, 1)
    }
    if(all(self$seed_method == "ECC")) {
      # todo
    }
  }

  if (!is.null(self$num_sets)) {
    eps_lm <- landmarks(x=fv, n=self$num_sets, seed_index=self$seed_index) # compute landmark set

    ## Get distance from each point to landmarks
    dist_to_lm <- proxy::dist(fv, fv[eps_lm,,drop=FALSE])
    pts_within_eps <- function(lm_dist){ which(lm_dist <= self$epsilon) }

    ## Construct the index set
    self$index_set <- as.character(eps_lm)

    ## Construct the preimages
    if(length(eps_lm) == 1){
      orderedIndices = apply(dist_to_lm,1, sort)
      max = which.max(orderedIndices)
      self$epsilon = dist_to_lm[max]    # ball radius should be distance of the farthest point from the landmark set so that all points are in at least one ball
      self$level_sets <- structure(as.list(list(t(apply(dist_to_lm, 2, pts_within_eps))[1,])), names=self$index_set)
    }else{
      orderedIndices = t(apply(dist_to_lm,1, sort))
      max = which.max(orderedIndices[,1])
      self$epsilon = dist_to_lm[max]    # ball radius should be distance of the farthest point from the landmark set so that all points are in at least one ball
      self$level_sets <- structure(as.list(apply(dist_to_lm, 2, pts_within_eps)), names=self$index_set)
    }
  }else if (!is.null(self$epsilon)) {
      eps_lm <- landmarks(x=fv, eps=self$epsilon, seed_index=self$seed_index) # compute landmark set
      dist_to_lm <- proxy::dist(fv, fv[eps_lm,,drop=FALSE]) # Get distance from each point to landmarks

      pts_within_eps <- function(lm_dist){ which(lm_dist <= self$epsilon) }

      ## Construct the index set
      self$index_set <- as.character(eps_lm)

      ## Construct the preimages
      if(length(eps_lm) == 1){
        self$level_sets <- structure(as.list(list(t(apply(dist_to_lm, 2, pts_within_eps))[1,])), names=self$index_set)
      }else{
        self$level_sets <- structure(as.list(apply(dist_to_lm, 2, pts_within_eps)), names=self$index_set)
      }
    }

  if (!missing(index)){
    return(self$level_sets[[index]]) }

  ## Always return self
  invisible(self)
})
