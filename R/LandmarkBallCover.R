#' Ball Cover
#'
#' @docType class
#' @description This class provides a cover whose open sets are formed by \deqn{\epsilon}-balls centered
#' about each point. Using this class requires the \code{RANN} package to be installed, and thus explicitly assumes
#' the filter space endowed with the euclidean metric.
#'
#' This differs from the DisjointBallCover in that it does NOT union intersecting cover sets.
#'
#' @field epsilon := radius of the ball to form around each point
#' @author Cory Brunsion, Yara Skaf
#' @export

library(proxy)

# TODO: make default spec with default seed 1l -> user can use default 1 or specify other seed
# Seed methods: SPEC, RAND, ECC (or specify seed_index, unspecified uses the first index)
LandmarkBallCover <- R6::R6Class(
  classname = "LandmarkBallCover",
  inherit = CoverRef,
  public = list(epsilon=NULL, num_sets=NULL, seed_index=1, seed_method="SPEC")
)

## initialize ------
LandmarkBallCover$set("public", "initialize", function(...){
  super$initialize(typename="ys_ball")
  params <- list(...)
  if ("epsilon" %in% names(params)){ self$epsilon <- params[["epsilon"]] }
  if ("num_sets" %in% names(params)){ self$num_sets <- params[["num_sets"]] }
  if ("seed_index" %in% names(params)){ self$seed_index <- params[["seed_index"]] }
  if ("seed_method" %in% names(params)){ self$seed_method <- params[["seed_method"]] }
})

## validate ------
LandmarkBallCover$set("public", "validate", function(filter){
  stopifnot(!is.null(self$epsilon) || !is.null(self$num_sets))
  writeLines("\n(eps, num_sets, seed_index, seed_method):\n")
  print(self$epsilon)
  print(self$num_sets)
  print(self$seed_index)
  print(self$seed_method)
})

## format ----
LandmarkBallCover$set("public", "format", function(...){
  titlecase <- function(x){
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }
  if (!is.null(self$num_sets)) {
    sprintf("%s Cover: (num_sets = %s, seed_index = %s)", titlecase(private$.typename), self$num_sets, self$seed_index)
  }
  if (!is.null(self$epsilon)) {
    sprintf("%s Cover: (epsilon = %.2f, seed_index = %s)", titlecase(private$.typename), self$epsilon, self$seed_index)
  }
})

## construct_cover ------
LandmarkBallCover$set("public", "construct_cover", function(filter, index=NULL){
  writeLines("Index:")
  print(index)
  print("constructing")
  if (!requireNamespace("RANN", quietly = TRUE)){
    stop("Package \"RANN\" is needed for to use this cover.", call. = FALSE)
  }
  print("validating")
  self$validate()

  print("filtering")
  ## Get filter values
  fv <- filter()
  f_dim <- ncol(fv)
  f_size <- nrow(fv)


  ## Set the seed index if necessary
  if(is.null(index)){
    if(all(self$seed_method == "RAND")) {
      self$seed_index = sample(1:f_size, 1)
      print(self$seed_index)
      print("rand")
    }
    if(all(self$seed_method == "ECC")) {
      #self$seed_index = sample(1:f_size, 1)
      print(self$seed_index)
      print("ecc")
    }

  }


  if (!is.null(self$num_sets)) {
    print("seed")
    print(self$seed_index)
    print(self$num_sets)
    eps_lm <- landmarks(x=fv, n=self$num_sets, seed_index=self$seed_index)
    print(eps_lm)
    ## Get distance from each point to landmarks
    dist_to_lm <- proxy::dist(fv, fv[eps_lm,,drop=FALSE])

    orderedIndices = t(apply(dist_to_lm,1, sort))
    max = which.max(orderedIndices[,1])
    self$epsilon = dist_to_lm[max]    # ball radius should be distance of the farthest point from the landmark set so that all points are in at least one ball
    pts_within_eps <- function(lm_dist){ which(lm_dist <= self$epsilon) }

    ## Construct the index set + the preimages
    self$index_set <- as.character(eps_lm)
    print(self$index_set)
    self$level_sets <- structure(as.list(apply(dist_to_lm, 2, pts_within_eps)), names=self$index_set)
    print("done level sets")
  }else{
    if (!is.null(self$epsilon)) {
      print("eps")

    }
  }




  ## Construct the balls
  # if spec, elif rand, elif ecc -> set self$seed_index
  # then if epsilon elif num_sets using self$seed_index
  # if (!is.null(self$seed_method)){
  #   if (all(self$seed_method == "RAND")) {
  #     self$seed_index = sample(1:f_size, 1)
  #     print(self$seed_index)
  #     print("rand")
  #     eps_lm <- landmarks(fv, self$num_sets, seed_index=self$seed_index)
  #
  #     ## Get distance from each point to landmarks
  #     dist_to_lm <- proxy::dist(fv, fv[eps_lm,,drop=FALSE])
  #     pts_within_eps <- function(lm_dist){ which(lm_dist < self$epsilon) } # TODO: does this make sense? we shouldnt be able to sepcify both eps and num_sets
  #
  #     ## Construct the index set + the preimages
  #     self$index_set <- as.character(eps_lm)
  #     writeLines("set index set")
  #     self$level_sets <- structure(as.list(apply(dist_to_lm, 2, pts_within_eps)), names=self$index_set)
  #     writeLines("set level sets")
  #   }
  #   if (all(self$seed_method == "ECC")) {
  #     print("ecc")
  #   }
  # }else{
  #   if (!is.null(self$num_sets)){
  #     print("seed")
  #     print(self$seed_index)
  #     print(self$num_sets)
  #     eps_lm <- landmarks(x=fv, n=self$num_sets, seed_index=self$seed_index)
  #     print(eps_lm)
  #     ## Get distance from each point to landmarks
  #     dist_to_lm <- proxy::dist(fv, fv[eps_lm,,drop=FALSE])
  #
  #     orderedIndices = t(apply(dist_to_lm,1, sort))
  #     max = which.max(orderedIndices[,1])
  #     self$epsilon = dist_to_lm[max]    # ball radius should be distance of the farthest point from the landmark set so that all points are in at least one ball
  #     pts_within_eps <- function(lm_dist){ which(lm_dist <= self$epsilon) }
  #
  #     ## Construct the index set + the preimages
  #     self$index_set <- as.character(eps_lm)
  #     print(self$index_set)
  #     self$level_sets <- structure(as.list(apply(dist_to_lm, 2, pts_within_eps)), names=self$index_set)
  #   }else {
  #     writeLines("else")
  #     eps_lm <- landmarks(fv, 2L)
  #
  #     ## Get distance from each point to landmarks
  #     dist_to_lm <- proxy::dist(fv, fv[eps_lm,,drop=FALSE])
  #     pts_within_eps <- function(lm_dist){ which(lm_dist < self$epsilon) }
  #
  #     ## Construct the index set + the preimages
  #     self$index_set <- as.character(eps_lm)
  #     self$level_sets <- structure(as.list(apply(dist_to_lm, 2, pts_within_eps)), names=self$index_set)
  #   }
  # }


  # print("making balls")
  # ## Construct the balls
  # ball_cover  <- RANN::nn2(fv, query = fv, searchtype = "radius", radius = self$epsilon)
  #
  # ## Assign the centers of the balls as the index set
  # self$index_set <- as.character(ball_cover[["nn.idx"]][,1])
  #
  # ## Calculate level (pullback) sets
  # ls <- lapply(seq_len(nrow(ball_cover[["nn.idx"]])), function(i) ball_cover[["nn.idx"]][i,])
  # ls <- lapply(ls, function(i) Filter(function(x) any(x != 0), i))
  #
  # self$level_sets <- structure(ls, names=self$index_set)
  print("test")
  print(index)
  if (!missing(index)){
    print("hello")
    return(self$level_sets[[index]]) }
  print("done")

  ## Always return self
  invisible(self)
})
