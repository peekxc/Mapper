#' Ball Cover
#'
#' @docType class
#' @description This class provides a cover whose open sets are formed by the union of \deqn{\epsilon}-balls centered 
#' about each point. Using this class requires the \code{RANN} package to be installed, and thus explicitly assumes
#' the filter space endowed with the euclidean metric. 
#' 
#' @field epsilon := radius of the ball to form around each point
#' @author Matt Piekenbrock
#' @export
BallCover <- R6::R6Class(
  classname = "BallCover",
  inherit = CoverRef,
  public = list(epsilon=NA)
)

BallCover$set("public", "initialize", function(filter_values, ...){
  super$initialize(filter_values, typename="ball")
  params <- list(...)
  if ("epsilon" %in% names(params)){ self$epsilon <- params[["epsilon"]] }
})

BallCover$set("public", "construct_cover", function(){
  stopifnot(!is.na(self$epsilon))
  if (!requireNamespace("RANN", quietly = TRUE)){
    stop("Package \"RANN\" is needed for to use this cover.", call. = FALSE)
  }
  
  ## Construct the balls
  ball_cover  <- RANN::nn2(self$filter_values, query = self$filter_values, searchtype = "radius", radius = self$epsilon)
  
  ## Union them together 
  ds <- union_find(private$.filter_size)
  apply(ball_cover$nn.idx, 1, function(idx){
    connected_idx <- idx[idx != 0] - 1L
    if (length(connected_idx) > 0){
      ds$union_all(connected_idx)
    }
  })
  
  ## Construct the intersections between the open sets and the data
  cc <- ds$connected_components()
  self$index_set <- as.character(unique(cc))
  ls <- lapply(self$index_set, function(idx){
    which(cc == as.integer(idx))
  })
  self$level_sets <- structure(ls, names=self$index_set)
  
  ## Always return self 
  invisible(self)
})
