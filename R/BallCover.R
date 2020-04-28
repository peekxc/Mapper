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
  public = list(epsilon=NULL)
)

## initialize ------
BallCover$set("public", "initialize", function(...){
  super$initialize(typename="ball")
  params <- list(...)
  if ("epsilon" %in% names(params)){ self$epsilon <- params[["epsilon"]] }
})

## validate ------
BallCover$set("public", "validate", function(filter){
  stopifnot(!is.null(self$epsilon))
})

## format ----
BallCover$set("public", "format", function(...){
  titlecase <- function(x){
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }
  sprintf("%s Cover: (epsilon = %.2f)", titlecase(private$.typename), self$epsilon)
})

## construct_cover ------
BallCover$set("public", "construct_cover", function(filter, index=NULL){
  if (!requireNamespace("RANN", quietly = TRUE)){
    stop("Package \"RANN\" is needed for to use this cover.", call. = FALSE)
  }
  self$validate()

  if(missing(index) || is.null(index)){
    ## Get filter values
    fv <- filter()
    f_dim <- ncol(fv)
    f_size <- nrow(fv)

    ## Construct the balls
    ball_cover  <- RANN::nn2(fv, query = fv, searchtype = "radius", radius = self$epsilon)

    ## Union them together
    ds <- union_find(f_size)
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
  }
  if (!missing(index)){ return(self$level_sets[[index]]) }

  ## Always return self
  invisible(self)
})
