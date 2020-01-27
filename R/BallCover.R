#' Ball Cover
#' @docType class
#' @description This class provides a cover whose open sets are formed by landmarks that 
#' are each at least \code{epsilon}-distance from each other. The landmarks are chosen 
#' by the \code{\link[landmarks]{maxmin}} method (1), in a similar way as e.g. witness 
#' complexes are constructed (2). One interpretation of this cover in the context 
#' of mapper is given by Dłotko (3).
#' 
#' @field epsilon radius of the ball to form the \eqn{\epsilon}-cover
#' @field method metric to create the \eqn{\epsilon}-cover with  
#' @references 
#' 1. Arafat, Naheed Anjum, Debabrota Basu, and Stéphane Bressan. "Topological Data Analysis with \eqn{\epsilon}-net Induced Lazy Witness Complex." arXiv preprint arXiv:1906.06122 (2019).
#' 
#' 2. De Silva, Vin, and Gunnar E. Carlsson. "Topological estimation using witness complexes." SPBG 4 (2004): 157-166.
#' 
#' 3. Dłotko, Paweł. "Ball mapper: a shape summary for topological data analysis." arXiv preprint arXiv:1901.07410 (2019).
#' @author Matt Piekenbrock
#' @family cover
#' @export
BallCover <- R6::R6Class(
  classname = "BallCover",
  inherit = CoverRef,
  public = list(epsilon=NULL, method="euclidean")
)

## initialize ------
BallCover$set("public", "initialize", function(...){
  super$initialize(typename="ball")
  params <- list(...)
  if ("epsilon" %in% names(params)){ self$epsilon <- params[["epsilon"]] }
  if ("method" %in% names(params)){ self$method <- params[["method"]] }
})

## validate ------
BallCover$set("public", "validate", function(filter){
  stopifnot(!is.null(self$epsilon))
  stopifnot(is.character(self$method) || is.function(self$method))
})

## format ----
BallCover$set("public", "format", function(...){
  titlecase <- function(x){
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }
  sprintf("%s Cover: (epsilon = %.2f)", titlecase(private$.typename), self$epsilon)
})

## construct ------
BallCover$set("public", "construct", function(filter, index=NULL, cache=TRUE){
  self$validate()
  
  ## Check to see if need to reconstruct cover or not. 
  if (private$.cached && !is.null(self$sets)){ 
    if (!missing(index)){ return(self$sets[index]) }
    return(self$sets) 
  }
  
  ## Get filter values 
  fv <- filter()
  
  ## Construct the balls
  eps_lm <- landmarks(fv, eps=self$epsilon, dist_method=self$method)
  
  ## Get distance from each point to landmarks 
  dist_to_lm <- proxy::dist(fv, fv[eps_lm,,drop=FALSE], method = self$method)
  pts_within_eps <- function(lm_dist){ which(lm_dist < self$epsilon) }
  
  ## Construct the index set + the preimages
  index_set <- as.character(eps_lm)
  sets <- structure(as.list(apply(dist_to_lm, 2, pts_within_eps)), names=index_set)
  
  ## Very dumb way to support index queries
  if (!missing(index)){
    stopifnot(index %in% index_set)
    sets <- sets[index]
  }
  
  ## Either cache and return, or just return 
  if (cache){
    private$.cached <- TRUE
    self$sets <- modifyList(self$sets, sets)
    return(invisible(self$sets))
  }
  return(sets)
})

## Returns the pairs of sets to check for potential non-empty intersection
BallCover$set("public", "neighborhood", function(filter, k){
  if (length(self$sets) > 0){
    fv <- filter()
    eps_lm <- as.integer(names(self$sets)) ## the names correspond to the landmark indices
    dist_to_lm <- as.matrix(proxy::dist(fv[eps_lm,,drop=FALSE], method = self$method)) ## distance between landmarks
    within_eps <- function(lm_dist){ unname(which((lm_dist/2.0) < self$epsilon)) }
    pairs <- do.call(rbind, lapply(1:length(eps_lm), function(i){
      cbind(i, within_eps(dist_to_lm[i,,drop=FALSE]))
    }))
    return(apply(pairs, 2, function(idx){ self$index_set[idx] }))
  } else {
    return(super$neighborhood(filter, k))
  }
})

BallCover$set("public", "plot", function(filter, add=FALSE, index=NULL, ball_opt=NULL, ...){

  ## Get filter values 
  fv <- filter()
  f_dim <- ncol(fv)
  f_size <- nrow(fv)
  
  stopifnot(self$index)
  as.integer(m$cover$index_set)
  if (missing(add) || !add){
    set_rng <- range(fv)
    graphics::plot.new()
    graphics::plot.window(xlim=set_rng[1], ylim=c(0,n_overlap+1))
    
    xyasp <- par("pin")
    xycr <- diff(par("usr"))[c(1, 3)]
    ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
    
    theta <- 2 * pi/100
    angles <- seq(0, 2 * pi - theta, by = theta)
    
    
    for (id in self$index_set) {
      xy <- as.vector(fv[as.integer(id),])
      xv <- cos(angles) * self$eps + xy[1]
      yv <- sin(angles) * self$eps * ymult + xy[2]
      polygon(xv, yv, border = border, col = col[circle], lty = lty, 
              density = density, angle = angle, lwd = lwd)
    }
  }

  
  
})
