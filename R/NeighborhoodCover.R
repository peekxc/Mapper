#' Neighborhood Cover
#'
#' @docType class
#' @description This class provides a cover whose open sets are formed by \code{k}-neighborhoods about a landmark
#' set. If no seed or seed_method is specified, default behavior uses the first data point as the seed. Cover sets
#' may contain more than \code{k} points if there are more than \code{k} points equidistant from the central point.
#' Using this class requires the \code{RANN} package to be installed, and thus explicitly assumes the filter space
#' endowed with the euclidean metric.
#'
#' @field k := desired number of neighbord to include in a cover set
#' @field seed_index := index of data point to use as the seed for landmark set calculation
#' @field seed_method := method to select a seed ("SPEC" : user specified index | "RAND" : random index
#'  | "ECC" : point with highest eccentricity in the filter space)
#' @author Yara Skaf, Cory Brunsion
#' @family cover
#' @export

library(proxy)

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
  sprintf("%s Cover: (k = %s, seed index = %s)", titlecase(private$.typename), self$k, self$seed_index)
})

## construct_cover ------
NeighborhoodCover$set("public", "construct_cover", function(filter, index=NULL, dist_method = "euclidean"){
  if (!requireNamespace("RANN", quietly = TRUE)){
    stop("Package \"RANN\" is needed for to use this cover.", call. = FALSE)
  }
  self$validate(filter)
  stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))

  ## Get filter values
  fv <- filter()
  f_dim <- ncol(fv)
  f_size <- nrow(fv)

  if(is.null(index)){
    ## Set the seed index if necessary
    if(all(self$seed_method == "RAND")) { self$seed_index = sample(1:f_size, 1) }
    if(all(self$seed_method == "ECC")) {  self$seed_index = which.max(eccentricity(from=fv, x=fv)) }

    ## Compute the set of k-nhds
    C = c() # create an empty list to store indices of centers
    k_nhds = list() # create empty list of lists to store points in each neighborhood
    ptsLeft = c(1:f_size) # keep track of indices that are still available to be chosen as a center

    nextC = self$seed_index # use the seed as the first center
    while(TRUE){
      # add the new center to the list and compute its k-neighborhood
      C = append(C, nextC)
      neighbors = proxy::dist(matrix(fv[nextC,], ncol=f_dim), fv, method = dist_method)
      sortedNeighbors = sort(neighbors)

      # include all points that are equidistant from the center, even if those points create a nhd size > k
      i = 1
      while( ((self$k+i) <= f_size) && (sortedNeighbors[self$k] == sortedNeighbors[self$k + i]) ){i = i + 1}
      nhd = order(neighbors)[1:(self$k + i - 1)]
      k_nhds = append(k_nhds, list(nhd))

      # points that are included in a nhd are no longer eligible to be chosen as centers
      ptsLeft = setdiff(ptsLeft, nhd)
      if(length(ptsLeft) <= 0){break}

      # select the next center
      dists = proxy::dist(matrix(fv[C,], ncol=f_dim), matrix(fv[ptsLeft,], ncol=f_dim), method = dist_method)
      sortedDists = matrix(apply(dists,2,sort),ncol=length(ptsLeft))
      max = which.max(sortedDists[1,])
      nextC = ptsLeft[max]
    }

    ## Construct the index set and level sets
    self$index_set = as.character(C)
    self$level_sets = structure(k_nhds, names=self$index_set)
  }
  if (!missing(index)){ return(self$level_sets[[index]]) }

  ## Always return self
  invisible(self)
})
