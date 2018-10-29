#' Adaptive Cover
#'
#' @docType class
#' @description This class provides cover which 'adapts' the lengths of the open sets to the density of 
#' the filter values. The cutoff points are chosen with pre-selected quantile method, see \code{?quantile}
#' for more details on the calculation. Note that the percentage that each bin overlaps is calculated before 
#' transforming the cutoff points, so the interpretation of the bins as overlapping by \code{percent_overlap}%
#' isn't quite accurate. 
#' 
#' @field number_intervals := vector of number of bins to cover the filter space with (per dimension)
#' @field percent_overlap := vector of overlap percentages
#' @field quantile_method := the method to use to compute the quantiles
#' @author Matt Piekenbrock
#' @export
AdaptiveCover <- R6::R6Class(
  classname = "AdaptiveCover",
  inherit = CoverRef,
  public = list(number_intervals=NA, percent_overlap=NA, quantile_method=NA)
)

AdaptiveCover$set("public", "initialize", function(filter_values, ...){
  super$initialize(filter_values, typename="adaptive cover")
  params <- list(...)
  if ("number_intervals" %in% names(params)){ self$number_intervals <- params[["number_intervals"]] }
  if ("percent_overlap" %in% names(params)){ self$percent_overlap <- params[["percent_overlap"]] }
  if ("quantile_method" %in% names(params)){ self$quantile_method <- params[["quantile_method"]] }
})

AdaptiveCover$set("public", "construct_cover", function(){
  stopifnot(!is.na(self$number_intervals))
  stopifnot(!is.na(self$percent_overlap))
  stopifnot(!is.na(self$quantile_method))
  
  ## Setup a valid index set (via cartesian product)
  cart_prod <- arrayInd(seq(prod(self$number_intervals)), .dim = self$number_intervals)
  self$index_set <- apply(cart_prod, 1, function(x){ sprintf("(%s)", paste0(x, collapse = " ")) })
  
  ## Get filter min and max ranges
  filter_rng <- apply(self$filter_values, 2, range)
  { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
  filter_len <- diff(filter_rng)
  
  ## Construct the level set bounds
  { k <- self$number_intervals; g <- self$percent_overlap/100 } 
  interval_bnds <- lapply(seq(private$.filter_dim), function(d_i){
    { ki <- k[d_i]; gi <- g[d_i] }  ## number intervals + overlap (this dimension)
    len <- (1 / (ki - (ki-1)*gi))   ## interval length
    e <- len*(1-gi)                 ## step size 
    
    ## All intervals here are between [0, 1]
    intervals <- sapply(seq(0L, ki-1L), function(i){
      c(e*i, e*i + len) + c(-1, 1)*sqrt(.Machine$double.eps)
    })
    
    ## Convert the proportions to quantiles based on data density 
    t(apply(intervals, 2, function(bnds){
      quantile(self$filter_values[, d_i], probs = bnds, type = self$quantile_method)
    }))
  })
  interval_bnds <- do.call(cbind, interval_bnds) 
  
  ## Construct the level sets
  self$level_sets <- constructIsoAlignedLevelSets(self$filter_values, as.matrix(interval_bnds))
  
  ## Always return self 
  invisible(self)
})
