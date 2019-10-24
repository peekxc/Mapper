#' landmarks
#' @description Computes a set of \code{n} landmarks of a point cloud.
#' @details This function uses the greedy maxmin procedure to produce a set of \code{n} evenly spaced landmark points in the data set 
#' \code{x}. Maxmin is a simple greedy algorithm that is relatively efficient, but it has a tendency to pick out extremal points. 
#' If the distance metric is euclidean, an efficient Rcpp implementation is used. If another metric is requested, 
#' the algorithm is performed in R.  
#' @param x a data matrix. 
#' @param n the number of landmarks requested. 
#' @param dist_method the distance metric to use. Any distance measure in the \code{proxy} package is supported.
#' @param seed_index the first landmark to seed the algorithm. 
#' @references De Silva, Vin, and Gunnar E. Carlsson. "Topological estimation using witness complexes." SPBG 4 (2004): 157-166.
#' @export
landmarks <- function(x, n, eps, dist_method = "euclidean", seed_index = 1){
  stopifnot(is.matrix(x))
  stopifnot(seed_index >= 1L && seed_index <= nrow(x))
  if (!xor(missing(n), missing(eps))){ stop("Either n or eps can be provided, but not both."); }
  
  ## Choose whether to fix n or fix eps
  use_n <- missing(eps)
  
  ## Choose the distance metric
  if (is.character(dist_method)){
    dist_choice <- switch(tolower(dist_method), "euclidean"=0L, "manhattan"=1L, "maximum"=2L, 0L)
    if (use_n){ idx <- maxmin_n(x, as.integer(n), metric=dist_choice, seed=seed_index-1L) }
    else { idx <- maxmin_eps(x, as.numeric(eps), metric=dist_choice, seed=seed_index-1L) }
    return(idx)
  } else if (is.function(dist_method)){
    if (use_n){ idx <- maxmin_n_f(x, as.integer(n), dist_f=dist_method, seed=seed_index-1L) }
    else { idx <- maxmin_eps_f(x, as.numeric(eps), dist_f=dist_method, seed=seed_index-1L) }
  } else { stop(sprintf("Unsupported distance method passed: %s\n", dist_method)) }
}