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
#' @param shuffle_data whether to first randomly shuffle the data.
#' @references De Silva, Vin, and Gunnar E. Carlsson. "Topological estimation using witness complexes." SPBG 4 (2004): 157-166.
#' @export
landmarks <- function(x, n, dist_method = "euclidean", seed_index = 1, shuffle_data=FALSE){
  stopifnot(is.matrix(x))
  stopifnot(seed_index >= 1L && seed_index < nrow(x))
  shuffle_idx <- NA
  if (shuffle_data){ 
    shuffle_idx <- sample(seq(nrow(x)))
    x <- x[,,drop=FALSE] 
  }
  if (missing(dist_method) || toupper(dist_method) == "EUCLIDEAN"){
    landmark_idx <- landmark_maxmin(x, n, seed_index-1L)
  } else if (requireNamespace("proxy", quietly = TRUE)){
    stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))
    landmark_idx <- vector(mode="integer", n)
    landmark_idx[1L] <- seed_index
    landmark_min_dist <- rep(Inf, nrow(x))
    for (i in 2L:n){
      landmark_dist <- proxy::dist(x, x[landmark_idx[i-1L],,drop=FALSE], method = dist_method)
      landmark_min_dist <- pmin(landmark_dist, landmark_min_dist)
      potential_idx <- setdiff(seq(nrow(x)), landmark_idx[c(1:i)])
      landmark_idx[i] <- potential_idx[which.max(landmark_min_dist[potential_idx])]
    }
  } else {
    stop(sprintf("Unsupported distance method passed: %s\n", dist_method))
  }
  if (is.na(shuffle_idx)){ landmark_idx } else { shuffle_idx[landmark_idx] }
}