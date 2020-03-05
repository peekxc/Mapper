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
landmarks <- function(x, n=NULL, eps=NULL, dist_method = "euclidean", seed_index = 1, shuffle_data=FALSE){
  stopifnot(is.matrix(x))
  stopifnot(seed_index >= 1L && seed_index <= nrow(x))
  stopifnot(!is.null(n) || !is.null(eps)) # must specify either a number of balls or a radius

  shuffle_idx <- NA
  if (shuffle_data){
    shuffle_idx <- sample(seq(nrow(x)))
    x <- x[,,drop=FALSE]
  }

  if(!is.null(n)){
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
  } else if(!is.null(eps)){
    # TODO: Currently written for euclidean distance in 1-dimensional lens space
    stopifnot(missing(dist_method) || toupper(dist_method) == "EUCLIDEAN")

    # STEP 1: Pick point in the space (seed) and add it to the list of centers/landmarks
    C = list(seed_index)
    f_C = list(x[seed_index])

    # STEP 2: Compute distance between landmark set and each point in the space
    dists = sapply(f_C, function(c) {
      abs(c-x)
    })
    max = which.max(dists)
    d = dists[max]

    # Continue if distance is greater than epsilon
    if(d >= eps){
      C = append(C, max)
      f_C = append(f_C, x[max])

      # STEP 2: Compute distance between landmark set and each point in the space
      while(TRUE){
        dists = sapply(f_C, function(c) {
          abs(c-x)
        })

        orderedIndices = t(apply(dists,1, sort))
        max = which.max(orderedIndices[,1])
        d = orderedIndices[max,1]

        # Continue until distance is less than epsilon (i.e. stop when all points are contained within an epsilon-ball)
        if(d >= eps){
          C = append(C,max)
          f_C = append(f_C,x[max])
        } else{ break }
      }
    }
    landmark_idx = unlist(C)
  }

  if (is.na(shuffle_idx)){ landmark_idx } else { shuffle_idx[landmark_idx] }
}
