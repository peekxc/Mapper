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
landmarks <- function(x, n=NULL, eps=NULL, k=NULL, dist_method = "euclidean", seed_index = 1, shuffle_data=FALSE){
  stopifnot(is.matrix(x))
  stopifnot(seed_index >= 1L && seed_index <= nrow(x))
  stopifnot(!is.null(n) || !is.null(eps) || !is.null(k)) # must specify a number of balls, a radius, or k-neighborhood

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
    stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))

    # STEP 1: Pick point in the space (seed) and add it to the list of centers/landmarks
    C = list(seed_index)
    f_C = list(x[seed_index])

    # STEP 2: Compute distance between landmark set and each point in the space
    dists = sapply(f_C, function(c) {
      proxy::dist(c, x, method = dist_method)
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
          proxy::dist(c, x, method = dist_method)
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
  } else if (!is.null(k)){
    stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))

    # STEP 1: Pick point in the space (seed) and add it to the list of centers/landmarks
    C = list(seed_index)
    f_C = list(x[seed_index])

    # STEP 2: Compute distance between landmark set and each point in the space
    dists = sapply(f_C, function(c) {
      proxy::dist(c, x, method = dist_method)
    })
    max = which.max(dists)
    d = dists[max]

    k_nhds = list(apply(dists,2,order)[1:k])
    pt_list = k_nhds[[1]]

    # Continue if distance is greater than epsilon
    if(!(max %in% pt_list)){
      # STEP 2: Compute distance between landmark set and each point in the space
      while(TRUE){
        dists = sapply(pt_list, function(c) {
          proxy::dist(x[c], x, method = dist_method)
        })
        orderedIndices = t(apply(dists,1, sort))
        max = which.max(orderedIndices[,1])
        d = orderedIndices[max,1]

        # Continue until all points are within a k-neighborhood
        if(!(max %in% pt_list)){
          C = append(C,max)
          new_pts = order(proxy::dist(x[max], x, method = dist_method))[1:k]
          k_nhds = append(k_nhds, list(new_pts))
          pt_list = append(pt_list, new_pts)
        } else{ break }
      }
    }
    idx_list = as.character(unlist(C))
    landmark_idx = structure(k_nhds, names=idx_list)
  }

  if (is.na(shuffle_idx)){ landmark_idx } else { shuffle_idx[landmark_idx] }
}
