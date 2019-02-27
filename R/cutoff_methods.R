#' cutoff_first_bin
#'
#' @description Uses the heuristic described in Section 3.1 of the original Mapper paper to cut a given hierarchical tree to produce a partitioning of the data.
#'
#' @param hcl hierarchical clustering in the form an 'hclust' object.
#' @param num_bins Controls how many bins there are in the histogram used to determine cutoff.
#' @param check_skew Whether to check if the distribution is left-skewed. See details. 
#' 
#' @details This method implements a cutting heuristic to determine the cut height in the given hierarchical tree \code{hcl}. 
#' The cut value is chosen as the lowest break point corresponding to the first empty bin of a histogram of the linkage distances.
#' The motivation for the heuristic is that its often observed empirically for 'nice' clustering situations that the majority of 
#' linkage distances representing inter-cluster distances are relatively smooth, whereas the intra-cluster distances are sparse. 
#' Intuitively, such a distribution may be thought as being right-skewed, and the first empty interval in a histogram of linkage
#' distances may be a decent splot to cut the hierarchy. If \code{check_skew} is TRUE (default), the linkage distances are checked 
#' that they are indeed right-skewed, and then the heuristic is used. If the distribution is left-skewed, the assumption for the 
#' heuristic is not true, and the trivial clustering is returned instead. 
#' 
#' @importFrom stats cutree
#' @export
cutoff_first_bin <- function(hcl, num_bins, check_skew=TRUE) {
  if (!is(hcl, "hclust")){ stop("'cutoff_first_bin' expects an 'hclust' object as input.") }
  breaks <- as.double(seq(head(hcl$height,1L), tail(hcl$height,1L), length.out = num_bins))
  bin_idx <- tabulate(findInterval(x = as.double(hcl$height), vec = breaks, rightmost.closed = FALSE, all.inside = TRUE, left.open = FALSE), nbins = num_bins)
  
  ## If the majority of the mass of linkage-distances are left-skewed (high), then the 
  ## motivation for the heuristic is not met, return everything in 1 cluster.
  if (check_skew){
    mid <- median(seq(num_bins))
    if (sum(bin_idx[1L:(mid-1L)]) < sum(bin_idx[mid:num_bins])){
      return(as.vector(cutree(hcl, k = 1L)))
    }
  }
  
  ## Check for the first empty
  cut_idx <- findFirstEqual(bin_idx, 0L) ## only need the first position
  
  ## If there's a uniform distribution of distances, assign every point to the same cluster
  if (cut_idx == length(bin_idx) + 1){ as.vector(cutree(hcl, k = 1L)) }
  return(as.vector(cutree(hcl, h = breaks[cut_idx])))
}


#' cutoff_first_threshold 
#' 
#' @description Uses a heuristic similar to the one described in Section 3.1 of the original Mapper paper to cut a given hierarchical tree to produce a partitioning of the data.
#' 
#' @param hcl hierarchical clustering in the form an 'hclust' object.
#' @param threshold density threshold.  
#' @param ... Additional parameters passed to \code{\link[stats]{density}}.
#' @details The cut value is chosen as the lowest (interpolated) height value where the density of linkage distances is lower than or equal to 
#' a given threshold value. Note that while this method offers a smooth alternative to \code{\link{cutoff_first_bin}} and may provide empirically better 
#' clusterings, it's noticeably slower due to the call the \code{\link[stats]{density}}. 
#' 
#' @importFrom stats cutree
#' @export
cutoff_first_threshold <- function(hcl, threshold = 0.0, ...){
  if (!is(hcl, "hclust")){ stop("'cutoff_first_bin' expects an 'hclust' object as input.") }
  n <- nrow(hcl$merge)+1L
  density_params <- list(...)
  if (is.null(density_params$n)){ density_params$n <- min(c(2^ceiling(log(n)/log(2)), 512L)) }
  density_params$x <- hcl$height
  f_h <- do.call(stats::density, density_params)
  cut_idx <- Position(function(x) abs(x - threshold) < sqrt(.Machine$double.eps), x = f_h$y)
  if (is.na(cut_idx) || cut_idx == 1L) { as.vector(cutree(hcl, k = 1L)) }
  else { as.vector(cutree(hcl, h = f_h$x[cut_idx])) }
}