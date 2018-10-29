#' cutoff_first_bin function
#'
#' Uses the heuristic from the original Mapper paper to cut a given hierarchical cluster tree at the
#' threshold corresponding to the first empty of a histogram of the linkage distances.
#'
#' @param hcl hierarchical clustering in the form an 'hclust' object.
#' @param diam largest pairwise distance between the objects used to create the hierarchy.  
#' @param num_bins Controls how many bins there are in the histogram used to determine cutoff.
#' @importFrom stats cutree
#' @export
cutoff_first_bin <- function(hcl, diam, num_bins) {
  if (!is(hcl, "hclust")){ stop("'cutoff_first_bin' expects an 'hclust' object as input.") }
  breaks <- as.double(seq(min(hcl$height), diam, length.out = num_bins))
  bin_idx <- tabulate(findInterval(x = as.double(hcl$height), vec = breaks, rightmost.closed = FALSE, all.inside = TRUE, left.open = FALSE), nbins = num_bins)
  cut_idx <- findFirstEqual(bin_idx, 0L) ## only need the first position
  ## If there's a uniform distribution of distances, assign every point to the same cluster
  if (cut_idx == length(bin_idx) + 1){ as.vector(cutree(hcl, k = 1)) }
  else { as.vector(cutree(hcl, h = breaks[cut_idx])) }
}
