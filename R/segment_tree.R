#' Segment Tree
#' @description Fast segment tree class exposed as an Rcpp Module. 
#' @details A segment tree is rooted-augment binary search tree useful for querying information across a fixed set of intervals.
#' This segment tree class is specialized to handle interval queries which return indices of points that intersect the given query interval.
#' The query interval must be composed using endpoints which were used to construct the tree. 
#' @param intervals an (n x 2) matrix of intervals.
#' @return A queryable segment tree. 
#' @references Mehta, Dinesh P., and Sartaj Sahni. "Chapter 18: Interval, Segment, Range, and Priority Search Trees." Handbook of data structures and applications. Chapman and Hall/CRC, 2004.
#' @import Rcpp
#' @importFrom Rcpp evalCpp Module cpp_object_initializer
#' @importFrom methods new
#' @export
segment_tree <- function(intervals){
  stopifnot(is.matrix(intervals))
  stopifnot(dim(intervals)[[2]] == 2)
  return(new(SegmentTree, intervals))
}

## Segment Tree module
Rcpp::loadModule("segment_tree_module", TRUE)