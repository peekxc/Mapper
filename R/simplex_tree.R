#' Simplex Tree
#' @description Simplex tree class exposed as an Rcpp Module. 
#' @details A simplex tree is a trie specialized for storing simplicial complexes. 
#' @field n_vertices Number of 0-simplexes to initialize the simplex tree with.
#' @return A queryable simplex tree. 
#' @references Mehta, Dinesh P., and Sartaj Sahni. "Chapter 18: Interval, Segment, Range, and Priority Search Trees." Handbook of data structures and applications. Chapman and Hall/CRC, 2004.
#' @import Rcpp
#' @importFrom Rcpp evalCpp Module cpp_object_initializer
#' @importFrom methods new
#' @export
simplex_tree <- function(n_vertices){
  # print(x)
  # st_module <- Module("st_module", PACKAGE = "Mapper")
  # return(new(st_module$segment_tree, x))
  return(new(SimplexTree))
}

## Load the exported SegmentTree class into the package namespace
Rcpp::loadModule("simplex_tree_module", TRUE)