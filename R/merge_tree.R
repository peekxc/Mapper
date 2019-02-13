#' merge_tree 
#' Creates a 'merge tree' 
#' @param smesh simplicial mesh describing the adjacency relations. 
#' @param h vector of 'height' values. Must be unique.
#' @export
merge_tree <- function(smesh, h){
  stopifnot(is(smesh, "Rcpp_SimplexTree"), stree$n_simplexes[1] == length(h), !anyDuplicated(h))
  merge_tree_res <- Mapper::simplex_tree()
  res <- Mapper:::construct_merge_tree2(stree$as_XPtr(), h, merge_tree_res$as_XPtr())
  return(merge_tree_res)
}
