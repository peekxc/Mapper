#' Simplex Tree
#' @description Simplex tree class exposed as an Rcpp Module. 
#' @details A simplex tree is a trie specialized for storing simplicial complexes. 
#' @return A queryable simplex tree. 
#' @references Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
#' @import Rcpp
#' @importFrom Rcpp evalCpp Module cpp_object_initializer
#' @importFrom methods new
#' @export
simplex_tree <- function(){
  return(new(SimplexTree))
}

## Load the exported SegmentTree class into the package namespace
Rcpp::loadModule("simplex_tree_module", TRUE)

SimplexTree.default_print <- setMethod("show", "Rcpp_SimplexTree", function (object) {
  max_k <- length(object$n_simplexes)
  if (max_k == 0){ cat("< empty simplex tree >\n") }
  else {
    cat(sprintf("Simplex Tree with (%s) (%s)-simplices\n", paste0(object$n_simplexes, collapse = ", "), paste0(0L:(max_k-1L), collapse = ", ")))
  }
})