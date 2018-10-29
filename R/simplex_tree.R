#' Simplex Tree
#' @description Simplex tree class exposed as an Rcpp Module. 
#' @details A simplex tree is a trie specialized for storing simplicial complexes. The current implementation provides 
#' a limited API and basic functionality. 
#' @return A queryable simplex tree. 
#' @references Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
#' @export
simplex_tree <- function(){
  return(new(SimplexTree))
}

setClass("Rcpp_SimplexTree")
.print_simplex_tree <- setMethod("show", "Rcpp_SimplexTree", function (object) {
  max_k <- length(object$n_simplexes)
  if (max_k == 0){ cat("< empty simplex tree >\n") }
  else {
    cat(sprintf("Simplex Tree with (%s) (%s)-simplices\n", paste0(object$n_simplexes, collapse = ", "), paste0(0L:(max_k-1L), collapse = ", ")))
  }
})

## Simplex Tree module (can be loaded anywhere)
Rcpp::loadModule("simplex_tree_module", TRUE)