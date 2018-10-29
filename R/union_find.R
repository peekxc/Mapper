#' Union Find
#' 
#' @description Exports a UnionFind data structure exposed as an Rcpp Module. 
#' @details UnionFind, also called a disjoint-set, is a data structure that allows efficient tracking and unioning of 
#' elements in terms of connected components. To use the structure, instantiate it with a given (fixed) size of components, and then 
#' union or find as needed. This implementation uses path compression to speed up successive find operations. 
#' @section Methods:
#' \itemize{
#'   \item{union}{ Unions two elements together. }
#'   \item{union_all}{ Unions all the elements in a given vector together. }
#'   \item{find}{ Find which components a given element is associated with. If the elements given doesn't exist, -1 is returned. }
#'   \item{find_all}{ Finds which component each element in a given vector is associated with. If the elements given doesn't exist, -1 is returned. }
#'   \item{print}{ Prints the connected components. }
#'   \item{connected_components}{ Retrieves the connected components as a integer vector. }
#' }
#' @return A queryable disjoint-set structure. 
#' @references Tarjan, Robert Endre. "Efficiency of a good but not linear set union algorithm." Journal of the ACM (JACM) 22.2 (1975): 215-225.
#' @export
union_find <- function(size){
  size <- as.integer(size)
  stopifnot(size > 0)
  return(new(UnionFind, size))
}

setClass("Rcpp_UnionFind")
.print_union_find <- setMethod("show", "Rcpp_UnionFind", function (object) {
  object$print()
})

## Simplex Tree module (can be loaded anywhere)
Rcpp::loadModule("union_find_module", TRUE)