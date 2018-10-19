#' Cover reference class - Reference Class (R6) implementation of a Cover
#'
#' @docType class
#' @field filter_values (n x d) matrix of filter values
#' @field type character string of the type of cover
#' @author: Matt Piekenbrock
#' @export CoverRef
CoverRef <- R6::R6Class("CoverRef", 
  public = list(filter_values=NA),
  private = list(
    .level_sets=NA,
    .index_set = NA,
    .filter_size = NA,
    .filter_dim = NA, 
    .type = character(0)
  )
)

## Cover initialization
CoverRef$set("public", "initialize", function(filter_values, type){
  if (is.null(dim(filter_values))) { filter_values <- array(filter_values, dim = c(length(filter_values), 1)) }
  self$filter_values <- filter_values
  private$.filter_size <- dim(filter_values)[[1]]
  private$.filter_dim <- dim(filter_values)[[2]]
  private$.type <- type
})

CoverRef$set("public", "format", function(...){
  browser()
  message <- c(sprintf("Open cover for %d objects (d = %d)", nrow(self$filter_values), private$.filter_dim))
  return(message)
})

## Type field
CoverRef$set("active", "type", 
  function(value){
    if (missing(value)){ private$.type } else {
      stop("Cover 'type' member is read-only.")
    }
})

## The index set may be composed of any data type, but the collection of indices must uniquely
## index the level sets list via the `[[` operator.
CoverRef$set("active", "index_set", 
 function(value){
   if (missing(value)){
     private$.index_set
   } else {
     stopifnot(is.vector(value))
     tmp <- structure(vector("list", length = length(value)), names = value)
     stopifnot(length(unique(names(tmp))) == length(value))
     private$.index_set <- names(tmp)
   }
})

## The level sets must be a list indexed by the index set. If the list is named, a check is performed to make sure the 
## names match the values of the index set, and in the proper order. Otherwise, the order is assumed to be correct. 
CoverRef$set("active", "level_sets", 
  function(value){
    if (missing(value)){
      private$.level_sets
    } else {
      stopifnot(is.list(value) && (length(value) == length(private$.index_set)))
      if (!is.null(value)){ stopifnot(all(names(value) == private$.index_set)) }
      private$.level_sets <- structure(value, names = private$.index_set)
    }
  }
)

## Default cover 
CoverRef$set("public", "construct_cover", function(){
  stop("Base class cover construction called. This method must be overridden.")
})

## Which level sets (in terms of their corresponding indices in the index set) should be compared? 
## This can be customized based on the cover to (dramatically) reduce the number of intersection checks
## needed to generate the k-skeletons, where k >= 1. Defaults to every pairwise combination of level sets. 
CoverRef$set("public", "level_sets_to_compare", function(){
  t(combn(1L:length(private$.index_set), 2))
})

