#' Cover abstract class
#' 
#' @description Reference Class (R6) implementation of a Cover. This class is meant to act as an abstract class to derive other 
#' types of covering generators with. Minimally, a derived covering class must implement the
#' 'construct_cover' method to populate the 'level_sets' list with point indices, and any parameters 
#' that the derived class requires. 
#' 
#' Additional methods may also be added to improve the efficiency of the cover. See the vignette on creating a custom 
#' cover for details. 
#' 
#' @section Private variables:
#' The following is a list of private variables available for derived classes. Each may be accessed 
#' by the \code{private$} access method, see \code{?R6} for more details. 
#' \itemize{
#'  \item{.level_sets}{names list of indices in the original data set to cluster over.}
#'  \item{.index_set}{character vector of keys that uniquely index the level sets.}
#'  \item{.filter_dim}{constant representing the filter dimension.}
#'  \item{.filter_size}{constant representing the number of points in the filter space.}
#'  \item{.typename}{string identifier of the covering method.}
#' }
#'  
#'    
#' @docType class
#' @field filter_values (n x d) matrix of filter values
#' @field typename Unique string identifier for the covering. 
#' @field index_set character vector used to index the 'level_sets' list 
#' @field level_sets list of the 
#' @format An \code{\link{R6Class}} generator object
#' 
#' @method level_sets_to_compare testing
#' 
#' @author Matt Piekenbrock
#' @export CoverRef
CoverRef <- R6::R6Class("CoverRef", 
  public = list(filter_values=NA),
  private = list(
    .level_sets=NA,
    .index_set = NA,
    .filter_size = NA,
    .filter_dim = NA, 
    .typename = character(0)
  )
)

## Cover initialization
CoverRef$set("public", "initialize", function(filter_values, typename){
  if (is.null(dim(filter_values))) { filter_values <- array(filter_values, dim = c(length(filter_values), 1)) }
  self$filter_values <- filter_values
  private$.filter_size <- dim(filter_values)[[1]]
  private$.filter_dim <- dim(filter_values)[[2]]
  private$.typename <- typename
})

CoverRef$set("public", "format", function(...){
  message <- c(sprintf("Open cover for %d objects (d = %d)", nrow(self$filter_values), private$.filter_dim))
  return(message)
})

## Typename field
CoverRef$set("active", "typename", 
  function(value){
    if (missing(value)){ private$.typename } else {
      stop("Cover 'typename' member is read-only.")
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
  all_combs <- t(combn(seq(length(private$.index_set)), 2))
  apply(all_combs, 2, function(x){ private$.index_set[x] })
})

## Validates that the constructed cover is indeed a valid cover. 
CoverRef$set("public", "validate", function(){
  if ( length(private$.index_set) != length(private$.level_sets) ){
    stop("Cover invalid: length of the index set does nto match length of the level sets")
  }
  if ( all(!names(private$.level_sets) %in% private$.index_set) ) {
    stop("Cover invalid: Not all the level sets have a corresponding index in the index set.")
  }
  idx <- unique(unlist(private$.level_sets))
  if ( length(idx) != private$.filter_size ){
    stop("Cover invalid: Not all point indices account for in the open sets!")
  }
})

# TODO
# .available_covers <- list()

#' Print available covers
#' @description Prints the covering generators available with with package, along with their parameters.
#' @export
covers_available <- function(){
  line_format <- " %-28s %-34s %-15s"
  writeLines(c(
    sprintf("Typename:%-20sGenerator:%-25sParameters:%-26s", "", "", ""),
    sprintf(line_format, "fixed rectangular", "FixedRectangularCover", paste0(c("number_intervals", "percent_overlap"), collapse = ", ")), 
    sprintf(line_format, "restrained rectangular", "RestrainedRectangularCover", paste0(c("number_intervals", "percent_overlap"), collapse = ", ")),
    sprintf(line_format, "adaptive", "AdaptiveCover", paste0(c("number_intervals", "percent_overlap", "quantile_method"), collapse = ", ")),
    sprintf(line_format, "ball", "BallCover", paste0("epsilon", collapse = ", "))
  ))
}

