#' Cover abstract class
#' @aliases cover
#' @description Reference Class (R6) implementation of a Cover. This class is meant to act as an abstract class to derive other 
#' types of covering generators with. Minimally, a derived covering class must implement the
#' 'construct_cover' method to populate the 'level_sets' list with point indices, and any parameters 
#' that the derived class requires. 
#' 
#' Additional methods may also be added to improve the efficiency of the cover. 
#' See the \href{https://peekxc.github.io/Mapper/articles/UsingCustomCover.html}{vignette} on creating a custom 
#' cover for details. 
#' 
#' @section Fields:
#' The following is a list of the fields available for derived classes. Each may be accessed 
#' by the \code{self} environment. See \code{?R6} for more details. 
#' \itemize{
#'  \item{\emph{typename}:}{ name of the covering method.}
#'  \item{\emph{sets}:}{ named list, indexed by \code{index_set}, representing the collection of sets.}
#'  \item{\emph{index_set}:}{ character vector of keys that uniquely index the \code{sets} of the cover.}
#' }
#' @format An \code{\link{R6Class}} generator object
#' @author Matt Piekenbrock
#' @export CoverRef
CoverRef <- R6::R6Class("CoverRef", 
  private = list(
    .typename = character(0),
    .sets = NULL
  ), 
  lock_class = FALSE,  ## Feel free to add your own members
  lock_objects = FALSE ## Or change existing ones 
)

## Cover initialization
CoverRef$set("public", "initialize", function(typename){
  private$.typename <- typename
})

## Typename field
## typename ----
CoverRef$set("active", "typename", function(value){
  if (missing(value)){ private$.typename } else {
    stop("Cover 'typename' member is read-only.")
  }
})

## The level sets must be a list indexed by the index set. If the list is named, a check is performed to make sure the 
## names match the values of the index set, and in the proper order. Otherwise, the order is assumed to be correct. 
## sets ----
CoverRef$set("active", "sets", 
  function(value){
    if (missing(value)){
      private$.sets
    } else {
      stopifnot(is.list(value), all(is.character(names(value))), length(names(value)) == length(value))
      private$.sets <- value
    }
  }
)

## The index set may be composed of any data type, but the collection of indices must uniquely
## index the level sets list via the `[[` operator.
## index_set ----
CoverRef$set("active", "index_set", function(value){
  if (missing(value)){ names(private$.sets) } 
  else {
    stopifnot(is.vector(value), length(value) == length(self$sets), all(is.character(value)))
    names(private$.set) <- value
  # tmp <- structure(vector("list", length = length(value)), names = value)
  # stopifnot(length(unique(names(tmp))) == length(value))
  # private$.index_set <- names(tmp)
  }
})

## The order of a cover is the smallest number n such that
## each point of the space belongs to at most n sets in the cover.
CoverRef$set("active", "order", function(value){
  
})

## Which indices of the index set should be compared in constructing the k-simplices? 
## This can be customized based on the cover to (dramatically) reduce 
## the number of intersection checks needed to generate the k-skeletons, where k >= 1. 
## Defaults to every pairwise combination of level sets. 
## neighborhood ----
CoverRef$set("public", "neighborhood", function(filter, k=1){
  # if (length(private$.index_set) <= 1L){ return(NULL) }
  self$index_set 
  # k_combs <- t(combn(length(private$.index_set), k+1))
  # relist(private$.index_set[k_combs], k_combs)
})

## Validates that the constructed cover is indeed a valid cover. 
## validate ----
CoverRef$set("public", "validate", function(){
  if ( length(private$.index_set) != length(private$.level_sets) ){
    stop("Cover invalid: length of the index set does nto match length of the level sets")
  }
  if ( all(!names(private$.level_sets) %in% private$.index_set) ) {
    stop("Cover invalid: Not all the level sets have a corresponding index in the index set.")
  }
  idx <- unique(unlist(private$.level_sets))
  # if ( length(idx) != private$.filter_size ){
  #   stop("Cover invalid: Not all point indices account for in the open sets!")
  # }
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
    sprintf(line_format, "fixed interval", "FixedIntervalCover", paste0(c("number_intervals", "percent_overlap"), collapse = ", ")), 
    sprintf(line_format, "restrained interval", "RestrainedIntervalCover", paste0(c("number_intervals", "percent_overlap"), collapse = ", ")),
    # sprintf(line_format, "adaptive", "AdaptiveCover", paste0(c("number_intervals", "percent_overlap", "quantile_method"), collapse = ", ")),
    sprintf(line_format, "ball", "BallCover", paste0("epsilon", collapse = ", "))
  ))
}

# TODO: add ability to convert cover to new type to preserve copying the filter points multiple times, or 
# move the filter points to the MapperRef object 

