#' Subspace Cover
#'
#' @docType class
#' @description The subspace cover is a special type of cover to represent subspace topologies. 
#' That is, if \eqn{X} is a topological space with topology \eqn{\tau}, and \eqn{Y} is a 
#' subset of \eqn{X}, the collection
#' \deqn{\delta = { Y ⋂ U | U ∈ \tau }}
#' is a topology on \eqn{Y}, called the \emph{subspace topology}. The open sets of the subspace 
#' topology consist of the intersections between the open sets of \eqn{X} with \eqn{Y}.\cr 
#' \cr 
#' The subspace cover represents subspaces by a vector of integer indices between \eqn{[1, n]}, where 
#' \eqn{n} is the number points in \eqn{X}. This specific class allows the creation of two types of subspaces: 
#' 1)  
#' @field number_intervals := vector of number of bins to cover the Z with (per dimension)
#' @field percent_overlap := vector of overlap percentages
#' @author Matt Piekenbrock
#' @export
SubspaceCover <- R6::R6Class(
  "SubspaceCover",
  inherit = CoverRef,
  private = list(.cover=NA, .subset_idx=NA), 
  lock_class = FALSE,  ## Feel free to add your own members
  lock_objects = FALSE ## Or change existing ones 
)

#' @export
SubspaceCover$set("public", "initialize", function(cover, ...){
  # browser()
  priv_env <- environment(cover$initialize)[["private"]]
  public_env <- environment(cover$initialize)[["self"]]
  public_sym <- Filter(function(x) !(x %in% c("clone", ".__enclos_env__", "initialize", "typename", "filter_values")), ls(public_env, all.names=TRUE))
  for(e in ls(priv_env, all.names=TRUE)) { assign(e, get(e, priv_env), private) }
  for(e in public_sym) { 
    if (e %in% ls(self)){ unlockBinding(sym = e, env = self) } 
    assign(x = e, value = get(e, public_env), envir = self) 
    if (is.function(get(e, public_env))){
      environment(self[[e]]) <- environment(self$initialize)
    }
  }
  # self$cover <- cover
  # as.environment(as.list(priv_env, all.names=TRUE))
  # ls(environment(wut$initialize)[["private"]], all.names = TRUE)
  # private$.filter_size <-  cover$.__enclos_env__[["private"]]$.filter_size 
  # private$.filter_dim <- cover$.__enclos_env__[["private"]]$.filter_dim 
  # private$.typename <- cover$.__enclos_env__[["private"]]$.typename 
})

## Set overlap/gain threshold
SubspaceCover$set("active", "cover",
  function(value){
    if (missing(value)){ private$.cover }
    else {
      stopifnot(is(value, "CoverRef"))
      private$.cover <- value
      self
    }
  }
)

## The subspace is represented as a vector of integers forming a subset of
## the 1-based indices used in the parent cover 
SubspaceCover$set("active", "subspace",
  function(value){
    if (missing(value)){ private$.subset_idx }
    else {
      stopifnot(all(value %in% 1L:private$.filter_size))
      private$.subset_idx <- value
    #   browser()
    #   rm(list = "filter_values", envir = self, inherits = TRUE)
    #   makeActiveBinding("filter_values", function(value){
    #     if (missing(value)){ self$filter_values[private$.subset_idx,,drop=FALSE] }
    #     else {
    #       self$filter_values <- value
    #       # stop("filter values in subspace cover are read-only.")
    #     }
    #   }, env = self)
    #   self
    }
  }
)

SubspaceCover$set("active", "filter_values", function(value){
  if (missing(value)){  return(self$filter_values) } 
  else {
    if (is.null(dim(value))){ value <- matrix(values) }
    if (nrow(value) )
    self$filter_values <- value
    
  }

})

SubspaceCover$set("public", "format", function(...){
  paste0("Subspace ", format(private$.cover))
})

## Maybe just inherit from prior cover type
## Given the current set of parameter values, construct the level sets whose union covers the filter space
# SubspaceCover$set("public", "construct_cover", function(...){
#   stopifnot(!is.na(private$.cover))
#   stopifnot(!is.na(private$.subset_idx))
#   if (is.null(self$level_sets)){ private$.cover$construct_cover(...) }
#   self$level_sets <- lapply(self$level_sets, function(ls_idx){
#     intersect(ls_idx, private$.subset_idx)
#   })
#   # do.call(private$.cover$construct_cover, list(), envir = self)
#   invisible(self)  ## Always return self 
# })

