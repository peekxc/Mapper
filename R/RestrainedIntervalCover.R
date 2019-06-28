#' Restrained Interval Cover
#'
#' @docType class
#' @description The interval cover is multidimensional, two-parameter family of covers. Given the number of
#' intervals and the overlap percentage between these intervals, this class constructs hyper-rectangular cover of the
#' given filter space, where individual level sets are distributed uniformly along a rectangular grid.
#' @field number_intervals vector of number of bins to cover the Z with (per dimension)
#' @field percent_overlap vector of overlap percentages
#' @author Matt Piekenbrock
#' @export
RestrainedIntervalCover <- R6::R6Class("RestrainedIntervalCover",
  inherit = CoverRef,
  private = list(.number_intervals=NULL, .percent_overlap=NULL), 
  lock_objects = TRUE
)

#' @export
RestrainedIntervalCover$set("public", "initialize", function(...){
  super$initialize(typename="Restrained Interval")
  params <- list(...)
  if ("number_intervals" %in% names(params)){ self$number_intervals <- params[["number_intervals"]] }
  if ("percent_overlap" %in% names(params)){ self$percent_overlap <- params[["percent_overlap"]] }
})

## Set overlap/gain threshold
## percent_overlap ----
RestrainedIntervalCover$set("active", "percent_overlap",
  function(value){
    if (missing(value)){ private$.percent_overlap }
    else {
      if (any(value < 0) || any(value >= 100)){ stop("The percent overlap must be a percentage between [0, 100).") }
      private$.percent_overlap <- value
      self
    }
  }
)

## Active binding to set the number of intervals to distribute along each dimension. 
## By default, if a scalar is given and the filter dimensionality is > 1, the scalar is 
## repeated along each dimension. 
## number_intervals ----
RestrainedIntervalCover$set("active", "number_intervals",
  function(value){
    if (missing(value)){ private$.number_intervals }
    else {
      stopifnot(all(value > 0))
      private$.number_intervals <- value
      self
    }
  }
)

## Validates the parameter settings
## validate ----
RestrainedIntervalCover$set("public", "validate", function(filter){
  stopifnot(!is.null(private$.percent_overlap))
  stopifnot(!is.null(private$.number_intervals))
  stopifnot(all(self$number_intervals > 0))
  stopifnot(all(self$percent_overlap >= 0), all(self$percent_overlap < 100))
  fv <- filter()
  f_dim <- ncol(fv)
  if (length(self$number_intervals) == 1 && f_dim > 1){
    self$number_intervals <- rep(self$number_intervals[1], f_dim)
  }
  if (length(self$percent_overlap) == 1 && f_dim > 1){
    self$percent_overlap <- rep(self$percent_overlap[1], f_dim)
  }
})


## format ----
RestrainedIntervalCover$set("public", "format", function(...){
  titlecase <- function(x){
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }
  sprintf("%s Cover: (number intervals = [%s], percent overlap = [%s]%%)",
          titlecase(private$.typename),
          paste0(self$number_intervals, collapse = ", "),
          paste0(format(self$percent_overlap, digits = 3), collapse = ", "))
})

## Setup a valid index set (via cartesian product)
## construct_index_set ----
RestrainedIntervalCover$set("public", "construct_index_set", function(...){
  cart_prod <- arrayInd(seq(prod(self$number_intervals)), .dim = self$number_intervals)
  self$index_set <- apply(cart_prod, 1, function(x){ sprintf("(%s)", paste0(x, collapse = " ")) })
})


## Given the current set of parameter values, construct the level sets whose union covers the filter space
## construct_cover ----
RestrainedIntervalCover$set("public", "construct_cover", function(filter, index=NULL){
  stopifnot(is.function(filter))
  self$validate(filter)
  
  ## Get filter values 
  fv <- filter()
  f_dim <- ncol(fv)
  
  ## If the index set hasn't been made yet, construct it.
  self$construct_index_set()
  
  ## Get filter min and max ranges
  filter_rng <- apply(fv, 2, range)
  { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
  filter_len <- diff(filter_rng)
  
  ## Construct the level sets
  { k <- self$number_intervals; g <- self$percent_overlap/100 } 
  r <- filter_len/(k - g*(k - 1L))
  e <- r * (1 - g)
  cart_prod <- arrayInd(seq(prod(self$number_intervals)), .dim = self$number_intervals)
  ls_bnds <- t(apply(cart_prod, 1, function(idx){
    s <- filter_min + (as.integer(idx)-1L) * e
    c(s, s + r)
  }))
  self$level_sets <- constructIsoAlignedLevelSets(fv, as.matrix(ls_bnds))
  if (!missing(index)){ return(self$level_sets[[index]]) }
  invisible(self)
})
