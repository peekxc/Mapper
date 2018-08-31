#' Fixed Rectangular Cover
#'
#' @docType class
#' @description The rectangular cover is multidimensional, two-parameter family of covers. Given the number of
#' intervals and the overlap percentage between these intervals, this class constructs hyper-rectangular cover of the
#' given filter space, where individual level sets are distributed uniformly along a rectangular grid.
#' @field number_intervals := vector of number of bins to cover the Z with (per dimension)
#' @field percent_overlap := vector of overlap percentages
#'
#' @author Matt Piekenbrock
#' @export
RestrainedRectangularCover <- R6Class("RestrainedRectangularCover",
  inherit = CoverRef,
  private = list(.number_intervals=NA, .percent_overlap=NA), 
  lock_objects = FALSE
)

#' @export
RestrainedRectangularCover$set("public", "initialize", function(number_intervals, percent_overlap, ...){
  private$number_intervals <- number_intervals
  private$percent_overlap <- percent_overlap
  constructCover() # Construct the actual cover
})

## Set overlap/gain threshold
RestrainedRectangularCover$set("active", "percent_overlap",
  function(value){
    if (missing(value)){ .private$percent_overlap }
    else {
      if (any(value < 0) || any(value > 1)){ stop("The percent overlap must be a percentage between [0, 1].") }
      if (length(value) != .filter_dim && length(value) != 1){ stop("The percent overlap must be a single scalar or a vector of scalars with length equal to the dimensionality of the filter space.") }
      if (length(value) == 1 && .filter_dim > 1){ value <- rep(value, .filter_dim) } ## create a vector
      private$.overlap <- value
      self
    }
  }
)

## Active binding to set the number of intervals to distribute along each dimension. 
## By default, if a scalar is given and the filter dimensionality is > 1, the scalar is 
## repeated along each dimension. 
RestrainedRectangularCover$set("active", "number_intervals",
  function(value){
    if (missing(value)){ private$.number_intervals }
    else {
      if (length(value) == 1 && .filter_dim > 1){ value <- rep(value, .filter_dim) } ## create a vector
      stopifnot(all(value > 0))
      stopifnot(length(value) == .filter_dim)
      number_intervals <- value
      self
    }
  }
)

RestrainedRectangularCover$set("public", "format", function(...){
  type_pretty <- paste0(toupper(substr(self$type, start = 1, stop = 1)), tolower(substr(self$type, start = 2, stop = nchar(self$type))))
  sprintf("Cover: (type = %s, number intervals = [%s], overlap = [%s])",
          type_pretty,
          paste0(private$.number_intervals, collapse = ", "),
          paste0(format(private$.percent_overlap, digits = 3), collapse = ", "))
})

## Given the current set of parameter values, construct the level sets whose union covers the filter space
FixedRectangularCover$set("public", "constructCover", function(){
  stopifnot(!is.na(private$.percent_overlap))
  stopifnot(!is.na(private$.number_intervals))
  
  ## Setup unexpanded and full indexing set (the full indexing set is the cartesian product of each unexpanded index set)
  indices <- lapply(self$number_intervals, function(k_l) 1:k_l) ## per-dimension possible indexes
  self$index_set <- structure(eval(as.call(append(quote(expand.grid), indices))), names = paste0("d", 1:private$.filter_dim)) ## full indexing set
  
  ## Get filter min and max ranges
  filter_rng <- apply(filter_values, 2, range)
  { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
  filter_len <- diff(filter_rng)
  
  ## Construct the level sets
  browser()
  self$level_sets <- constructFixedLevelSets(
    filter_values = matrix(self$filter_values, ncol = private$.filter_dim),
    index_set = as.matrix(self$index_set),
    overlap = as.numeric(self$percent_overlap),
    number_intervals = as.numeric(self$number_intervals),
    filter_range = matrix(filter_rng, ncol = private$.filter_dim),
    filter_len = as.numeric(filter_len)
  )
  invisible(self)
})
