#' Fixed Rectangular Cover
#'
#' @docType class
#' @description The rectangular cover is multidimensional, two-parameter family of covers. Given the number of
#' intervals and the overlap percentage between these intervals, this class constructs hyper-rectangular cover of the
#' given filter space, where individual level sets are distributed uniformly along a rectangular grid.
#' @field number_intervals := vector of number of bins to cover the Z with (per dimension)
#' @field percent_overlap := vector of overlap percentages
#' @author Matt Piekenbrock
#' @export
RestrainedRectangularCover <- R6::R6Class("RestrainedRectangularCover",
  inherit = CoverRef,
  private = list(.number_intervals=NA, .percent_overlap=NA), 
  lock_objects = TRUE
)

#' @export
RestrainedRectangularCover$set("public", "initialize", function(filter_values, ...){
  super$initialize(filter_values, type="Restrained Rectangular")
  params <- list(...)
  if ("number_intervals" %in% names(params)){ self$number_intervals <- params[["number_intervals"]] }
  if ("percent_overlap" %in% names(params)){ self$percent_overlap <- params[["percent_overlap"]] }
})

## Set overlap/gain threshold
RestrainedRectangularCover$set("active", "percent_overlap",
  function(value){
    if (missing(value)){ private$.percent_overlap }
    else {
      if (any(value < 0) || any(value > 1)){ stop("The percent overlap must be a percentage between [0, 1].") }
      if (length(value) != private$.filter_dim && length(value) != 1){ stop("The percent overlap must be a single scalar or a vector of scalars with length equal to the dimensionality of the filter space.") }
      if (length(value) == 1 && private$.filter_dim > 1){ value <- rep(value, private$.filter_dim) } ## create a vector
      private$.percent_overlap <- value
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
      if (length(value) == 1 && private$.filter_dim > 1){ value <- rep(value, private$.filter_dim) } ## create a vector
      stopifnot(all(value > 0))
      stopifnot(length(value) == private$.filter_dim)
      private$.number_intervals <- value
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
RestrainedRectangularCover$set("public", "construct_cover", function(){
  stopifnot(!is.na(private$.percent_overlap))
  stopifnot(!is.na(private$.number_intervals))
  
  ## Setup a valid index set (e.g. cartesian product)
  indices <- lapply(self$number_intervals, function(k) seq(k)) ## per-dimension possible indexes
  cart_prod <- as.matrix(do.call(expand.grid, structure(indices, names = paste0("d", 1:private$.filter_dim))))
  self$index_set <- apply(cart_prod, 1, function(x){ sprintf("(%s)", paste0(x, collapse = " ")) })
  
  ## Get filter min and max ranges
  filter_rng <- apply(self$filter_values, 2, range)
  { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
  filter_len <- diff(filter_rng)
  
  ## Construct the level sets
  # browser()
  { k <- self$number_intervals; g <- self$percent_overlap } 
  r <- filter_len/(k - g*(k - 1L))
  e <- r * (1 - g)
  ls_bnds <- t(apply(cart_prod, 1, function(idx){
    s <- filter_min + (as.integer(idx)-1L) * e
    c(s, s + r)
  }))
  self$level_sets <- constructIsoAlignedLevelSets(self$filter_values, as.matrix(ls_bnds))
  invisible(self)
})
