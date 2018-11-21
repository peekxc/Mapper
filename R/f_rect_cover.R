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
FixedRectangularCover <- R6::R6Class("FixedRectangularCover",
  inherit = CoverRef,
  private = list(.number_intervals=NA, .percent_overlap=NA)
)

#' @export
FixedRectangularCover$set("public", "initialize", function(filter_values, ...){
  super$initialize(filter_values, typename="Fixed Rectangular")
  params <- list(...)
  if ("number_intervals" %in% names(params)){ self$number_intervals <- params[["number_intervals"]] }
  if ("percent_overlap" %in% names(params)){ self$percent_overlap <- params[["percent_overlap"]] }
})

## Set overlap/gain threshold
FixedRectangularCover$set("active", "percent_overlap",
  function(value){
    if (missing(value)){ private$.percent_overlap }
    else {
      if (any(value < 0) || any(value >= 100)){ stop("The percent overlap must be a percentage between [0, 100).") }
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
FixedRectangularCover$set("active", "number_intervals", 
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

FixedRectangularCover$set("public", "format", function(...){
  # type_pretty <- paste0(toupper(substr(self$typename, start = 1, stop = 1)), tolower(substr(self$typename, start = 2, stop = nchar(self$typename))))
  sprintf("Cover: (typename = %s, number intervals = [%s], percent overlap = [%s]%%)",
          private$.typename,
          paste0(private$.number_intervals, collapse = ", "),
          paste0(format(private$.percent_overlap, digits = 3), collapse = ", "))
})

## 
## This Function is specific to the rectangular-type covers
FixedRectangularCover$set("public", "level_set_bounds", function(){
  stopifnot(!is.na(private$.percent_overlap))
  stopifnot(!is.na(private$.number_intervals))
  
  ## LS Multi-indices == Cartesian product of the intervals
  cart_prod <- arrayInd(seq(prod(self$number_intervals)), .dim = self$number_intervals)
 
  ## Get filter min and max ranges
  filter_rng <- apply(self$filter_values, 2, range)
  { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
  filter_len <- diff(filter_rng)
  
  ## Construct the interval bounds for each level set 
  prop_overlap <- self$percent_overlap/100
  base_interval_length <- filter_len/self$number_intervals
  interval_length <- base_interval_length + (base_interval_length * prop_overlap)/(1.0 - prop_overlap)
  eps <- (interval_length/2.0) + sqrt(.Machine$double.eps) ## ensures each point is in the cover
  ls_bounds <- t(apply(cart_prod, 1, function(idx){
    centroid <- filter_min + ((as.integer(idx)-1L)*base_interval_length) + base_interval_length/2.0
    c(centroid - eps, centroid + eps)
  }))
  
  ## Return 
  return(ls_bounds)
})

## Given the current set of parameter values, construct the level sets whose union covers the filter space
FixedRectangularCover$set("public", "construct_cover", function(...){
  stopifnot(!is.na(private$.percent_overlap))
  stopifnot(!is.na(private$.number_intervals))
  
  ## Setup a valid index set (via cartesian product)
  cart_prod <- arrayInd(seq(prod(self$number_intervals)), .dim = self$number_intervals)
  self$index_set <- apply(cart_prod, 1, function(x){ sprintf("(%s)", paste0(x, collapse = " ")) })
  
  ## Retrieve the level set bounds
  ls_bnds <- self$level_set_bounds()
  self$level_sets <- constructIsoAlignedLevelSets(self$filter_values, as.matrix(ls_bnds))
  
  ## Always return self 
  invisible(self)
})

FixedRectangularCover$set("public", "level_sets_to_compare", function(){
  # browser()
  all_pairs <- t(combn(1L:length(private$.index_set), 2))
  multi_index <- arrayInd(seq(prod(self$number_intervals)), .dim = self$number_intervals)
  
  ## Get filter min and max ranges
  filter_rng <- apply(self$filter_values, 2, range)
  { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
  filter_len <- diff(filter_rng)
  d_rng <- 1:ncol(self$filter_values)
  
  ## Compute the critical distances wherein if the 
  base_interval_length <- filter_len/self$number_intervals
  prop_overlap <- self$percent_overlap/100
  critical_dist <- lapply(d_rng, function(d_i) { base_interval_length[d_i] + ((base_interval_length[d_i]/2)*seq(1, self$number_intervals[d_i] - 1))*2 })
  c_interval_length <- base_interval_length + (base_interval_length * prop_overlap)/(1.0 - prop_overlap)
  
  ## Get the maximum index deviation allowed between level sets
  max_dev <- sapply(d_rng, function(d_i) { findInterval(c_interval_length[d_i], critical_dist[[d_i]])+1L })
  
  ## Filter based on percent overlap  
  which_pairs <- apply(all_pairs, 1, function(ls_pair){
    m1 <- multi_index[ls_pair[1],]
    m2 <- multi_index[ls_pair[2],]
    all(sapply(d_rng, function(d_i){ abs(m1[d_i] - m2[d_i]) <= max_dev[d_i] }))
  })
  
  ## Return the bounded pairs to compute
  res <- apply(all_pairs[which_pairs,], 2, function(x) { self$index_set[x] })
  return(res)
})

## Converts percent overlap to interval length for a fixed number of intervals
FixedRectangularCover$set("public", "overlap_to_interval_len", function(percent_overlap){
  stopifnot(all(is.numeric(self$number_intervals)))
  filter_rng <- apply(self$filter_values, 2, range)
  { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
  filter_len <- diff(filter_rng)
  base_interval_length <- filter_len/self$number_intervals
  prop_overlap <- percent_overlap/100
  return(base_interval_length + (base_interval_length*prop_overlap)/(1.0 - prop_overlap))
})

