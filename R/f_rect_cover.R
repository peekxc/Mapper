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
FixedRectangularCover <- R6Class("FixedRectangularCover",
  inherit = CoverRef,
  private = list(.number_intervals=NA, .percent_overlap=NA)
)

#' @export
FixedRectangularCover$set("public", "initialize", function(filter_values, ...){
  super$initialize(filter_values, type="Fixed Rectangular")
  params <- list(...)
  if ("number_intervals" %in% names(params)){ self$number_intervals <- params[["number_intervals"]] }
  if ("percent_overlap" %in% names(params)){ self$percent_overlap <- params[["percent_overlap"]] }
})

## Set overlap/gain threshold
FixedRectangularCover$set("active", "percent_overlap",
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
  type_pretty <- paste0(toupper(substr(self$type, start = 1, stop = 1)), tolower(substr(self$type, start = 2, stop = nchar(self$type))))
  sprintf("Cover: (type = %s, number intervals = [%s], overlap = [%s])",
          type_pretty,
          paste0(private$.number_intervals, collapse = ", "),
          paste0(format(private$.percent_overlap, digits = 3), collapse = ", "))
})

## This Function is specific to the rectangular-type covers
FixedRectangularCover$set("public", "get_level_set_bounds", function(which_levels){
  stopifnot(!is.na(private$.percent_overlap))
  stopifnot(!is.na(private$.number_intervals))
  
  ## Get filter min and max ranges
  filter_rng <- apply(self$filter_values, 2, range)
  { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
  filter_len <- diff(filter_rng)
  
  ## Construct the interval bounds for each level set 
  base_interval_length <- filter_len/self$number_intervals
  interval_length <- base_interval_length + (base_interval_length * self$percent_overlap)/(1.0 - self$percent_overlap)
  eps <- interval_length/2.0
  ls_bounds <- t(apply(cart_prod, 1, function(idx){
    centroid <- filter_min + ((as.integer(idx)-1L)*base_interval_length) + base_interval_length/2.0
    c(centroid - eps, centroid + eps)
  }))
})

## Given the current set of parameter values, construct the level sets whose union covers the filter space
FixedRectangularCover$set("public", "construct_cover", function(...){
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
  base_interval_length <- filter_len/self$number_intervals
  interval_length <- base_interval_length + (base_interval_length * self$percent_overlap)/(1.0 - self$percent_overlap)
  eps <- interval_length/2.0
  ls_bnds <- t(apply(cart_prod, 1, function(idx){
    centroid <- filter_min + ((as.integer(idx)-1L)*base_interval_length) + base_interval_length/2.0
    c(centroid - eps, centroid + eps)
  }))
  self$level_sets <- constructIsoAlignedLevelSets(self$filter_values, as.matrix(ls_bnds))
  # all(which(filter_values[, 1] >= 3.247658 & filter_values[, 2] >= 3.247658 & filter_values[, 1] <= 4.295289 & filter_values[, 2] <= 4.295289) == level_sets[[25]]$points_in_level_set)

  invisible(self)
  #  ## length(unique(unlist(lapply(level_sets, function(ls) ls$points_in_level_set))))
})

## Given a single LSFI 'from' and a vector of LSFIs 'to' orthogonal to the level set mapper by 'from',
## retrieves the LSFIs that lie between 'from' and 'to', in both orthogonal and non-orthogonal directions (inclusively)
# FixedRectangularCover$methods(getRectLSFI = function(from, to){
#   filter_dim <- ncol(filter_values) ## filter dimensionality
#   if (length(from) != 1 || length(to) != filter_dim){ stop("'from' must be length 1 and 'to' must be length < filter dim >") }
#   sbase <- getLSMI(from)
#   to_lsmi <- lapply(1:filter_dim, function(d_i) getLSMI(to[d_i]))
#   for (d_i in 1:filter_dim){
#     if (any(to_lsmi[[d_i]][, -d_i] != sbase[, -d_i])){ stop("'to' level sets must be orthogonal to originating set") }
#   }
#   to_indices <- lapply(1:filter_dim, function(d_i){ sbase[,d_i]:to_lsmi[[d_i]][,d_i] })
#   target_lvl_sets <- as.matrix(eval(as.call(append(quote(expand.grid), to_indices))))
#   getLSFI(target_lvl_sets)
# })



## Which level set (flat indices) are valid to compare against? When the overlap is <= 50%, it's required that the
## level sets are checked only against immediate neighbors of that set (when the index is at maximum one
## index away). If the overlap is > 50%, then check every combination.
# FixedRectangularCover$methods(valid_pairs = function(){
#   all_ls_pairs <- t(combn(1:length(level_sets), 2))
#   idx_set_pairs <- as.matrix(index_set[all_ls_pairs[, 1],] - index_set[all_ls_pairs[, 2],])
#   if (all(overlap <= 0.50)){
#     valid_ls_pairs <- as.vector(abs(apply(idx_set_pairs, 1, max))) <= 1
#     return(all_ls_pairs[valid_ls_pairs,])
#   } else {
#     ## TODO: Improve this, and extend mapper to compute more than the 1-skeleton efficiently
#     return(all_ls_pairs)
#   }
# })
