#' Restrained Rectangular Cover
#'
#' Restrained rectangular cover reference class - part of Reference Class (R5) implementation for the Mapper package in R
#' @description The rectangular cover is multidimensional, two-parameter family of covers. Given the number of
#' intervals and the overlap percentage between these intervals, this class constructs hyperrectangular cover of the
#' given filter space, where individual level sets are distributed uniformly along a rectangular grid.
#' @field num_intervals := vector of number of bins to cover the Z with (per dimension)
#' @field percent_overlap := vector of overlap percentages
#'
#' @author Matt Piekenbrock
#' @export
r_rect_cover <- setRefClass("RestrainedRectangularCover", fields = c("number_intervals", "overlap"), contains = "Cover")

#' @export
r_rect_cover$methods(initialize = function(number_intervals, overlap, ...){
  callSuper(...)
  setResolution(number_intervals); setOverlap(overlap)
  constructCover() # Construct the actual cover
})

## Set overlap/gain threshold
r_rect_cover$methods(setOverlap = function(percent_overlap){
  dim_Z <- ncol(filter_values)
  if (any(percent_overlap < 0) || any(percent_overlap > 1)){ stop("The percent overlap must be a percentage between [0, 1].") }
  if (length(percent_overlap) != dim_Z && length(percent_overlap) != 1){ stop("The percent overlap must be a single scalar or a vector of scalars with length equal to the dimensionality of the filter space.") }
  if (length(percent_overlap) == 1 && dim_Z > 1){ percent_overlap <- rep(percent_overlap, dim_Z) } ## create a vector
  overlap <<- percent_overlap
})

## Set resolution
r_rect_cover$methods(setResolution = function(num_intervals){
  dim_Z <- ncol(filter_values)
  if (length(num_intervals) != dim_Z && length(num_intervals) != 1){ stop("The resolution must be a single scalar or a vector of scalars with length equal to the dimensionality of the filter space.") }
  if (length(num_intervals) == 1 && dim_Z > 1){ num_intervals <- rep(num_intervals, dim_Z) } ## create a vector
  number_intervals <<- num_intervals
})

## Given the current set of parameter values, construct the level sets whose union covers the filter space
r_rect_cover$methods(constructCover = function(){
  if ("uninitializedField" %in% class(overlap)){ stop("'overlap' must be set prior to constructing a valid cover.") }
  if ("uninitializedField" %in% class(number_intervals)){ stop("'number_intervals' must be set prior to constructing a valid cover.") }

  ## Setup unexpanded and full indexing set (the full indexing set is the cartesian product of each unexpanded index set)
  filter_dim <- ncol(filter_values) ## filter dimensionality
  indices <- lapply(number_intervals, function(k_l) 1:k_l) ## per-dimension possible indexes
  index_set <<- structure(eval(as.call(append(quote(expand.grid), indices))), names = paste0("d", 1:filter_dim)) ## full indexing set

  ## Get filter min and max ranges
  filter_rng <- apply(filter_values, 2, range)
  { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }

  ## Using the current gain, setup the cover. This computes the min and max bounds of each level set.
  interval_length <- (filter_max - filter_min)/(number_intervals - overlap * (number_intervals - 1)) ## interval length
  step_size <- interval_length * (1 - overlap) ## step size

  ## Construct the level sets
  level_sets <<- constructRestrainedLevelSets(filter_values, as.matrix(index_set), as.numeric(interval_length),  as.numeric(step_size), as.numeric(filter_min))
  # browser() ## length(unique(unlist(lapply(level_sets, function(ls) ls$points_in_level_set))))
})

## Given a single LSFI 'from' and a vector of LSFIs 'to' orthogonal to the level set mapper by 'from',
## retrieves the LSFIs that lie between 'from' and 'to', in both orthogonal and non-orthogonal directions (inclusively)
r_rect_cover$methods(getRectLSFI = function(from, to){
  filter_dim <- ncol(filter_values) ## filter dimensionality
  if (length(from) != 1 || length(to) != filter_dim){ stop("'from' must be length 1 and 'to' must be length < filter dim >") }
  sbase <- getLSMI(from)
  to_lsmi <- lapply(1:filter_dim, function(d_i) getLSMI(to[d_i]))
  for (d_i in 1:filter_dim){
    if (any(to_lsmi[[d_i]][, -d_i] != sbase[, -d_i])){ stop("'to' level sets must be orthogonal to originating set") }
  }
  to_indices <- lapply(1:filter_dim, function(d_i){ sbase[,d_i]:to_lsmi[[d_i]][,d_i] })
  target_lvl_sets <- as.matrix(eval(as.call(append(quote(expand.grid), to_indices))))
  getLSFI(target_lvl_sets)
})

r_rect_cover$methods(summary = function(){
  type_pretty <- paste0(toupper(substr(type, start = 1, stop = 1)), tolower(substr(type, start = 2, stop = nchar(type))))
  return(sprintf("Cover: (type = %s, number intervals = [%s], overlap = [%s])",
                 type_pretty,
                 paste0(number_intervals, collapse = ", "),
                 paste0(format(overlap, digits = 3), collapse = ", ")))
})

## Which level set (flat indices) are valid to compare against? When the overlap is <= 50%, it's required that the
## level sets are checked only against immediate neighbors of that set (when the index is at maximum one
## index away). If the overlap is > 50%, then check every combination.
r_rect_cover$methods(valid_pairs = function(){
  all_ls_pairs <- t(combn(1:length(level_sets), 2))
  idx_set_pairs <- as.matrix(index_set[all_ls_pairs[, 1],] - index_set[all_ls_pairs[, 2],])
  if (all(overlap <= 0.50)){
    valid_ls_pairs <- as.vector(abs(apply(idx_set_pairs, 1, max))) <= 1
    return(all_ls_pairs[valid_ls_pairs,])
  } else {
    ## TODO: Improve this, and extend mapper to compute more than the 1-skeleton efficiently
    return(all_ls_pairs)
  }
})
