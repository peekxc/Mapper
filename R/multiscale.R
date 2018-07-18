# Computes all unique Mapper constructions
# Only supports rectangular cover
multiscale <- function(mapper_obj){
  if (class(mapper_obj) != "MapperRef"){ stop("'multiscale' expects a Mapper reference object.") }
  if (mapper_obj$cover$type != "rectnagular"){ stop("'multiscale' is only compatable with rectangular covers.") }
  n <- ifelse("dist" %in% class(X), attr(X, "Size"), nrow(X))
  d <- ncol(mapper_obj$cover$filter_values)

  ## Start with the base cover
  mapper_obj$cover$setOverlap(0.0)
  mapper_obj$cover$constructCover()

  ## Get the range of the indices per dimension
  idx_rng <- apply(mapper_obj$cover$index_set, 2, range)

  ## Get the LSFI's for each point
  pt_lsfi <- vector(mode = "numeric", length = n)
  for (i in 1:length(mapper_obj$cover$level_sets)){
    ls <- mapper_obj$cover$level_sets[[i]]
    pt_lsfi[ls$points_in_level_set] <- i
  }

  ## Translate to LSMI
  pt_lsmi <- as.matrix(mapper_obj$cover$getLSMI(pt_lsfi))

  for (d_i in 1:d){

  }

  ## For each point, calculate the distance to the neighboring intervals
  mapper_obj$cover$filter_values
  for (i in 1:n){

  }
  mapper_obj$cover$getLSMI(1:10)

  mapper_obj$cover$index_set
}

## Returns whether a lsmi is on the border
is_border <- function(idx, index_set){ return(idx %in% range(index_set)) }

## Compute the interval size(s) needed for a given point p to intersect its neighboring level sets
## Accepts as inputs:
## p := point in question
## k := number of intervals
## z := range of the filter domain
## p_lsmi := the level set multi index (lsmi) of point p
## index_set := the index set of the filter space spanned by z
computeR <- function(p, k, z, p_lsmi, index_set){
  lower_lsfi <- index_set[index_set < p_lsmi]
  upper_lsfi <- index_set[index_set > p_lsmi]
  lower_is_border <- is_border(lower_lsfi, index_set)
  upper_is_border <- is_border(upper_lsfi, index_set)
  lower_R <- sapply(1:length(lower_is_border), function(i){
    is_edge_case <- lower_is_border[i]
    base_R <- (z / k)
    if (is_edge_case){ return(base_R + abs(base_R - p)) }
    else { return((-k * p + p + z)/(-k*lower_lsfi[i] + lower_lsfi[i] + 1)) }
  })
  upper_R <- sapply(1:length(upper_is_border), function(i){
    is_edge_case <- upper_is_border[i]
    base_R <- (z / k)
    if (is_edge_case){ return(base_R + abs(base_R - p)) }
    else { return(z - ((k - 1)/(upper_lsfi[i]))*p) }
  })
  c(lower_R, upper_R)
}
