#' Cover reference class - Reference Class (R6) implementation of a Cover
#'
#' @docType class
#' @field filter_values (n x d) matrix of filter values
#' @field type character string of the type of cover
#' @author: Matt Piekenbrock
#' @import R6
#' @export CoverRef
CoverRef <- R6Class("CoverRef", 
  public = list(filter_values=NA),
  private = list(
    .level_sets=NA,
    .index_set = NA,
    .filter_size = NA,
    .filter_dim = NA, 
    .type = character(0)
  )
)

## Cover initialization
CoverRef$set("public", "initialize", function(filter_values, type){
  if (is.null(dim(filter_values))) { filter_values <- array(filter_values, dim = c(length(filter_values), 1)) }
  self$filter_values <- filter_values
  private$.filter_size <- dim(filter_values)[[1]]
  private$.filter_dim <- dim(filter_values)[[2]]
  private$.type <- type
})

CoverRef$set("public", "format", function(...){
  browser()
  message <- c(sprintf("Open cover for %d objects (d = %d)", nrow(self$filter_values), private$.filter_dim))
  return(message)
})

## Type field
CoverRef$set("active", "type", 
  function(value){
    if (missing(value)){ private$.type } else {
      stop("Cover 'type' member is read-only.")
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
  t(combn(1L:length(private$.index_set), 2))
})

# ## Generic method which 
# CoverRef$set("public", "plot_cover", function(){
#   
# })


# CoverRef$methods(getLSMI = function(lsfi){
#   index_set[lsfi,]
# })

# CoverRef$methods(getLSFI = function(lsmi){
#   .valid_check <- function(lsfi){ ifelse(length(lsfi) == 0, NA, lsfi) }
#   if (is.null(lsmi) || (length(lsmi) == ncol(filter_values))){
#     .valid_check(which(apply(index_set, 1, function(ls_idx) all(ls_idx == lsmi))))
#   } else {
#     apply(lsmi, 1, function(lsmi_i){
#       .valid_check(which(apply(index_set, 1, function(ls_idx) all(ls_idx == lsmi_i))))
#     })
#   }
# })

## Constructs a base cover where the Lebesgue covering dimension = 0 (i.e. a valid cover, but where each set is disjoint)
## TODO: Add support for other types of covering methods
# CoverRef$methods(constructBaseCover = function(cover_type = c("rectangular")){
#   if ("uninitializedField" %in% class(percent_overlap)){ stop("'percent_overlap' must be set prior to constructing a valid cover.") }
#   if ("uninitializedField" %in% class(num_intervals)){ stop("'num_intervals' must be set prior to constructing a valid cover.") }
# 
#   ## Setup
#   filter_dim <- ncol(filter_values) ## filter dimensionality
#   indices <- lapply(num_intervals, function(k) 1:k) ## per-dimension possible indexes
#   index_set <<- structure(eval(as.call(append(quote(expand.grid), indices))), names = paste0("d", 1:filter_dim)) ## full indexing set
# 
#   ## Find which level set each point lies within under the case the given filter space is
#   ## zero-dimensional w.r.t. its Lebesgue covering dimension
#   points_in_level_set <- sapply(1:filter_dim, function(i) {
#     x <- filter_values[, i]
#     findInterval(x = x, seq(min(x), max(x), length.out = num_intervals[[i]]+1), all.inside = TRUE)
#   })
#   browser()
#   ## Create the level sets
#   index_set_str <- apply(index_set, 1, function(alpha) paste0(as.character(alpha), collapse = ","))
#   level_sets <<- structure(points_in_level_set, names = index_set_str)
# })

## Compute and store the interpoint distances for each subset
# CoverRef$methods(computeLevelSets = function(){
#   for (i in 1:length(level_sets)){
#     if ("dist" %in% class(X)){
#       level_sets[[i]]$dist <<- dist_subset(dist = X, idx = level_sets[[i]]$points_in_level_set)
#     } else { level_sets[[i]]$dist <<- dist(X[level_sets[[i]]$points_in_level_set,], method = "euclidean") }
#   }
# })


## Given an n x d set of box sizes / interval lengths, computes the overlap percentage corresponding to the given box,
## per dimension, fixed according to the current 'num_intervals' and filter space range settings.
# CoverRef$methods(computeOverlap = function(box_sizes){
#   if (is.null(dim(box_sizes)) && length(box_sizes != filter_dim)){ stop("'compute_overlap' expects a matrix of box sizes as input.") }
#   if (ncol(box_sizes) != ncol(filter_values)){ stop("The number of columns in 'box_sizes' must match the dimensionality of the filter space.")}
#   if (!is.matrix(box_sizes)){ box_sizes <- as.matrix(box_sizes) }
#   box_sizes <- matrix(box_sizes, ncol = filter_dim) ## ensure dimensionality is correct
# 
#   ## Get filter settings
#   filter_dim <- ncol(filter_values) ## filter dimensionality
#   filter_rng <- apply(filter_values, 2, range)
#   { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
#   filter_len <- filter_max - filter_min
# 
#   ## Compute the overlap and return
#   t_box_sizes <- t(box_sizes)
#   # browser()
#   overlap <- (filter_len - t_box_sizes*num_intervals)/(t_box_sizes*(-num_intervals + 1))
#   if (any(overlap < 0)) { warning("Detected negative overlap values. Some of the given box sizes do not cover the space under the current interval settings.") }
#   if (any(overlap > 1)) { warning("Detected overlap values > 1. Some of the given box sizes exceed the range of the filter space.") }
#   return(t(overlap))
# })

## Given the current set of filter values and a fixed number of intervals to distribute along each dimension,
## this method compute the halfspace distance associated with each point in the filter space, per dimension.
## A d-length list of n x 2 matrices are returned, each matrix at index i representing the lower and upper
## halfspace distances along dimension d_i
# CoverRef$methods(computeHalfSpaceDistances = function(){
#   if ("uninitializedField" %in% class(num_intervals)){ stop("'num_intervals' must be set prior to computing all overlap values.") }
# 
#   ## Get filter settings
#   filter_dim <- ncol(filter_values) ## filter dimensionality
#   filter_rng <- apply(filter_values, 2, range)
#   { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
#   filter_len <- filter_max - filter_min
# 
#   ## Base interval length (step size == 0), i.e. interval length corresponding to
#   ## the cover that is zero-dimensional w.r.t. its Lebesgue covering dimension
#   base_interval_length <- (filter_len)/(num_intervals) ## base interval length
# 
#   ## Find which level set each point lies within under the case the given filter space is
#   ## Note that is equivalent to integer division of each dimension by the base interval length given above,
#   ## but 'findInterval' is a bit safer for handling edge cases
#   original_level_sets <- sapply(1:filter_dim, function(i) {
#     x <- filter_values[, i]
#     findInterval(x = x, seq(min(x), max(x), length.out = num_intervals[[i]]+1), all.inside = TRUE)
#   })
# 
#   ## Translate + Project values onto base level set
#   proj_fv <- t((t(filter_values) - filter_min) - t(original_level_sets - 1)*base_interval_length)
# 
#   ## Compute the lower and upper halfspace distances
#   halfspace_dist <- lapply(1:filter_dim, function(d_i){
#     res <- cbind(proj_fv[, d_i], t(base_interval_length[d_i] - t(proj_fv[, d_i])))
#     colnames(res) <- c("lower", "upper")
#     res
#   })
# 
#   ## Return the resulting distances
#   list(base_interval = base_interval_length, base_lsmi = original_level_sets, halfspace_dist = halfspace_dist)
# })
# 
# ## Computes the target level set flat indices
# CoverRef$methods(computeTargetLSFI = function(){
#   if ("uninitializedField" %in% class(num_intervals)){ stop("'num_intervals' must be set prior to computing all overlap values.") }
# 
#   ## Get filter settings
#   filter_dim <- ncol(filter_values) ## filter dimensionality
#   filter_rng <- apply(filter_values, 2, range)
#   { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
#   filter_len <- filter_max - filter_min
# 
#   ## Base interval length (step size == 0), i.e. interval length corresponding to
#   ## the cover that is zero-dimensional w.r.t. its Lebesgue covering dimension
#   base_interval_length <- (filter_len)/(num_intervals) ## base interval length
# 
#   ## Find which level set each point lies within under the case the given filter space is
#   ## Note that is equivalent to integer division of each dimension by the base interval length given above,
#   ## but 'findInterval' is a bit safer for handling edge cases
#   original_level_sets <- sapply(1:filter_dim, function(i) {
#     x <- filter_values[, i]
#     findInterval(x = x, seq(min(x), max(x), length.out = num_intervals[[i]]+1), all.inside = TRUE)
#   })
# 
#   ## For each dimension of Z, retrieve the box sizes required to intersect the k^th nearest
#   ## level sets, for all k_i \in {1, 2, ..., < number of intervals >}
#   res <- vector(mode = "list", length=filter_dim)
#   orig_lsmi <- original_level_sets
#   for (d_i in 1:filter_dim){
#     k_di <- 1:num_intervals[d_i]
#     target_delta <- t(sapply(orig_lsmi[,d_i], function(orig_idx) { k_di[which(k_di != orig_idx)] - orig_idx }))
#     target_lsfi <- t(sapply(1:nrow(target_delta), function(i){
#       right_dim <- if(d_i == filter_dim){  NULL } else { (d_i+1):filter_dim }
#       target_lsmi <- cbind(orig_lsmi[i, 0:(d_i-1)], orig_lsmi[i, d_i] + target_delta[i,], orig_lsmi[i, right_dim])
#       getLSFI(target_lsmi)
#     }))
#     res[[d_i]] <- list(target_lsfi=target_lsfi, delta = target_delta)
#   }
# 
#   ## Return resulting level set flat indices
#   return(res)
# })

# base - 3
# target_dir <- t(sapply(orig_lsmi[,d_i], function(orig_idx){
#   tmp <- target_delta
#   tmp[1:orig_idx] <- -tmp[1:orig_idx]
#   tmp <- Filter(f = function(lsfi) abs(lsfi) != orig_idx, tmp)
#   ifelse(tmp < 0, -(orig_idx + tmp), tmp - orig_idx)
# }))
# target_lsfi <- t(sapply(1:nrow(target_dir), function(i){
#   getLSFI(sapply(1:filter_dim, function(d_ii) {
#     if (d_ii == d_i) { orig_lsmi[i, d_i] + target_dir[i,] }
#     else { rep(orig_lsmi[i, d_ii], filter_dim) }
#   }))
# }))
  # browser()

  # ## Get the per-point interval length which would cause a level-set membership change
  # if (max_overlap <= 0.50){
  #   ## Get the distance to the nearest box
  #   intersect_dist <- sapply(1:filter_dim, function(i){
  #     R_tmp <- cbind(proj_fv[, i], (base_interval_length - proj_fv)[, i])
  #     apply(R_tmp, 1, min)
  #   })
  #
  #   ## The corresponding overlap occurs when each box expands to at least this distance
  #   box_size <- base_interval_length + intersect_dist
  #
  #   ## Final overlap percentages (per dimension), for only orthogonal level sets
  #   overlap_params <- (filter_len - box_size * num_intervals)/(box_size * (-num_intervals + 1))
  #   overlap_params <- round(overlap_params, digits = 14)
  #
  #   ## Retrieve the 'target' index sets
  #   ## This becomes easier with some magical/unreadable NSE used to expand the index set in every dimension
  #   index_diff <- t(apply((base_interval_length - proj_fv) < (box_size/2), 1, function(diff) ifelse(diff, 1L, -1L)))
  #   target_sets <- sapply(1:filter_dim, function(d_i){
  #     gamma <- eval(as.call(append(quote(cbind), c(rep(0, length(1:d_i) - 1L), quote(index_diff[, d_i]), rep(0, length(d_i:filter_dim) - 1L)))))
  #     getLSFI(original_level_sets+gamma)
  #   })
  #   ls_idx <- cbind(orig=getLSFI(original_level_sets), target_sets)
  #
  #   ## Return relevant info, including:
  #   ## overlap := the overlap percentage required to observe a level set intersection, per point and per dimension
  #   ## level_set_intersects := the LSFI's of the level sets intersected with the computed overlap values
  #   ## box_sizes := the minimum size of the boxes needed to observe the intersections
  #   return(list(overlap = overlap_params, level_set_intersects = ls_idx, box_sizes = box_size))
  # } else {

    # halfspace_direction <- sapply(1:filter_dim, function(d_i){
    #   ifelse(apply(halfspace_dist[[d_i]], 1, which.min) == 1, -1, 1)
    # })


    ## For each dimension of Z, retrieve the box sizes required to intersect the k^th nearest
    ## level sets, for all k_i \in {1, 2, ..., < number of intervals >}
    # overlap_info <- vector(mode = "list", length=filter_dim)
    # orig_lsmi <- original_level_sets
    # for (d_i in 1:filter_dim){
    #   target_delta <- 1:num_intervals[d_i]
    #   # pt_idx <- match(orig_lsfi[,d_i], target_lsfi)
    #   target_dir <- t(sapply(orig_lsmi[,d_i], function(orig_idx){
    #     tmp <- target_delta
    #     tmp[1:orig_idx] <- -tmp[1:orig_idx]
    #     tmp <- Filter(f = function(lsfi) abs(lsfi) != orig_idx, tmp)
    #     ifelse(tmp < 0, -(orig_idx + tmp), tmp - orig_idx)
    #   }))
    #   target_lsfi <- t(sapply(1:nrow(target_dir), function(i){
    #     getLSFI(sapply(1:filter_dim, function(d_ii) {
    #       if (d_ii == d_i) { orig_lsmi[i, d_i] + target_dir[i,] }
    #       else { rep(orig_lsmi[i, d_ii], length(target_dir[i,])) }
    #     }))
    #   }))
    #   box_sizes <- t(sapply(1:nrow(target_dir), function(i){
    #     p <- target_dir[i,]
    #     # browser()
    #
    #     ## If the original point lies on the border, then just add the halfspace distances to
    #     ## integer multiples of the base interval length; else, the "box" the point resides in
    #     ## must be have its halfspace distance doubled to account for growing in both directions.
    #     is_border <- sapply(orig_lsmi[,d_i], function(lsfi_di) if (lsfi_di %in% c(head(num_intervals, 1L), tail(num_intervals, 1L))) 1 else 2)
    #     lower <- halfspace_dist[[d_i]][i,1]*is_border
    #     upper <- halfspace_dist[[d_i]][i,2]*is_border
    #     ifelse(p < 0, abs(p*base_interval_length) + lower,
    #                   abs(p*base_interval_length) + upper)
    #   }))
    #   overlap <- (filter_len - box_sizes*num_intervals[d_i])/(box_sizes*(-num_intervals[d_i] + 1))
    #   overlap_info[[d_i]] <- list(overlap=overlap, target_lsfi=target_lsfi, box_sizes = box_sizes)
    # }

    ## The final overlap information is a list of d components, where each components stores:
    ## 1) overlap := an n x k matrix of overlap values specifying the overlap needed, in that dimension, to intersect another level set
    ## 2) target_lsfi := an n x k matrix of the level set flat indices that the point at index i will intersect with
    ## 3) box_sizes := an n x k matrix of the interval lengths (box sizes) needed to for the point at index i to intersect another level set
    # return(overlap_info)
#
#
#     target_lsfi <- Filter(f = function(lsfi) abs(lsfi) != original_level_sets[1,1], target_lsfi)
#
#     box_sizes <- sapply(target_lsfi, function(p) {
#       ifelse(p < 0, abs(p*base_interval_length) + halfspace_dist[[1]][1,1],
#                     abs(p*base_interval_length) + halfspace_dist[[1]][1,2])
#     })
#     # target_lsfi[pt_idx:num_intervals[1]]-seq(1L, original_level_sets[1,1], length.out = original_level_sets[1,1]-1L),
#     #                  seq(original_level_sets[1,1], num_intervals[1], length.out = num_intervals[1]-original_level_sets[1,1]))
#     # target_lsfi[which(target_lsfi )]
#     ## Given a point LSFI, compute
#     ## Given the distances to the lower and upper halfspaces, compute the scalar multiple
#     for (d_i in 1:filter_dim){
#
#     }
#
#
#     ## Given a points directions (-1 or +1) and its intersection distance, compute the box distances
#     ## and LSFI's for all orthogonal directions
#     function(pt_direction, int_dist){
#       sapply(1:filter_dim, function(d_i){
#
#       })
#     }
#     max_box_size <- filter_len/(num_intervals - c(0.5, 0.5)*(num_intervals - 1))
#     ## Halfspace calculation
#     intersect_dist <- list(near=round(pmin(proj_fv, t(base_interval_length - t(proj_fv))), digits = 14),
#                            far=round(pmax(proj_fv, t(base_interval_length - t(proj_fv))), digits = 14))
#     ## Assert near halfspace computation: all(t(intersect_dist$near) - (base_interval_length/2) < 0)
#     ## Assert far halfspace computation: all(t(intersect_dist$far) - (base_interval_length/2) >= 0)
#
#
#     ## Box size calculation
#     near_box_sizes <-  t(t(intersect_dist$near) + base_interval_length)
#     far_box_sizes <- t(t(intersect_dist$far) + base_interval_length)
#     ## Assert: all(far_box_sizes[, 2] > base_interval_length[2]*3/2)
#
#     ## Overlap calculation
#     num_intervals
#     num <- (filter_len - (t(near_box_sizes)*num_intervals))
#     denom <- (-t(near_box_sizes)*num_intervals) + t(near_box_sizes)
#     apply(t(num/denom), 2, range)
#
#     num <- (filter_len - (t(far_box_sizes)*num_intervals))
#     denom <- (-t(far_box_sizes)*num_intervals) + t(far_box_sizes)
#     far_overlap <- t(num/denom)
#     apply(far_overlap, 2, range)
#
#     # crit_box <- filter_len/(num_intervals - c(0.50, 0.50)*(num_intervals - 1L))
#     near_overlap <- t((filter_len - (t(near_box_sizes)*num_intervals))/(t(near_box_sizes)*(-num_intervals + 1)))
#     far_overlap <-  t((filter_len - (t(far_box_sizes)*num_intervals))/(t(far_box_sizes)*(-num_intervals + 1)))
  ## "Expansion Factor" interval length
  #int_len_exp_factor <- base_interval_length + artificial_int_len
# })


## Function to plot the filter space, including the level set bounds
# CoverRef$methods(plotFilterSpace = function(show_ls_bounds = TRUE, show_lsfi = FALSE, ...){
#   filter_dim <- ncol(filter_values)
#   cpar <- par("mar")
#   par(mar=rep(0.25, 4)) ## adjust margins prior to plotting
#   min_idx <- sum(head(index_set, 1))
#   max_idx <- sum(tail(index_set, 1))
#   rbw_pal <- rev(rainbow(100, start = 0, end = 4/6))
#   binned_color <- substr(rbw_pal[cut(apply(index_set, 1, sum), breaks = 100, labels = FALSE)], start = 0, stop = 7)
#   if (filter_dim > 2){ stop("Cannot plot a filter space with more than 3 dimensions") }
#   if (filter_dim == 1){
#     filter_sz <- nrow(filter_values)
#     alpha_scale <- ifelse(filter_sz >= 1000, 0.01, 1 - (log(filter_sz)-1)/(log(1000)-1))
#     plot(cbind(filter_values, rep(0L, length(filter_values))), ylim = c(-1, 1), pch = "|", cex = 5,
#          col = adjustcolor("black", alpha.f = alpha_scale),
#          ylab = "", xlab = "Filter Values")
#     i <- 1L
#     for (ls in level_sets){
#       lines(x = c(ls$bounds[1, 1], ls$bounds[2, 1]), y = c(-0.5, -0.5), col = binned_color[i], ...)
#       lines(x = c(ls$bounds[2, 1], ls$bounds[2, 1]), y = c(-0.5, 0.5), col = binned_color[i], ...)
#       lines(x = c(ls$bounds[2, 1], ls$bounds[1, 1]), y = c(0.5, 0.5), col = binned_color[i], ...)
#       lines(x = c(ls$bounds[1, 1], ls$bounds[1, 1]), y = c(0.5, -0.5), col = binned_color[i], ...)
#       i <- i + 1
#     }
#   } else if (filter_dim == 2){
#     min_bounds <- head(level_sets, 1)[[1]]$bounds[1, ]
#     max_bounds <- tail(level_sets, 1)[[1]]$bounds[2, ]
#     plot(filter_values, xlim = c(min_bounds[1], max_bounds[1]), ylim = c(min_bounds[2], max_bounds[2]), axes = FALSE)
#     i <- 1L
#     for (ls in level_sets){
#       lines(x = c(ls$bounds[1, 1], ls$bounds[2, 1]), y = c(ls$bounds[1, 2], ls$bounds[1, 2]), col = binned_color[i], ...)
#       lines(x = c(ls$bounds[2, 1], ls$bounds[2, 1]), y = c(ls$bounds[1, 2], ls$bounds[2, 2]), col = binned_color[i], ...)
#       lines(x = c(ls$bounds[2, 1], ls$bounds[1, 1]), y = c(ls$bounds[2, 2], ls$bounds[2, 2]), col = binned_color[i], ...)
#       lines(x = c(ls$bounds[1, 1], ls$bounds[1, 1]), y = c(ls$bounds[2, 2], ls$bounds[1, 2]), col = binned_color[i], ...)
#       if (show_lsfi){
#         text(x = mean(ls$bounds[, 1]), y = mean(ls$bounds[, 2]), labels = i)
#       }
#       i <- i + 1
#     }
#   }
#   par(mar = cpar) ## reset
# })

