#' @title Multiscale Mapper
#' @description Computes all unique Mapper constructions. Only supports rectangular-like covers.
#' @examples 
#'   data("noisy_circle")
#'   ll <- noisy_circle[which.min(apply(noisy_circle, 1, sum)),]
#'   f_x <- t(apply(noisy_circle, 1, function(x) abs(x - ll)))
#'   m_ref <- Mapper::mapper(noisy_circle, filter_values = f_x)
#' @export
MapperRef$set("public", "multiscale", 
  function(mapper_obj, type = c("restrained rectangular"), ...){
    if (class(mapper_obj) != "MapperRef"){ stop("'multiscale' expects a Mapper reference object.") }
    if (mapper_obj$cover$type != "rectangular"){ stop("'multiscale' is only compatible with rectangular covers.") }
    n <- ifelse("dist" %in% class(X), attr(X, "Size"), nrow(X))
    d <- ncol(mapper_obj$cover$filter_values)
  
  }
)
#   ## Start with the base cover
#   mapper_obj$cover$setOverlap(0.0)
#   mapper_obj$cover$constructCover()
# 
#   ## Get the range of the indices per dimension
#   idx_rng <- apply(mapper_obj$cover$index_set, 2, range)
# 
#   ## Get the LSFI's for each point
#   pt_lsfi <- vector(mode = "numeric", length = n)
#   for (i in 1:length(mapper_obj$cover$level_sets)){
#     ls <- mapper_obj$cover$level_sets[[i]]
#     pt_lsfi[ls$points_in_level_set] <- i
#   }
# 
#   ## Translate to LSMI
#   pt_lsmi <- as.matrix(mapper_obj$cover$getLSMI(pt_lsfi))
# 
#   for (d_i in 1:d){
# 
#   }
# 
#     # overlap_params <- seq(0.1, 0.6, by = 0.1)
#     # overlap_params <- cbind(overlap_params, overlap_params)
#     ## Order the overlap parameters in non-decreasing order; save the ordering for return
#     # input_order <- order(overlap_params)
#     n_pts <- nrow(m$X)
#     n_lvl_sets <- length(m$cover$level_sets)
#     n_blocks <- nrow(overlap_params) + 1
#     browser()
# 
#     ## First step: construct a base cover (where overlap == 0)
#     m$cover$setOverlap(0)
#     m$cover$setupCover()
#   #   # 
#   #   # ## Get the halfspace distances
#   #   # hf_dist <- m$cover$computeHalfSpaceDistances()
#   #   # 
#   #   # ## Get which points lie in a boundary set
#   #   # boundary_pts <- apply(sapply(1:ncol(m$cover$filter_values), function(d_i) {
#   #   #   hf_dist$base_lsmi[, d_i] %in% c(1, m$cover$num_intervals[d_i])
#   #   # }), 1, any)
#   #   # 
#   #   # ## For each new overlap parameter, update:
#   #   # ## 1) The points within the level set
#   #   # ## 2) The pairwise distances of the points (in the data space) from (1)
#   #   # ## 3) The nodes and edges in the graph
#   #   # target_info <- m$cover$computeTargetLSFI()
#   #   # 
#   #   # ## Helper function to determine if a given level set lies on the boundary
#   #   # isBoundaryLS <- function(lsfi, specific_dim = NULL){
#   #   #   ls_lsmi <- m$cover$getLSMI(lsfi)
#   #   #   if (missing(specific_dim) || is.null(specific_dim)){
#   #   #     res <- apply(ls_lsmi, 1, function(lsmi) { any(sapply(1:filter_dim, function(d_i){ lsmi[d_i] %in% c(1, m$cover$num_intervals[d_i]) })) })
#   #   #   } else {
#   #   #     res <- apply(ls_lsmi, 1, function(lsmi) { lsmi[specific_dim] %in% c(1, m$cover$num_intervals[specific_dim]) })
#   #   #   }
#   #   #   return(res)
#   #   # }
#   #   # 
#   #   # ## Compute the overlap percentage necessary to merge with each of the target level sets
#   #   # ## This effectively, per point, the overlap needed for the target level set to intersect the current point
#   #   # target_box_sizes <- lapply(1:filter_dim, function(d_i){
#   #   #   t(sapply(1:n_pts, function(i){
#   #   #     p <- target_info[[d_i]]$delta[i,] ## offset vector
#   #   #     # boundary_constant <- m$cover$num_intervals[d_i]-1 # ifelse(isBoundaryLS(target_info[[d_i]]$target_lsfi[1,], specific_dim = d_i), 1, m$cover$num_intervals[d_i]-1) ## internal boxes expand k-times as fast
#   #   #     base_box_sizes <- abs(p)*hf_dist$base_interval[[d_i]]
#   #   #     hf_dir_dist <- ifelse(p < 0, hf_dist$halfspace_dist[[d_i]][i, 1], hf_dist$halfspace_dist[[d_i]][i, 2])
#   #   #     boundary_constant <- num_intervals[d_i] - abs(p)
#   #   #     box_sizes <- base_box_sizes + hf_dir_dist*boundary_constant ## Multiply by boundary constant
#   #   #     box_sizes
#   #   #   }))
#   #   # })
#   #   # 
#   #   # all_box_sizes <- cbind(target_box_sizes[[1]][10, ], target_box_sizes[[2]][10, ])
#   #   # wut <- m$cover$computeOverlap(all_box_sizes)
#   #   # m$cover$setOverlap(wut[2,2])
#   #   # m$cover$setupCover()
#   #   # m$cover$plotFilterSpace(show_ls_bounds = TRUE, show_lsfi = TRUE)
#   #   # points(m$cover$filter_values[10, 1], m$cover$filter_values[10, 2], col = "red") ## needs to be near 0.955
#   #   # 
#   #   # diff(m$cover$level_sets[[25]]$bounds[, 2])
#   #   # 
#   #   # ## Say we want to compute the distances to the nearest box
#   #   # nearer_hf_dist <- cbind(apply(hf_dist$halfspace_dist[[1]], 1, min), apply(hf_dist$halfspace_dist[[2]], 1, min))
#   #   # nearest_box_sizes <- lower_hf_dist + ifelse(boundary_pts, hf_dist$base_interval, hf_dist$base_interval/2)
#   #   # 
#   #   # m$cover$setOverlap(0)
#   #   # m$cover$setupCover()
#   #   # m$cover$plotFilterSpace(show_ls_bounds = TRUE, show_lsfi = TRUE)
#   #   # points(m$cover$filter_values[10, 1], m$cover$filter_values[10, 2], col = "red")
#   #   # 
#   #   # lines(x = c(m$cover$filter_values[1, 1], m$cover$filter_values[1, 1]+nearer_hf_dist[1,1]),
#   #   #      y = c(m$cover$filter_values[1, 2], m$cover$filter_values[1, 2]),
#   #   #      col = "red")
#   #   # 
#   #   # test_box_size <- matrix(c(hf_dist$base_interval[1] + nearer_hf_dist[1,1]*4, hf_dist$base_interval[2]), ncol = 2)
#   #   # wut <- m$cover$computeOverlap(box_sizes = test_box_size)
#   #   # m$cover$setOverlap(wut)
#   #   # m$cover$setupCover()
#   #   # m$cover$plotFilterSpace(show_ls_bounds = TRUE, show_lsfi = TRUE)
#   #   # points(m$cover$filter_values[1, 1], m$cover$filter_values[1, 2], col = "red")
#   #   # abline(v = m$cover$level_sets[[2]]$bounds[1, 1], col = "red")
#   #   # abline(v = m$cover$level_sets[[2]]$bounds[2, 1], col = "red")
#   #   # 
#   #   # ## Get the base level set flat indices
#   #   # base_lsfi <- m$cover$getLSFI(hf_dist$base_lsmi)
#   #   # g <- m$cover$computeOverlap()           ## Get the set of possible overlap values
#   #   # # overlap_order <- order(g$overlap)       ## Order them in increasing order
#   #   # # g <- lapply(g, function(g_comp) matrix(g_comp[overlap_order,])) ## Reorder all components of g
#   #   # 
#   #   # ## Detect which points are going to intersect each level set, indexed by the level set flat index (lsfi)
#   #   # max_overlap <- apply(overlap_params, 2, max)
#   #   # filter_dim <- ncol(m$cover$filter_values)
#   #   # 
#   #   # ## Sort G in increasing order
#   #   # for (d_i in 1:filter_dim){
#   #   #   overlap_order <- apply(g[[d_i]]$overlap, 1, order)
#   #   #   g[[d_i]]$target_lsfi <- t(sapply(1:n_pts, function(i){
#   #   #     g[[d_i]]$target_lsfi[i, overlap_order[, i]]
#   #   #   }))
#   #   #   g[[d_i]]$overlap <- t(apply(g[[d_i]]$overlap, 1, sort))
#   #   # }
#   #   # 
#   #   # 
#   #   # ## TODO: wrap this in C++
#   #   # ## Create  a list of point indices, per level set, which will intersect the given level set at some overlap value
#   #   # ## les than the current maximum overlap requested
#   #   # new_ls_pt_intersects <- vector(mode = "list", length = n_lvl_sets)
#   #   # for (d_i in 1:filter_dim){
#   #   #   for (i in 1:n_pts){
#   #   #     targets <- g[[d_i]]$target_lsfi[i, (g[[d_i]]$overlap[i,] < max_overlap[d_i]) ] ## only consider points below the max overlap
#   #   #     for (lsfi in targets){
#   #   #       new_ls_pt_intersects[[lsfi]] <- c(new_ls_pt_intersects[[lsfi]], i)
#   #   #     }
#   #   #   }
#   #   # }
#   #   # 
#   #   # ## For each level set index, get the points that *will* intersect that level set
#   #   # ## before or at the maximum overlap value needed.
#   #   # # new_ls_pt_intersects <- lapply(1:length(m$cover$level_sets), function(lsfi) {
#   #   # #   as.vector(unlist(sapply(1:filter_dim, function(d_i) {
#   #   # #     sapply(1:nrow(g[[d_i]]$overlap, function(i){
#   #   # #       any(g[[d_i]]$target_lsfi[i,] == lsfi)
#   #   # #     })
#   #   # #   apply(g[[1]]$target_lsfi, 1, function(lsfis) any(lsfis == lsfi))
#   #   # #
#   #   # #     which(g$ls_idx[, d_i + 1] == lsfi & g$overlap[, d_i] <= max_overlap[d_i])
#   #   # #   })))
#   #   # # })
#   #   # 
#   #   # # new_ls_pt_intersects <- lapply(1:length(m$cover$level_sets), function(lsfi) {
#   #   # #   ## For each level set index, get the points that *will* intersect that level set
#   #   # #   ## before or at the maximum overlap value needed.
#   #   # #   as.vector(unlist(sapply(1:filter_dim, function(d_i) {
#   #   # #     which(g$ls_idx[, d_i + 1] == lsfi & g$overlap[, d_i] <= max_overlap[d_i])
#   #   # #   })))
#   #   # # })
#   #   # # ls_intersections <- lapply(1:length(m$cover$level_sets), function(lsfi) overlap_order[which(g$ls_idx[, 2] == lsfi & g$overlap < max_overlap)]) ## which points
#   #   # # full_intersections <- mapply(function(ls, additional_pts) { union(ls$points_in_level_set, additional_pts) }, m$cover$level_sets, ls_intersections)
#   #   # 
#   #   # ## Pre-compute the level set distances under the *highest* overlap
#   #   # ## Although there may be an asymptotically lower runtime associated with appending distances by reference
#   #   # ## as points begin to intersect level set, computationally it's much more efficient to simply pre-compute the
#   #   # ## dissimilarities for the highest overlap needed, then subset by index
#   #   # ls_dists <- mapply(function(lvl_set, pt_idx){
#   #   #   max_pt_idx <- union(lvl_set$points_in_level_set, pt_idx)
#   #   #   if ("dist" %in% class(X)){ dist_subset(dist = X, idx = max_pt_idx) }
#   #   #   else { dist(X[max_pt_idx,], method = "euclidean") }
#   #   # }, m$cover$level_sets, new_ls_pt_intersects)
#   #   # 
#   #   # ## To get the Mapper construction per given parameter, first need to isolate which points and corresponding level sets
#   #   # ## need updating *per construction*. I call these groups of information which need to be incrementally updated 'blocks'
#   #   # # overlap_cuts <- sapply(1:filter_dim, function(d_i) findInterval(g$overlap[, d_i], vec = overlap_params[,d_i])) ## block by given overlap parameters
#   #   # block_ids <- sort(unique(as.vector(overlap_blocks)))
#   #   # overlap_cuts <- lapply(1:filter_dim, function(d_i){
#   #   # 
#   #   #   ## Find which block the computed overlap values correspond to (earliest)
#   #   #   t(apply(g[[d_i]]$overlap, 1, function(overlap) findInterval(overlap, vec = overlap_params[, d_i])))
#   #   # })
#   #   # 
#   #   #   # pts_to_update <- lapply(block_ids, function(b_i){
#   #   #   #   which(apply(overlap_blocks, 1, function(pt_mem){ b_i %in% pt_mem }))
#   #   #   # })
#   #   #   # names(pts_to_update) <- as.character(block_ids)
#   #   # 
#   #   #   ## For each block, extract:
#   #   #   ## 1) A named list, one for each level set, of point indices that are being unioned with that level set in this block update
#   #   #   update_tables <- lapply(1:filter_dim, function(d_i){
#   #   #     k_di <- ncol(g[[d_i]]$target_lsfi)
#   #   #     pt_idx <- t(sapply(1:n_pts, function(i){ rep(i, k_di) }))
#   #   #     data.table::data.table(pt_idx = as.vector(t(pt_idx)), lsfi = as.vector(t(g[[d_i]]$target_lsfi)), block_id = as.vector(t(overlap_blocks)))
#   #   #   })
#   #   # 
#   #   #   ## Find which level set each point lies within under the case the given filter space is
#   #   #   ## zero-dimensional w.r.t. its Lebesgue covering dimension
#   #   #   num_intervals <- m$cover$num_intervals
#   #   #   original_level_sets <- sapply(1:filter_dim, function(i) {
#   #   #     x <- m$cover$filter_values[, i]
#   #   #     findInterval(x = x, seq(min(x), max(x), length.out = num_intervals[[i]]+1), all.inside = TRUE)
#   #   #   })
#   #   #   orig_lsfi <- m$cover$getLSFI(original_level_sets)
#   #   # 
#   #   # 
#   #   #   to_update <- lapply(1:length(block_ids), function(i) vector(mode = "list", length = n_lvl_sets))
#   #   #   names(to_update) <- as.character(block_ids)
#   #   #   added_so_far <- vector(mode = "list", length = n_lvl_sets)
#   #   #   for (b_i in block_ids){
#   #   #     for (pt_id in 1:n_pts){
#   #   #       to_lsfi <- sapply(1:filter_dim, function(d_i){ update_tables[[d_i]][pt_idx == pt_id & block_id == b_i]$lsfi })
#   #   #       if (length(unlist(to_lsfi)) == 0){ next; }
#   #   #       to_lsfi <- matrix(to_lsfi, ncol = filter_dim)
#   #   #       update_lsfi <- lapply(1:nrow(to_lsfi), function(i){  m$cover$getRectLSFI(from = orig_lsfi[pt_id], to = to_lsfi[i,]) })
#   #   #       for (lsfi in setdiff(unique(unlist(update_lsfi)), orig_lsfi[pt_id])){
#   #   #         c_pts <- to_update[[as.character(b_i)]][[lsfi]]
#   #   #         pts_to_consider <- union(c(c_pts, pt_id)) ## Add the point to level set
#   #   #         new_pts <- setdiff(pts_to_consider, added_so_far[[lsfi]]) ## make sure it hasn't been added before
#   #   #         to_update[[as.character(b_i)]][[lsfi]] <- new_pts
#   #   #       }
#   #   #     }
#   #   #     for (lsfi in 1:n_lvl_sets){
#   #   #       update_pts <- union(added_so_far[[lsfi]], to_update[[as.character(b_i)]][[lsfi]])
#   #   #       if (!is.null(update_pts)){ added_so_far[[lsfi]] <- update_pts }
#   #   #     }
#   #   #   }
#   #   # 
#   #   #   lvl_set_of_interest <- 25
#   #   #   pt_colors <- rainbow(length(block_ids))
#   #   #   m$cover$plotFilterSpace(show_ls_bounds = TRUE, show_lsfi = TRUE)
#   #   #   for (b_i in 1:length(block_ids)){
#   #   #     new_pts <- to_update[[b_i]][[lvl_set_of_interest]]
#   #   #     if (!is.null(new_pts)){
#   #   #       points(m$cover$filter_values[new_pts, 1], m$cover$filter_values[new_pts, 2], col = pt_colors[b_i])
#   #   #     }
#   #   #     invisible(readline(prompt="Press [enter] to continue"))
#   #   #   }
#   #   # 
#   #   # 
#   #   #   to_update <- lapply(block_ids, function(b_i){
#   #   #     tmp <- ls_intersects_by_block[block_id == b_i]
#   #   #     uniq_lsfi <- unique(tmp$lsfi)
#   #   #     res <- lapply(uniq_lsfi, function(lsfi_i){ tmp[lsfi == lsfi_i]$pt_idx })
#   #   #     names(res) <- as.character(uniq_lsfi)
#   #   #     return(res)
#   #   #   })
#   #   #      # lvl_sets <- g[[d_i]]$target_lsfi[overlap_blocks == b_i] ## Which level sets need to be updated in this block (dimension-specific)
#   #   #     # pts_to_update <- which(apply(overlap_blocks == b_i, 1, any))
#   #   #     # print(length(lvl_sets))
#   #   #     # print(length(pts_to_update))
#   #   #     # data.table::data.table(lsfi = lvl_sets, pt_idx = pts_to_update)
#   #   #     # pts_in_lvl_sets <- lapply(lvl_sets, function(lsfi){
#   #   #     #   which(apply(overlap_blocks == lsfi, 1, any))
#   #   #     # })
#   #   #     # names(pts_in_lvl_sets) <- lvl_sets
#   #   #     # Filter(function(pt_membership) length(pt_membership) > 0, pts_in_lvl_sets)
#   #   # 
#   #   # 
#   #   # 
#   #   # #   orig_lsfi <- m$cover$getLSFI(original_level_sets)
#   #   # #   sapply(1:length(orig_lsfi), function(i){
#   #   # #     targets <- cbind(g[[1]]$target_lsfi[i,], g[[2]]$target_lsfi[i,])
#   #   # #     apply(targets, 1, function(pair) m$cover$getRectLSFI(orig_lsfi[i], pair))
#   #   # #   })
#   #   # #   wut <- cbind(g[[1]]$target_lsfi[1,], g[[2]]$target_lsfi[1,])
#   #   # #
#   #   # #   # m$cover$getRectLSFI(orig_lsfi, wut[1, ])
#   #   # #
#   #   # #   to_update
#   #   # # })
#   #   # 
#   #   # 
#   #   # block_ids <- (1:n_blocks) - 1L
#   #   # sapply(1:filter_dim, function(d_i){ apply(g[[d_i]]$overlap, 1, function(overlap) findInterval(overlap, vec = overlap_params[, d_i])) })
#   #   # lapply(block_ids, function(){
#   #   #   sapply(1:filter_dim, function(d_i){
#   #   #     apply(g[[d_i]]$overlap, 1, function(overlap) findInterval(overlap, vec = overlap_params[, d_i]))
#   #   #   })
#   #   # })
#   #   # 
#   #   # ## Remove overlaps corresponding to invalid target level sets
#   #   # overlap_cuts[which(is.na(g$ls_idx[, -1]))] <- NA
#   #   # block_ids <- sort(unique(as.integer(overlap_cuts)))
#   #   # 
#   #   # ## Utility function to determine which points need to be updated for a given block id
#   #   # need_updated <- function(block_id){
#   #   #   which(apply(overlap_cuts, 1, function(active_blocks) any(block_id %in% active_blocks)))
#   #   # }
#   #   # 
#   #   # ## Aggregate all of the information needed to incrementally construct Mapper constructions. This includes
#   #   # ## 1) A named list for each level set containing the point ids which will intersect that level set for the given block
#   #   # ## 2) The upper bound on the overlap associated with this block
#   #   # to_update <- lapply(block_ids, function(block_id) {
#   #   #   point_ids <- need_updated(block_id) ## which points must be updated in this block?
#   #   #   level_set_subset <- g$ls_idx[point_ids,-1] ## which level sets will these new points intersect in this block?
#   #   #   levels_to_update <- unique(as.vector(na.omit(as.vector(level_set_subset)))) ## what level sets (lsfi's) need to be updated this block?
#   #   #   pts_to_update_by_ls <- lapply(levels_to_update, function(lsfi){
#   #   #     as.vector(unlist(apply(level_set_subset, 2, function(target_ls_indexes) {
#   #   #       point_ids[which(target_ls_indexes == lsfi)]
#   #   #     })))
#   #   #   })
#   #   #   names(pts_to_update_by_ls) <- levels_to_update
#   #   #   list(points_to_update = pts_to_update_by_ls,
#   #   #        overlap_upper_bound = overlap_params[block_id+1,])
#   #   # })
#   #   # 
#   #   # ## Compute the new set of mapper instances
#   #   # all_mappers <- vector(mode = "list", length = length(to_update))
#   #   # i <- 0L
#   #   # for (block in to_update){
#   #   #   for (level in names(block$points_to_update)){
#   #   #     c_pts <- m$cover$level_sets[[as.integer(level)]]$points_in_level_set
#   #   # 
#   #   #     ## Update which points lie in the level set
#   #   #     new_pts <- union(c_pts, block$points_to_update[[level]])
#   #   #     m$cover$level_sets[[as.integer(level)]]$points_in_level_set <- new_pts
#   #   # 
#   #   #     ## Update the dist object in the level set with the (relatively indexed) newly intersected points
#   #   #     m$cover$level_sets[[as.integer(level)]]$dist <- TDAmapper:::dist_subset(dist = ls_dists[[as.integer(level)]],
#   #   #                                                                             idx = which(new_pts %in% new_ls_pt_intersects[[as.integer(level)]]))
#   #   #   }
#   #   #   ## Update nodes (Since dissimilarities have been updated, just rerun clustering)
#   #   #   m$computeNodes(... =  block$levels_to_update)
#   #   # 
#   #   #   ## Update edges (but only the by comparing the level sets which need updating)
#   #   #   levels_to_update <- TDAmapper:::valid_pairs(g$ls_idx[unname(unlist(block$points_to_update)),])
#   #   #   m$computeEdges(level_sets = block$levels_to_update)
#   #   # 
#   #   #   ## Story copy of object
#   #   #   all_mappers[[i]] <- m$G # m$copy(shallow = FALSE)
#   #   #   i <- i + 1
#   #   # }
#   # 
#   # })
#   
#  
# 
# })
# 
# # Returns whether a lsmi is on the border
# is_border <- function(idx, index_set){ return(idx %in% range(index_set)) }
# 
# # Compute the interval size(s) needed for a given point p to intersect its neighboring level sets
# # Accepts as inputs:
# # p := point in question
# # k := number of intervals
# # z := range of the filter domain
# # p_lsmi := the level set multi index (lsmi) of point p
# # index_set := the index set of the filter space spanned by z
# computeR <- function(p, k, z, p_lsmi, index_set){
#   lower_lsfi <- index_set[index_set < p_lsmi]
#   upper_lsfi <- index_set[index_set > p_lsmi]
#   lower_is_border <- is_border(lower_lsfi, index_set)
#   upper_is_border <- is_border(upper_lsfi, index_set)
#   lower_R <- sapply(1:length(lower_is_border), function(i){
#     is_edge_case <- lower_is_border[i]
#     base_R <- (z / k)
#     if (is_edge_case){ return(base_R + abs(base_R - p)) }
#     else { return((-k * p + p + z)/(-k*lower_lsfi[i] + lower_lsfi[i] + 1)) }
#   })
#   upper_R <- sapply(1:length(upper_is_border), function(i){
#     is_edge_case <- upper_is_border[i]
#     base_R <- (z / k)
#     if (is_edge_case){ return(base_R + abs(base_R - p)) }
#     else { return(z - ((k - 1)/(upper_lsfi[i]))*p) }
#   })
#   c(lower_R, upper_R)
# }
# 
# 
