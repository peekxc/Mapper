#' Reference Class (R5) implementation of Mapper for the Mapper package in R
#' Composes a set of classes to perform Mapper efficiently.
#' Author: Matt Piekenbrock

## Create the reference class, along with required field types
#' @export
mapper_ref <- setRefClass("MapperRef", fields = list(X="ANY", cover="ANY", clustering_algorithm="ANY", G="ANY", config="list"))

## The cover stores the filter values
mapper_ref$methods(setCover = function(fv, type = c("restrained rectangular"), ...){
  if (is.null(dim(fv)) && is.numeric(fv)){ fv <- matrix(fv, nrow = length(fv), ncol = 1) }
  if (is.null(dim(fv)) || (!"matrix" %in% class(fv))) { stop("Filter values must be numeric and in matrix form") }
  if (missing(type) || type == "restrained rectangular"){
    cover <<- r_rect_cover$new(filter_values = fv, type = "restrained rectangular", ...) # additional arguments to pass to the rectangular cover initialization method
  } else if (type == "fixed rectangular"){
    cover <<- f_rect_cover$new(filter_values = fv, type = "fixed rectangular", ...) # additional arguments to pass to the rectangular cover initialization method
  } else {
    stop(sprintf("%s cover not supported. Perhaps consider contributing it to the package w/ a pull request?", type))
  }
})

## Changes the clustering algorithm used by mapper.
## Must accept a 'dist' object and return a static or 'flat' clustering result
mapper_ref$methods(setClusteringAlgorithm = function(cl = c("single", "ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid"), num_bins = 10L, ...){

  ## Use a linakge criterion + cutting rule
  if (missing(cl)) { cl <- "single" }
  if (class(cl) == "character"){
    hclust_opts <- c("single", "ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
    if (!cl %in% hclust_opts){ stop(sprint("Unknown linkage method passed. Please use one of (%s). See ?hclust for details.", paste0(hclust_opts, collapse = ", "))) }
    clustering_algorithm <<- function(dist_x, num_bins = num_bins, ...){
      if (is.null(dist_x) || attr(dist_x, "Size") <= 1){ return(1L) } ## Handle cases where there are <= 1 points
      hcl <- hclust(dist_x, method = cl)
      cutoff_first_bin(hcl, num_bins)
    }
  }
  ## Allow the user to use their own clustering algorithm if they desire. Must accept a dist object and return a flat clustering result.
  else if (class(cl) == "function"){
    clustering_algorithm <<- cl
  } else {
    stop("'cl' must be either a linkage criterion (character) or a clustering algorithm (function)")
  }
})

## Computes the distance matrix for each level set
mapper_ref$methods(computeLevelSetDist = function(measure = "euclidean", ...){
  for (i in 1:length(cover$level_sets)){
    pt_idx <- cover$level_sets[[i]]$points_in_level_set
    if (length(pt_idx) <= 1){ cover$level_sets[[i]]$dist <<- dist(0L) }
    else {
      if ("dist" %in% class(X)){
        cover$level_sets[[i]]$dist <<- dist_subset(dist = X, idx = pt_idx)
      } else {
        if (is.character(measure)){
          cover$level_sets[[i]]$dist <<- parallelDist::parallelDist(X[pt_idx,], method = measure)
        } else if (is.function(measure)){
          cover$level_sets[[i]]$dist <<- measure(X[pt_idx,], ...)
        } else {
          stop("Unknown 'measure' argument. Must be either character string or function.")
        }
      }
    }
  }
})

## Computes the K-skeleton in a similar fashion as described by Singh et. al, section 3.2. Note that this
## procedure implicitly assumes the full abstract simplicial complex produced by Mapper can be treated as
## a flag complex.
mapper_ref$methods(computeKSimplex = function(k, ...){
  if (k < 1L){ stop("'computeKSimplex' expects a positive integer representing the highest dimension simplex to compute.") }
  if (k >= 1L){ computeNodes(...) }
  if (k >= 2L){ computeEdges(...) }
  if (k >= 3L){
    k_simplex <- vector(mode = "list", length = k - 2L)
    for (k_i in 3L:k){
      pt_map <- lapply(1:length(m$nodes), function(lsfi) cbind(pt_id = m$nodes[[lsfi]], node_id = rep(lsfi)))
      pt_map <- data.table::data.table(do.call(rbind, pt_map))
      intersection_pts <- pt_map[, .(n_intersects = nrow(.SD)), by = pt_id][n_intersects == k_i]$pt_id
      simplex_pts <- pt_map[pt_id %in% intersection_pts, .SD, by = pt_id]
      index <- (seq(nrow(simplex_pts))-1) %% k_i
      k_simplex_mat <- unique(do.call(cbind, lapply(0:(k_i - 1L), function(idx){ simplex_pts$node_id[index == idx] })))
      k_simplex[[k_i - 2L]] <- lapply(1:nrow(k_simplex_mat), function(i) k_simplex_mat[i,])
    }
  }
})

## The cover stores the filter values
mapper_ref$methods(computeNodes = function(...){
  if ("uninitializedField" %in% class(cover)){ stop("Cover must be constructed prior to computing nodes.") }

  ## Initialize the graph as an empty list
  G <<- list()

  ## Iterate through the 'dist' objects stored for each level set, perform the clustering
  cl_res <- lapply(cover$level_sets, function(ls) {
    if (length(ls$dist) > 0) clustering_algorithm(ls$dist, ...)
  })

  ## Precompute useful variables to know
  n_nodes <- sum(sapply(cl_res, function(cl) length(unique(cl))))
  node_idx <- unlist(mapply(function(cl, ls_i) if (length(cl) > 0) paste0(ls_i, ".", unique(cl)), cl_res, 1:length(cl_res)))
  n_lvlsets <- length(cover$level_sets)

  ## Agglomerate the nodes into a list. This matches up the original indexes of the filter values with the
  ## the clustering results, such that each node stores the original filter index values as well as the
  ## creating a correspondence between the node and it's corresponding level set flat index (lsfi)
  ## TODO: Cleanup and either vectorize or send down to C++
  G$nodes <<- vector(mode = "list", length = n_nodes)
  n_i <- 1L
  for (lsfi in 1:n_lvlsets){
    cl_i <- cl_res[[lsfi]]
    if (!is.null(cl_i)){
      ## Extract the node point indices for each cluster
      node_pt_idx <- lapply(unique(cl_i), function(cl_idx) cover$level_sets[[lsfi]]$points_in_level_set[which(cl_i == cl_idx)])
      for (node in node_pt_idx){
        attr(node, "level_set") <- lsfi
        if (any(is.na(node))){ browser() }
        G$nodes[[n_i]] <<- node
        n_i <- n_i + 1L
      }
    }
  }
})

mapper_ref$methods(computeEdges = function(level_sets = NULL){
  if ("uninitializedField" %in% class(G)){ stop("Graph with nodes must be constructed before computing edges.") }

  ## Retrieve the level set flat indices (LSFI) for each corresponding node
  node_lsfi <- sapply(G$nodes, function(node) attr(node, "level_set")) # which level set (by value) each node (by index) is in

  ## Create map from the level set flat index (by index) to the node indices the level set stores
  ## Note in this map empty level sets are NULL
  ls_node_map <- lapply(seq(length(cover$level_sets)), function(lvl_set_idx) {
    node_indices <- which(node_lsfi == lvl_set_idx)
    if (length(node_indices) == 0){ return(NULL) } else { return(node_indices) }
  })

  # browser()

  ## Retrieve the valid level set index pairs to compare. In the worst case, with no cover-specific optimization
  ## or 1-skeleton assumption, this may just be all pairwise combinations of LSFI's for the full simplicial complex.
  ## If the specific set of LSFI's were given
  if (missing(level_sets) || is.null(level_sets)){
    ## Let Rcpp handle the O(n^2) non-empty intersection checks
    ls_to_compare <- cover$valid_pairs()
    G$adjacency <<- adjacencyCpp(ls_pairs = ls_to_compare, nodes = G$nodes, ls_node_map = ls_node_map);
  } else {
    G$edgelist <<- Mapper:::edgeList_int(ls_pairs = level_sets, nodes = G$nodes, ls_node_map = ls_node_map)
  }
})

## Updates the mapper construction
mapper_ref$methods(update = function(){
  cover$constructCover()
  computeLevelSetDist()
  computeNodes(num_bins = 10L)
  computeEdges()
})


## plotNetwork uses the 'network' package to plot the Mapper construction, w/ suitable defaults
## corresponding to what's commonly used in practice, all of which can be overridden.
mapper_ref$methods(plotNetwork = function(...){

  ## Turn any given parameters into list
  params <- list(...)

  ## If vertex color supplied, great! If not, use rainbow scale from blue --> red with values
  ## corresponding to average filter values
  if (is.null(params[["vertex.col"]])){
    agg_pt_fv <- sapply(G$nodes, function(n_idx){ apply(matrix(cover$filter_values[n_idx,], nrow = length(n_idx)), 1, mean)})
    agg_node_fv <- sapply(agg_pt_fv, mean)
    rbw_pal <- rev(rainbow(100, start = 0, end = 4/6))
    binned_idx <- cut(agg_node_fv, breaks = 100, labels = F)
    params[["vertex.col"]] <- rbw_pal[binned_idx]
  }

  ## If vertex size specified, great! If not, scale node size logarithmically with the number of points each contains
  if (is.null(params[["vertex.cex"]])){
    node_sizes <- sapply(G$nodes, function(n_idx){ length(n_idx) })
    node_cex <- (6 - 1)*normalize(log(pmax(2, node_sizes)))+1 ## Scale node size logarithmically between [0.1, 4]
    params[["vertex.cex"]] <- node_cex
  }

  ## Other defaults
  if (is.null(params[["label"]])){ params[["label"]] <- 1:length(G$nodes) }
  if (is.null(params[["displaylabels"]])){ params[["displaylabels"]] <- TRUE }

  ## Construct the core network w/ network package
  base_net <- network::network.initialize(nrow(G$adjacency), directed = F)
  mapper_graph <- network::network.adjacency(x = G$adjacency, g = base_net)
  params[["x"]] <- mapper_graph

  ## Evaluate in terms of passed in parameters
  do.call(network::plot.network, params)
})

## S3-like print override
mapper_ref$methods(show = function(){
  if ("dist" %in% class(X)){ n <- attr(X, "Size") } else { n <- nrow(X) }
  if (!is.null(config$call))
    cat("\nCall:\n", deparse(config$call), "\n\n", sep = "")
  message <- sprintf("Mapper construction for %d objects", n)
  if (!"uninitializedField" %in% class(cover)){ message <- append(message, cover$summary()) }
  # if (!"uninitializedField" %in% class(clustering_algorithm)){ message <- append(message, cover$summary()) }
  writeLines(message)
})

only_combinations <- function(mat){
  mn <- pmin(mat[, 2], mat[, 1])
  mx <- pmax(mat[, 2], mat[, 1])
  int <- as.numeric(interaction(mn, mx))
  mat[match(unique(int), int),]
}

## Update by reference the current the mapper object w/ a new parameterization
## Given a set of overlap parameters, compute the mapper instance for each incrementally
# mapper_ref$methods(multiscale = function(overlap_params){
# 
#   # # overlap_params <- seq(0.1, 0.6, by = 0.1)
#   # # overlap_params <- cbind(overlap_params, overlap_params)
#   # ## Order the overlap parameters in non-decreasing order; save the ordering for return
#   # # input_order <- order(overlap_params)
#   # n_pts <- nrow(m$X)
#   # n_lvl_sets <- length(m$cover$level_sets)
#   # n_blocks <- nrow(overlap_params) + 1
#   # browser()
#   # 
#   # ## First step: construct a base cover (where overlap == 0)
#   # m$cover$setOverlap(0)
#   # m$cover$setupCover()
#   # 
#   # ## Get the halfspace distances
#   # hf_dist <- m$cover$computeHalfSpaceDistances()
#   # 
#   # ## Get which points lie in a boundary set
#   # boundary_pts <- apply(sapply(1:ncol(m$cover$filter_values), function(d_i) {
#   #   hf_dist$base_lsmi[, d_i] %in% c(1, m$cover$num_intervals[d_i])
#   # }), 1, any)
#   # 
#   # ## For each new overlap parameter, update:
#   # ## 1) The points within the level set
#   # ## 2) The pairwise distances of the points (in the data space) from (1)
#   # ## 3) The nodes and edges in the graph
#   # target_info <- m$cover$computeTargetLSFI()
#   # 
#   # ## Helper function to determine if a given level set lies on the boundary
#   # isBoundaryLS <- function(lsfi, specific_dim = NULL){
#   #   ls_lsmi <- m$cover$getLSMI(lsfi)
#   #   if (missing(specific_dim) || is.null(specific_dim)){
#   #     res <- apply(ls_lsmi, 1, function(lsmi) { any(sapply(1:filter_dim, function(d_i){ lsmi[d_i] %in% c(1, m$cover$num_intervals[d_i]) })) })
#   #   } else {
#   #     res <- apply(ls_lsmi, 1, function(lsmi) { lsmi[specific_dim] %in% c(1, m$cover$num_intervals[specific_dim]) })
#   #   }
#   #   return(res)
#   # }
#   # 
#   # ## Compute the overlap percentage necessary to merge with each of the target level sets
#   # ## This effectively, per point, the overlap needed for the target level set to intersect the current point
#   # target_box_sizes <- lapply(1:filter_dim, function(d_i){
#   #   t(sapply(1:n_pts, function(i){
#   #     p <- target_info[[d_i]]$delta[i,] ## offset vector
#   #     # boundary_constant <- m$cover$num_intervals[d_i]-1 # ifelse(isBoundaryLS(target_info[[d_i]]$target_lsfi[1,], specific_dim = d_i), 1, m$cover$num_intervals[d_i]-1) ## internal boxes expand k-times as fast
#   #     base_box_sizes <- abs(p)*hf_dist$base_interval[[d_i]]
#   #     hf_dir_dist <- ifelse(p < 0, hf_dist$halfspace_dist[[d_i]][i, 1], hf_dist$halfspace_dist[[d_i]][i, 2])
#   #     boundary_constant <- num_intervals[d_i] - abs(p)
#   #     box_sizes <- base_box_sizes + hf_dir_dist*boundary_constant ## Multiply by boundary constant
#   #     box_sizes
#   #   }))
#   # })
#   # 
#   # all_box_sizes <- cbind(target_box_sizes[[1]][10, ], target_box_sizes[[2]][10, ])
#   # wut <- m$cover$computeOverlap(all_box_sizes)
#   # m$cover$setOverlap(wut[2,2])
#   # m$cover$setupCover()
#   # m$cover$plotFilterSpace(show_ls_bounds = TRUE, show_lsfi = TRUE)
#   # points(m$cover$filter_values[10, 1], m$cover$filter_values[10, 2], col = "red") ## needs to be near 0.955
#   # 
#   # diff(m$cover$level_sets[[25]]$bounds[, 2])
#   # 
#   # ## Say we want to compute the distances to the nearest box
#   # nearer_hf_dist <- cbind(apply(hf_dist$halfspace_dist[[1]], 1, min), apply(hf_dist$halfspace_dist[[2]], 1, min))
#   # nearest_box_sizes <- lower_hf_dist + ifelse(boundary_pts, hf_dist$base_interval, hf_dist$base_interval/2)
#   # 
#   # m$cover$setOverlap(0)
#   # m$cover$setupCover()
#   # m$cover$plotFilterSpace(show_ls_bounds = TRUE, show_lsfi = TRUE)
#   # points(m$cover$filter_values[10, 1], m$cover$filter_values[10, 2], col = "red")
#   # 
#   # lines(x = c(m$cover$filter_values[1, 1], m$cover$filter_values[1, 1]+nearer_hf_dist[1,1]),
#   #      y = c(m$cover$filter_values[1, 2], m$cover$filter_values[1, 2]),
#   #      col = "red")
#   # 
#   # test_box_size <- matrix(c(hf_dist$base_interval[1] + nearer_hf_dist[1,1]*4, hf_dist$base_interval[2]), ncol = 2)
#   # wut <- m$cover$computeOverlap(box_sizes = test_box_size)
#   # m$cover$setOverlap(wut)
#   # m$cover$setupCover()
#   # m$cover$plotFilterSpace(show_ls_bounds = TRUE, show_lsfi = TRUE)
#   # points(m$cover$filter_values[1, 1], m$cover$filter_values[1, 2], col = "red")
#   # abline(v = m$cover$level_sets[[2]]$bounds[1, 1], col = "red")
#   # abline(v = m$cover$level_sets[[2]]$bounds[2, 1], col = "red")
#   # 
#   # ## Get the base level set flat indices
#   # base_lsfi <- m$cover$getLSFI(hf_dist$base_lsmi)
#   # g <- m$cover$computeOverlap()           ## Get the set of possible overlap values
#   # # overlap_order <- order(g$overlap)       ## Order them in increasing order
#   # # g <- lapply(g, function(g_comp) matrix(g_comp[overlap_order,])) ## Reorder all components of g
#   # 
#   # ## Detect which points are going to intersect each level set, indexed by the level set flat index (lsfi)
#   # max_overlap <- apply(overlap_params, 2, max)
#   # filter_dim <- ncol(m$cover$filter_values)
#   # 
#   # ## Sort G in increasing order
#   # for (d_i in 1:filter_dim){
#   #   overlap_order <- apply(g[[d_i]]$overlap, 1, order)
#   #   g[[d_i]]$target_lsfi <- t(sapply(1:n_pts, function(i){
#   #     g[[d_i]]$target_lsfi[i, overlap_order[, i]]
#   #   }))
#   #   g[[d_i]]$overlap <- t(apply(g[[d_i]]$overlap, 1, sort))
#   # }
#   # 
#   # 
#   # ## TODO: wrap this in C++
#   # ## Create  a list of point indices, per level set, which will intersect the given level set at some overlap value
#   # ## les than the current maximum overlap requested
#   # new_ls_pt_intersects <- vector(mode = "list", length = n_lvl_sets)
#   # for (d_i in 1:filter_dim){
#   #   for (i in 1:n_pts){
#   #     targets <- g[[d_i]]$target_lsfi[i, (g[[d_i]]$overlap[i,] < max_overlap[d_i]) ] ## only consider points below the max overlap
#   #     for (lsfi in targets){
#   #       new_ls_pt_intersects[[lsfi]] <- c(new_ls_pt_intersects[[lsfi]], i)
#   #     }
#   #   }
#   # }
#   # 
#   # ## For each level set index, get the points that *will* intersect that level set
#   # ## before or at the maximum overlap value needed.
#   # # new_ls_pt_intersects <- lapply(1:length(m$cover$level_sets), function(lsfi) {
#   # #   as.vector(unlist(sapply(1:filter_dim, function(d_i) {
#   # #     sapply(1:nrow(g[[d_i]]$overlap, function(i){
#   # #       any(g[[d_i]]$target_lsfi[i,] == lsfi)
#   # #     })
#   # #   apply(g[[1]]$target_lsfi, 1, function(lsfis) any(lsfis == lsfi))
#   # #
#   # #     which(g$ls_idx[, d_i + 1] == lsfi & g$overlap[, d_i] <= max_overlap[d_i])
#   # #   })))
#   # # })
#   # 
#   # # new_ls_pt_intersects <- lapply(1:length(m$cover$level_sets), function(lsfi) {
#   # #   ## For each level set index, get the points that *will* intersect that level set
#   # #   ## before or at the maximum overlap value needed.
#   # #   as.vector(unlist(sapply(1:filter_dim, function(d_i) {
#   # #     which(g$ls_idx[, d_i + 1] == lsfi & g$overlap[, d_i] <= max_overlap[d_i])
#   # #   })))
#   # # })
#   # # ls_intersections <- lapply(1:length(m$cover$level_sets), function(lsfi) overlap_order[which(g$ls_idx[, 2] == lsfi & g$overlap < max_overlap)]) ## which points
#   # # full_intersections <- mapply(function(ls, additional_pts) { union(ls$points_in_level_set, additional_pts) }, m$cover$level_sets, ls_intersections)
#   # 
#   # ## Pre-compute the level set distances under the *highest* overlap
#   # ## Although there may be an asymptotically lower runtime associated with appending distances by reference
#   # ## as points begin to intersect level set, computationally it's much more efficient to simply pre-compute the
#   # ## dissimilarities for the highest overlap needed, then subset by index
#   # ls_dists <- mapply(function(lvl_set, pt_idx){
#   #   max_pt_idx <- union(lvl_set$points_in_level_set, pt_idx)
#   #   if ("dist" %in% class(X)){ dist_subset(dist = X, idx = max_pt_idx) }
#   #   else { dist(X[max_pt_idx,], method = "euclidean") }
#   # }, m$cover$level_sets, new_ls_pt_intersects)
#   # 
#   # ## To get the Mapper construction per given parameter, first need to isolate which points and corresponding level sets
#   # ## need updating *per construction*. I call these groups of information which need to be incrementally updated 'blocks'
#   # # overlap_cuts <- sapply(1:filter_dim, function(d_i) findInterval(g$overlap[, d_i], vec = overlap_params[,d_i])) ## block by given overlap parameters
#   # block_ids <- sort(unique(as.vector(overlap_blocks)))
#   # overlap_cuts <- lapply(1:filter_dim, function(d_i){
#   # 
#   #   ## Find which block the computed overlap values correspond to (earliest)
#   #   t(apply(g[[d_i]]$overlap, 1, function(overlap) findInterval(overlap, vec = overlap_params[, d_i])))
#   # })
#   # 
#   #   # pts_to_update <- lapply(block_ids, function(b_i){
#   #   #   which(apply(overlap_blocks, 1, function(pt_mem){ b_i %in% pt_mem }))
#   #   # })
#   #   # names(pts_to_update) <- as.character(block_ids)
#   # 
#   #   ## For each block, extract:
#   #   ## 1) A named list, one for each level set, of point indices that are being unioned with that level set in this block update
#   #   update_tables <- lapply(1:filter_dim, function(d_i){
#   #     k_di <- ncol(g[[d_i]]$target_lsfi)
#   #     pt_idx <- t(sapply(1:n_pts, function(i){ rep(i, k_di) }))
#   #     data.table::data.table(pt_idx = as.vector(t(pt_idx)), lsfi = as.vector(t(g[[d_i]]$target_lsfi)), block_id = as.vector(t(overlap_blocks)))
#   #   })
#   # 
#   #   ## Find which level set each point lies within under the case the given filter space is
#   #   ## zero-dimensional w.r.t. its Lebesgue covering dimension
#   #   num_intervals <- m$cover$num_intervals
#   #   original_level_sets <- sapply(1:filter_dim, function(i) {
#   #     x <- m$cover$filter_values[, i]
#   #     findInterval(x = x, seq(min(x), max(x), length.out = num_intervals[[i]]+1), all.inside = TRUE)
#   #   })
#   #   orig_lsfi <- m$cover$getLSFI(original_level_sets)
#   # 
#   # 
#   #   to_update <- lapply(1:length(block_ids), function(i) vector(mode = "list", length = n_lvl_sets))
#   #   names(to_update) <- as.character(block_ids)
#   #   added_so_far <- vector(mode = "list", length = n_lvl_sets)
#   #   for (b_i in block_ids){
#   #     for (pt_id in 1:n_pts){
#   #       to_lsfi <- sapply(1:filter_dim, function(d_i){ update_tables[[d_i]][pt_idx == pt_id & block_id == b_i]$lsfi })
#   #       if (length(unlist(to_lsfi)) == 0){ next; }
#   #       to_lsfi <- matrix(to_lsfi, ncol = filter_dim)
#   #       update_lsfi <- lapply(1:nrow(to_lsfi), function(i){  m$cover$getRectLSFI(from = orig_lsfi[pt_id], to = to_lsfi[i,]) })
#   #       for (lsfi in setdiff(unique(unlist(update_lsfi)), orig_lsfi[pt_id])){
#   #         c_pts <- to_update[[as.character(b_i)]][[lsfi]]
#   #         pts_to_consider <- union(c(c_pts, pt_id)) ## Add the point to level set
#   #         new_pts <- setdiff(pts_to_consider, added_so_far[[lsfi]]) ## make sure it hasn't been added before
#   #         to_update[[as.character(b_i)]][[lsfi]] <- new_pts
#   #       }
#   #     }
#   #     for (lsfi in 1:n_lvl_sets){
#   #       update_pts <- union(added_so_far[[lsfi]], to_update[[as.character(b_i)]][[lsfi]])
#   #       if (!is.null(update_pts)){ added_so_far[[lsfi]] <- update_pts }
#   #     }
#   #   }
#   # 
#   #   lvl_set_of_interest <- 25
#   #   pt_colors <- rainbow(length(block_ids))
#   #   m$cover$plotFilterSpace(show_ls_bounds = TRUE, show_lsfi = TRUE)
#   #   for (b_i in 1:length(block_ids)){
#   #     new_pts <- to_update[[b_i]][[lvl_set_of_interest]]
#   #     if (!is.null(new_pts)){
#   #       points(m$cover$filter_values[new_pts, 1], m$cover$filter_values[new_pts, 2], col = pt_colors[b_i])
#   #     }
#   #     invisible(readline(prompt="Press [enter] to continue"))
#   #   }
#   # 
#   # 
#   #   to_update <- lapply(block_ids, function(b_i){
#   #     tmp <- ls_intersects_by_block[block_id == b_i]
#   #     uniq_lsfi <- unique(tmp$lsfi)
#   #     res <- lapply(uniq_lsfi, function(lsfi_i){ tmp[lsfi == lsfi_i]$pt_idx })
#   #     names(res) <- as.character(uniq_lsfi)
#   #     return(res)
#   #   })
#   #      # lvl_sets <- g[[d_i]]$target_lsfi[overlap_blocks == b_i] ## Which level sets need to be updated in this block (dimension-specific)
#   #     # pts_to_update <- which(apply(overlap_blocks == b_i, 1, any))
#   #     # print(length(lvl_sets))
#   #     # print(length(pts_to_update))
#   #     # data.table::data.table(lsfi = lvl_sets, pt_idx = pts_to_update)
#   #     # pts_in_lvl_sets <- lapply(lvl_sets, function(lsfi){
#   #     #   which(apply(overlap_blocks == lsfi, 1, any))
#   #     # })
#   #     # names(pts_in_lvl_sets) <- lvl_sets
#   #     # Filter(function(pt_membership) length(pt_membership) > 0, pts_in_lvl_sets)
#   # 
#   # 
#   # 
#   # #   orig_lsfi <- m$cover$getLSFI(original_level_sets)
#   # #   sapply(1:length(orig_lsfi), function(i){
#   # #     targets <- cbind(g[[1]]$target_lsfi[i,], g[[2]]$target_lsfi[i,])
#   # #     apply(targets, 1, function(pair) m$cover$getRectLSFI(orig_lsfi[i], pair))
#   # #   })
#   # #   wut <- cbind(g[[1]]$target_lsfi[1,], g[[2]]$target_lsfi[1,])
#   # #
#   # #   # m$cover$getRectLSFI(orig_lsfi, wut[1, ])
#   # #
#   # #   to_update
#   # # })
#   # 
#   # 
#   # block_ids <- (1:n_blocks) - 1L
#   # sapply(1:filter_dim, function(d_i){ apply(g[[d_i]]$overlap, 1, function(overlap) findInterval(overlap, vec = overlap_params[, d_i])) })
#   # lapply(block_ids, function(){
#   #   sapply(1:filter_dim, function(d_i){
#   #     apply(g[[d_i]]$overlap, 1, function(overlap) findInterval(overlap, vec = overlap_params[, d_i]))
#   #   })
#   # })
#   # 
#   # ## Remove overlaps corresponding to invalid target level sets
#   # overlap_cuts[which(is.na(g$ls_idx[, -1]))] <- NA
#   # block_ids <- sort(unique(as.integer(overlap_cuts)))
#   # 
#   # ## Utility function to determine which points need to be updated for a given block id
#   # need_updated <- function(block_id){
#   #   which(apply(overlap_cuts, 1, function(active_blocks) any(block_id %in% active_blocks)))
#   # }
#   # 
#   # ## Aggregate all of the information needed to incrementally construct Mapper constructions. This includes
#   # ## 1) A named list for each level set containing the point ids which will intersect that level set for the given block
#   # ## 2) The upper bound on the overlap associated with this block
#   # to_update <- lapply(block_ids, function(block_id) {
#   #   point_ids <- need_updated(block_id) ## which points must be updated in this block?
#   #   level_set_subset <- g$ls_idx[point_ids,-1] ## which level sets will these new points intersect in this block?
#   #   levels_to_update <- unique(as.vector(na.omit(as.vector(level_set_subset)))) ## what level sets (lsfi's) need to be updated this block?
#   #   pts_to_update_by_ls <- lapply(levels_to_update, function(lsfi){
#   #     as.vector(unlist(apply(level_set_subset, 2, function(target_ls_indexes) {
#   #       point_ids[which(target_ls_indexes == lsfi)]
#   #     })))
#   #   })
#   #   names(pts_to_update_by_ls) <- levels_to_update
#   #   list(points_to_update = pts_to_update_by_ls,
#   #        overlap_upper_bound = overlap_params[block_id+1,])
#   # })
#   # 
#   # ## Compute the new set of mapper instances
#   # all_mappers <- vector(mode = "list", length = length(to_update))
#   # i <- 0L
#   # for (block in to_update){
#   #   for (level in names(block$points_to_update)){
#   #     c_pts <- m$cover$level_sets[[as.integer(level)]]$points_in_level_set
#   # 
#   #     ## Update which points lie in the level set
#   #     new_pts <- union(c_pts, block$points_to_update[[level]])
#   #     m$cover$level_sets[[as.integer(level)]]$points_in_level_set <- new_pts
#   # 
#   #     ## Update the dist object in the level set with the (relatively indexed) newly intersected points
#   #     m$cover$level_sets[[as.integer(level)]]$dist <- TDAmapper:::dist_subset(dist = ls_dists[[as.integer(level)]],
#   #                                                                             idx = which(new_pts %in% new_ls_pt_intersects[[as.integer(level)]]))
#   #   }
#   #   ## Update nodes (Since dissimilarities have been updated, just rerun clustering)
#   #   m$computeNodes(... =  block$levels_to_update)
#   # 
#   #   ## Update edges (but only the by comparing the level sets which need updating)
#   #   levels_to_update <- TDAmapper:::valid_pairs(g$ls_idx[unname(unlist(block$points_to_update)),])
#   #   m$computeEdges(level_sets = block$levels_to_update)
#   # 
#   #   ## Story copy of object
#   #   all_mappers[[i]] <- m$G # m$copy(shallow = FALSE)
#   #   i <- i + 1
#   # }
# 
# })

## Exports the internal mapper core structures to a TDAmapper output
mapper_ref$methods(exportTDAmapper = function(){
  level_of_vertex <- sapply(G$nodes, function(ni) attr(ni, "level_set"))
  structure(list( adjacency = G$adjacency,
                  num_vertices = length(G$nodes),
                  level_of_vertex = level_of_vertex,
                  points_in_vertex = lapply(G$nodes, as.vector),
                  points_in_level = unname(lapply(cover$level_sets, function(ls) ls$points_in_level_set)),
                  vertices_in_level = lapply(1:length(cover$level_sets), function(ls_idx) {
                    tmp <- which(level_of_vertex == ls_idx)
                    if (length(tmp) == 0){ return(-1) }
                    else return(tmp)
                  })),
                  class = "TDAmapper")
})

## Exports the internal mapper core structures to a Mapper object, which is a list containing:
## 1. nodes member := list of integer vectors representing the indices of the original data the current node intersects with. Also
##                    contains attribute data storing the level set flat index of the level set the node is in.
## 2. adjacency := adjacency matrix of the resulting graph.
mapper_ref$methods(exportMapper = function(){
  result <- .self$G
  node_lsfi <- sapply(result$nodes, function(node) attr(node, "level_set"))
  result$level_sets <- lapply(1:length(cover$level_sets), function(lsfi){ which(node_lsfi == lsfi) })
  names(result$level_sets) <- apply(cover$index_set, 1, function(lsmi) paste0("(", paste(lsmi, collapse = ","), ")"))
  cover_type <- paste0(toupper(substr(cover$type, start = 1, stop = 1)), tolower(substr(cover$type, start = 2, stop = nchar(cover$type))))
  z_d <- ncol(cover$filter_values)
  attr(result, ".summary") <- c(sprintf("Mapper object with filter function f: %s -> %s",
                                        ifelse(is(X, "dist"), "dist(X)",
                                               ifelse(ncol(X) > 1, sprintf("X^%d", ncol(X)), "X")),
                                        ifelse(z_d > 1, sprintf("Z^%d", z_d), "Z")),
                                sprintf("The codomain of f is equipped with a %s cover composed of %d open sets", cover$type, length(cover$level_sets)),
                                sprintf("The 1-skeleton contains %d nodes and %d edges", length(G$nodes), sum(G$adjacency == 1L)/2L))
  class(result) <- "Mapper"
  return(result)
})



