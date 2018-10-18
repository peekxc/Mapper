#' Enable multiscale 
#' @description 
#' After calling this function, a new function 'update_mapper' becomes available as part of the MapperRef instance,
#' which accepts as its argument a single parameter 'percent_overlap', and returns the corresponding Mapper for that overlap parameterization. 
#' 
#' @details 
#' Calling this function produces an internal indexing structure capable of computing Mappers efficiently
#' between overlap parameterizations. Only supported for the the fixed rectangular cover. It's assumed the 
#' number of intervals is fixed--if the number of intervals to distribute along any dimension changes, this 
#' function must be recalled.  
#' 
#' @examples 
#'   data("noisy_circle")
#'   ll <- noisy_circle[which.min(apply(noisy_circle, 1, sum)),]
#'   f_x <- t(apply(noisy_circle, 1, function(x) abs(x - ll)))
#'   m_ref <- Mapper::mapper(noisy_circle, filter_values = f_x)
#' @export
MapperRef$set("public", "enable_multiscale", 
  function(){
    # if (self$cover$type != "rectangular"){ stop("'multiscale' is only compatible with rectangular covers.") }
    n <- ifelse("dist" %in% class(private$.X), attr(private$.X, "Size"), nrow(private$.X))
    d <- ncol(self$cover$filter_values)
    
    ## Construct base cover
    self$cover$percent_overlap <- 0.0
    self$cover$construct_cover()

    ## Retrieve the bounds on each level set 
    ls_bnds <- self$cover$level_set_bounds()
    
    ## Get the level set flat indices associated with each point
    lsfi <- Mapper:::constructLevelSetIndex(self$cover$filter_values, ls_bnds) ## returns 1-based
    
    ## Translate the flat index into a multi-index to encode the points positioning in the level set
    cart_prod <- arrayInd(seq(prod(self$cover$number_intervals)), .dim = self$cover$number_intervals)
    A <- matrix(cart_prod[lsfi,], ncol = d)
    
    ## Get filter min and max ranges
    filter_rng <- apply(self$cover$filter_values, 2, range)
    { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
    filter_len <- diff(filter_rng)
    base_interval_length <- filter_len/self$cover$number_intervals
    
    ## Use the encoding to translate or 'stack' each point onto a single level set 
    A_tmp <- apply(A, 1, function(a) as.vector(as.integer(a) - 1L) * as.vector(base_interval_length))
    Z_tmp <- apply(self$cover$filter_values, 1, function(z_i) as.numeric(z_i) - filter_min)
    Z_tilde <- matrix(Z_tmp - A_tmp, nrow = d)
    
    ## Compute the distances to the lower and upper level sets
    ## This is equivalent to detecting which halfspace each point z_i lies in, and getting the distance to the boundary
    dist_to_upper_ls <- matrix(apply(Z_tilde, 2, function(z_i) abs(as.vector(z_i) - as.vector(base_interval_length))), nrow = d)
    dist_to_lower_ls <- Z_tilde
    
    ## Use these distances to compute the distances to the rest of the 'target' level sets in each dimension
    dist_to_ls <- lapply(1:d, function(d_i){
      Mapper:::dist_to_boxes(
        positions = as.integer(A[, d_i]), ## column indexed
        interval_length = as.numeric(base_interval_length[d_i]), 
        num_intervals = as.integer(self$cover$number_intervals[d_i]), 
        dist_to_lower = as.numeric(dist_to_lower_ls[d_i,]),  ## row indexed
        dist_to_upper = as.numeric(dist_to_upper_ls[d_i,]) ## row indexed
      )
    })
    
    ## Order the distances to the target level sets
    dist_order <- lapply(1:d, function(d_i) order(dist_to_ls[[d_i]]$target_dist))
    
    ## Get the corresponding interval sizes 
    ## TODO: generalize this to other rectangular covers
    interval_sizes <- lapply(1:d, function(d_i) dist_to_ls[[d_i]]$target_dist[dist_order[[d_i]]]*2 + base_interval_length[d_i])
    
    ## Checks
    if (any(is.na(unlist(interval_sizes)))){ 
      stop("Something went wrong with the interval size calculation.")
    }
    
    ## At what interval size(s) do the level sets themselves intersect other level sets?
    intersection_cuts <- lapply(1:d, function(d_i) {
      swap_dist <- ((base_interval_length[d_i]/2) * seq(2, self$cover$number_intervals[d_i] - 1L))*2
      swap_idx <- findInterval(interval_sizes[[d_i]], vec = swap_dist)
      tmp <- rle(swap_idx)[["lengths"]]
      head(c(0L, cumsum(tmp)), length(tmp))
    })
    
    ## Extract all the unique level set paths for each point
    ls_paths <- vector(mode = "list", length = d)
    uniq_ls_paths <- vector(mode = "list", length = d)
    for (d_i in 1L:d){
      pt_ls_paths <- do.call(rbind, lapply(1:n, function(i) {
        c(A[i, d_i], with(dist_to_ls[[d_i]], { target_pos[i, order(target_dist[i,])] }))
      }))
      uniq_ls_paths[[d_i]] <- unique(pt_ls_paths)
      ls_paths[[d_i]] <- rowmatch(pt_ls_paths, uniq_ls_paths[[d_i]])
    }
    
    ## Create the indexing structure
    private$.multiscale <- Mapper:::MultiScale$new(n, as.vector(self$cover$number_intervals))
    private$.multiscale$insert_pts(A, do.call(cbind, ls_paths))
    for (d_i in 1L:d){
      private$.multiscale$create_filtration(dist_order[[d_i]]-1L, interval_sizes[[d_i]], d_i-1)
      private$.multiscale$set_filtration_rle(intersection_cuts[[d_i]], d_i-1)
      private$.multiscale$create_ls_paths(uniq_ls_paths[[d_i]]-1L, d_i-1)
    }
    
    ## Add the ability to 'update' the Mapper by changing its overlap.
    ## If stats=TRUE, returns information related to how the Mapper changes w.r.t the previous Mapper
    ## TODO: Move this out of the current environment
    self$update_mapper <- function(percent_overlap, stats=FALSE){
      stopifnot(all(percent_overlap >= 0 && percent_overlap < 1))
      stopifnot(length(percent_overlap) == ncol(self$cover$filter_values))
      
      ## If not populated yet, create a default mapping from the covers index set to the vertex ids
      if (length(private$.cl_map) == 0){
        private$.cl_map <- lapply(self$cover$index_set, function(x) integer(0L))
        names(private$.cl_map) <- self$cover$index_set
      }
      
      # browser()
      ## TODO generalize this
      R <- base_interval_length + (base_interval_length*percent_overlap)/(1.0 - percent_overlap)
      idx <- private$.multiscale$get_nearest_filtration_index(R)
      
      ## Update the segments 
      res <- private$.multiscale$update_segments(idx)
      
      ## Etxract the simplex tree 
      stree_ptr <- private$.simplicial_complex$as_XPtr()
      
      ## Detect the edges between level sets affected by the change
      connected_ls <- lapply(res$ls_to_update, function(ls){
        # browser()
        adj_lsfis <- check_connected(ls, private$.cl_map, self$vertices, stree_ptr)
        if (length(adj_lsfis) > 0){
          cbind(pmin(ls, adj_lsfis-1L), pmax(ls, adj_lsfis-1L))
        } else { NULL }
      })
      
      if (stats){
        old_vertices <- self$ls_vertex_map[(res$ls_to_update+1L)]
        old_ls_map <- self$ls_vertex_map
      }
      
      ## Update the vertices directly, without storing the level sets explicitly. 
      if (length(res$ls_to_update) > 0){
        ms_ptr <- private$.multiscale$as_XPtr()
        private$.vertices <- update_level_sets(
          which_levels = res$ls_to_update,
          ms = ms_ptr, 
          X = private$.X, 
          f = self$clustering_algorithm, 
          vertices = private$.vertices, 
          ls_vertex_map = private$.cl_map, 
          stree = stree_ptr
        )
      }
      
      ## The pairs to update
      ls_pairs_to_update <- unique(rbind(res$ls_pairs_to_update, do.call(rbind, connected_ls)))
      
      ## Update the 1-skeleton
      if (nrow(ls_pairs_to_update) > 0){
        self$compute_edges(ls_pairs_to_update+1)
      }
      
      ## Returns statistics about how the Mapper changed if requested, self o.w.
      if (stats){ 
        res_stats <- list(
          old_vertices=unname(unlist(old_vertices)), 
          new_vertices=unname(unlist(self$ls_vertex_map[(res$ls_to_update+1L)])), 
          old_ls_map = old_ls_map, 
          new_ls_map = self$ls_vertex_map, 
          updated_ls = res$ls_to_update, 
          updated_ls_pairs = ls_pairs_to_update
        )
        return(res_stats)
      } 
      else { return(self) }
    }
    
    return(self)
  }
) # multiscale 
    
    ## ---------------------------------------
    ## Below this point, all non-point *indexes* should be 0-based
    ## ---------------------------------------
    
    ## TODO: Can the memory requirements with this be improved?
    ## Use the ordering to extract the point indices, the level set indices the points are in, the level set indices the 
    ## points will be intersecting, and the distance they'll be starting the intersection
    # 
    # pt_idx <- lapply(1:filter_dim, function(d_i) row(dist_to_ls[[d_i]]$target_pos)[dist_order[[d_i]]]) ## pts are still 1-based
    # from_ls <- lapply(1:filter_dim, function(d_i) matrix(rep(A[, d_i], each = number_intervals[d_i]-1), ncol = number_intervals[d_i] - 1, byrow = TRUE)[dist_order[[d_i]]] - 1L)
    # to_ls <- lapply(1:filter_dim, function(d_i) dist_to_ls[[d_i]]$target_pos[dist_order[[d_i]]] - 1L)
    # 
    # #unique(do.call(rbind, lapply(1:nrow(dist_to_ls[[d_i]]$target_pos), function(i) dist_to_ls[[d_i]]$target_pos[i,order(dist_to_ls[[d_i]]$target_dist[i,])])))
    # 
    # ## == The following equations require knowledge of how the cover boxes are constructed ==
    # 
    # 
    # 
    # 
    # ## Stop here 
    # 
    # 
    # 
    # library("Mapper")
    # data("noisy_circle")
    # X <- noisy_circle[1:10,]
    # # X <- cbind(rnorm(1000), rnorm(1000), rnorm(1000))
    # left_pt <- X[which.min(X[, 1]),]
    # f_x <- apply(X, 1, function(pt) (pt - left_pt)[1])
    # # filter_values <- matrix(f_x)
    # filter_values <- cbind(f_x, (max(f_x) - f_x^2) + rnorm(length(f_x), sd = 2))
    # 
    # ## Extract the cartesian product of the level sets (index set)
    # indices <- lapply(number_intervals, seq) ## per-dimension possible indexes
    # cart_prod <- structure(matrix(as.integer(unlist(do.call(expand.grid, indices))), ncol = filter_dim), 
    #                        dimnames = list(NULL, paste0("d", 1:filter_dim)))
    # 
    # ## Get filter min and max ranges
    # filter_rng <- apply(filter_values, 2, range)
    # { filter_min <- filter_rng[1,]; filter_max <- filter_rng[2,] }
    # filter_len <- diff(filter_rng)
    # 
    # percent_overlap <- rep(0.0, filter_dim)
    # 
    # ## Construct the level sets
    # base_interval_length <- filter_len/number_intervals
    # eps <- base_interval_length/2.0
    # ls_bnds <- t(apply(cart_prod, 1, function(idx){
    #   centroid <- filter_min + ((as.integer(idx)-1L)*base_interval_length) + base_interval_length/2.0
    #   c(centroid - eps - sqrt(.Machine$double.eps), centroid + eps + sqrt(.Machine$double.eps))
    # }))
    # endpts <- lapply(1L:filter_dim, function(d_i){
    #   sort(as.vector(sapply(0L:(number_intervals[d_i] - 1L), function(idx){
    #     centroid <- filter_min[d_i] + (as.integer(idx)*base_interval_length[d_i]) + base_interval_length[d_i]/2.0
    #     c(centroid - eps[d_i], centroid + eps[d_i])
    #   })))
    # })
    # 
    # 
    # 
    # lsfi <- Mapper:::constructLevelSetIndex(filter_values, ls_bnds) ## returns 1-based
    # A <- matrix(cart_prod[lsfi,], ncol = filter_dim)
    # A_tmp <- apply(A, 1, function(a) (as.integer(a) - 1L) * base_interval_length)
    # Z_tmp <- apply(filter_values, 1, function(z_i) as.numeric(z_i) - filter_min)
    # Z_tilde <- matrix(Z_tmp - A_tmp, nrow = filter_dim)
    # 
    # ## Distance to upper and lower level sets (halfspace distances)
    # dist_to_upper_ls <- matrix(apply(Z_tilde, 2, function(z_i) abs(as.vector(z_i) - as.vector(base_interval_length))), nrow = filter_dim)
    # dist_to_lower_ls <- Z_tilde
    # 
    # ## Applied to all level sets, per dimension 
    # dist_to_ls <- lapply(1:filter_dim, function(d_i){
    #   Mapper:::dist_to_boxes(
    #     positions = as.integer(A[, d_i]), ## column indexed
    #     interval_length = as.numeric(base_interval_length[d_i]), 
    #     num_intervals = as.integer(number_intervals[d_i]), 
    #     dist_to_lower = as.numeric(dist_to_lower_ls[d_i,]),  ## row indexed
    #     dist_to_upper = as.numeric(dist_to_upper_ls[d_i,]) ## row indexed
    #   )
    # })
    # 
    # ## Order the distances to the target level sets
    # dist_order <- lapply(1:filter_dim, function(d_i) order(dist_to_ls[[d_i]]$target_dist))
    # 
    # ## ---------------------------------------
    # ## Below this point, all non-point *indexes* should be 0-based
    # ## ---------------------------------------
    # 
    # ## TODO: Can the memory requirements with this be improved?
    # ## Use the ordering to extract the point indices, the level set indices the points are in, the level set indices the 
    # ## points will be intersecting, and the distance they'll be starting the intersection
    # 
    # pt_idx <- lapply(1:filter_dim, function(d_i) row(dist_to_ls[[d_i]]$target_pos)[dist_order[[d_i]]]) ## pts are still 1-based
    # from_ls <- lapply(1:filter_dim, function(d_i) matrix(rep(A[, d_i], each = number_intervals[d_i]-1), ncol = number_intervals[d_i] - 1, byrow = TRUE)[dist_order[[d_i]]] - 1L)
    # to_ls <- lapply(1:filter_dim, function(d_i) dist_to_ls[[d_i]]$target_pos[dist_order[[d_i]]] - 1L)
    # 
    # #unique(do.call(rbind, lapply(1:nrow(dist_to_ls[[d_i]]$target_pos), function(i) dist_to_ls[[d_i]]$target_pos[i,order(dist_to_ls[[d_i]]$target_dist[i,])])))
    # 
    # ## == The following equations require knowledge of how the cover boxes are constructed ==
    # ## What interval size does each parameterization represent?
    # interval_size <- lapply(1:filter_dim, function(d_i) dist_to_ls[[d_i]]$target_dist[dist_order[[d_i]]]*2 + base_interval_length[d_i]) ## boxes intersect 
    # 
    # ## What interval sizes do the level sets change in their relative segment order?
    # swap_dist <- lapply(1:filter_dim, function(d_i) ((base_interval_length[d_i]/2) * seq(2, number_intervals[d_i] - 1L))*2)
    # swap_idx <- lapply(1:filter_dim, function(d_i) findInterval(interval_size[[d_i]], vec = swap_dist[[d_i]]) + 1L)
    # 
    # ## More preprocessing...
    # rowmatch <- function(A,B) { 
    #   # Rows in A that match the rows in B
    #   f <- function(...) paste(..., sep=":")
    #   if(!is.matrix(B)) B <- matrix(B, 1, length(B))
    #   a <- do.call("f", as.data.frame(A))
    #   b <- do.call("f", as.data.frame(B))
    #   match(a, b)
    # }
    # 
    # n <- nrow(A)
    # ls_paths <- vector(mode = "list", length = filter_dim)
    # uniq_ls_paths <- vector(mode = "list", length = filter_dim)
    # for (d_i in 1L:filter_dim){
    # 
    #   ## Extract all the unique level set paths for dimension d_i
    #   pt_ls_paths <- do.call(rbind, lapply(1:n, function(i) {
    #     c(A[i, d_i], with(dist_to_ls[[d_i]], { target_pos[i, order(target_dist[i,])] }))
    #   }))
    #   uniq_ls_paths[[d_i]] <- unique(pt_ls_paths)
    #   ls_paths[[d_i]] <- rowmatch(pt_ls_paths, uniq_ls_paths[[d_i]])
    # }
    # 
    # ## Testing second slimmed down version 
    # ms_mapper2 <- Mapper:::MultiScale2$new(n, as.vector(number_intervals))
    # for (d_i in 1L:filter_dim){
    #   ms_mapper2$create_filtration(dist_order[[d_i]]-1L, d_i-1)
    #   ls_idx_swaps <- rle(swap_idx[[d_i]]-1L)[["lengths"]]
    #   ls_idx_swaps <- head(c(0L, cumsum(ls_idx_swaps)), length(ls_idx_swaps))
    #   ms_mapper2$set_filtration_rle(ls_idx_swaps, d_i-1)
    #   ms_mapper2$create_ls_paths(uniq_ls_paths[[d_i]]-1L, d_i-1)
    # }
    # ms_mapper2$insert_pts(A, do.call(cbind, ls_paths))
    # 
    # plot_2d(idx = c(4L, 4L), pt_ids = TRUE)
    # 
    # plot_configuration(d_i = 1L, idx = 6L)
    # ms_mapper2$update_segments2(5L)
    # Mapper:::index_test(c(7L, 7L))
    # 
    # ms_mapper2$update_segments2(c(5,5))
    # plot_configuration(idx = 1, d_i = 1)
    # 
    # 
    # ## Given a parameterization, determine the index of the filtration, then pass that index down to the multiscale structure. 
    # ## The structure should then use its structure to compute which of the level sets are changing. 
    # f <- function(X, idx){ return(rep(1L, length(idx))) }
    # ms_mapper <- Mapper:::MultiScale$new(X = X, f = identity, k = number_intervals)
    # for (d_i in 1L:filter_dim){
    #   ms_mapper$set_everything(pt_idx[[d_i]], from_ls[[d_i]], to_ls[[d_i]], swap_idx[[d_i]], interval_size[[d_i]], swap_dist[[d_i]], d_i-1L)
    # }
    # ms_mapper$build_segment_trees(endpts = endpts, filter_pts = filter_values)
    # ms_mapper$build_multiscale_configuration((A-1L)*2, filter_values);
    # # ms_mapper$get_segment_map()
    # # ms_mapper$update_multi_cover(c(0L, 0L))
    # ms_mapper$compute_mapper(c(0, -1))
    # ms_mapper$compute_mapper(c(0, 0))
    # ms_mapper$compute_mapper(c(1, 0))
    # ms_mapper$compute_mapper(c(1, 1))
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # plot_2d(idx = c(1, 0), pt_ids = TRUE)
    # 
    # ms_mapper$get_ls_that_change(c(0, 0), c(1, 0))
    # 
    # wut <- t(outer(4:0, 3:0, Vectorize(function(x, y){ ms_mapper$query_ls_map(c(x, y)) })))
    # t(apply(wut, 1, rev))
    # 
    # # mapply(function(d_i, k_i) swap_idx[[d_i]][k_i], 1:filter_dim, c(1, 1))
    # # plot_2d(c(1L, 1L))
    # # data.frame(A, sep="", ms_mapper$get_ls_that_change(c(0, 0), A))
    # 
    # ## Build a example multiscale filtration of mappers 
    # g <- seq(0+0.005, 1-0.005, by = 0.005) ## overlap parameters 
    # target_idx <- sapply(1:filter_dim, function(d_i){
    #   R_hat <- base_interval_length[[d_i]] + (base_interval_length[[d_i]]*g)/(1 - g)
    #   tmp <- sapply(R_hat, function(r_hat) {
    #     idx_lt_thresh <- interval_size[[d_i]] <= r_hat
    #     if (all(idx_lt_thresh == FALSE)){ 0L }
    #     else { max(which(idx_lt_thresh))}
    #   })
    #   tmp
    # })
    # # Iso::biviso(sapply(1:filter_dim, function(d_i){
    # #   base_interval_length[[d_i]] + (base_interval_length[[d_i]]*g)/(1 - g)
    # # }))
    # 
    # 
    # ## Start with a simple Mapper reference object with a constructed, disjoint cover and trivial clustering function
    # m <- MapperRef$new(X)$
    #   use_cover(filter_values = filter_values, type = "fixed rectangular", number_intervals = number_intervals, percent_overlap = 0)$
    #   use_distance_measure("euclidean")
    # m$clustering_algorithm <- function(X, idx, ...){ rep(1L, length(idx)) }
    # 
    # microbenchmark::microbenchmark({ invisible(ms_mapper$get_ls_that_change(target_idx[1,], target_idx[2,])) }, times = 15L)
    # microbenchmark::microbenchmark({
    #   Gs <- vector(mode = "list", length = nrow(target_idx))
    #   for (i in 1:nrow(target_idx)){
    #     t_idx <- target_idx[i,] - 1L
    #     # // ms_mapper$get_ls_that_change(c(1, 0))
    #     ms_mapper$update_multi_cover(t_idx)
    #     ls_tmp <- ms_mapper$get_level_sets()
    #     m$cover$level_sets <- ls_tmp[match(m$cover$index_set, names(ls_tmp))]
    #     ls_config <- sapply(1:filter_dim, function(d_i) ifelse(t_idx[[d_i]] == -1L, 0, swap_idx[[d_i]][t_idx[[d_i]]+1L]))
    #     
    #     ## TODO: either control fro shrinking doing new ls updates, or fix the damn LS updater!
    #     ls_pairs <- ms_mapper$ls_to_change(ls_config)+1L
    #     ls_to_update <- unique(as.vector(ls_pairs))
    #     if (length(ls_to_update) == 0){
    #       m$compute_vertices()
    #     } else {
    #       m$compute_vertices(which_levels = m$cover$index_set[ls_to_update])
    #       m$compute_edges(which_level_pairs = ls_pairs)
    #     }      
    #     Gs[[i]] <- m$simplicial_complex$as_adjacency_matrix()
    #   }
    # }, times = 15L)
    # 
    # ## -- vs. --
    # microbenchmark::microbenchmark({
    #   m2 <- MapperRef$new(X)$
    #     use_cover(filter_values = filter_values, type = "fixed rectangular", number_intervals = number_intervals, percent_overlap = 0)$
    #     use_distance_measure("euclidean")
    #   m2$clustering_algorithm <- function(X, idx, ...){ rep(1L, length(idx)) }
    #   Gs2 <- vector(mode = "list", length = length(g))
    #   for (i in 1:length(g)){
    #     m2$cover$percent_overlap <- g[i]
    #     m2$cover$construct_cover()
    #     m2$compute_k_skeleton(k = 1L)
    #     Gs2[[i]] <- m2$simplicial_complex$as_adjacency_matrix()
    #   }
    # }, times = 15L)
    # 
    # # all.equal(Gs, Gs2)
    # sapply(1:length(g), function(i){
    #   g1 <- igraph::graph_from_adjacency_matrix(Gs[[i]])
    #   g2 <- igraph::graph_from_adjacency_matrix(Gs2[[i]])
    #   all(c(
    #     igraph::vcount(g1) == igraph::vcount(g2), 
    #     igraph::ecount(g1) == igraph::ecount(g2),
    #     sort(igraph::local_scan(g1, k = 1)) == sort(igraph::local_scan(g2, k = 1)), 
    #     sort(igraph::local_scan(g1, k = 2)) == sort(igraph::local_scan(g2, k = 2)), 
    #     sort(igraph::local_scan(g1, k = 3)) == sort(igraph::local_scan(g2, k = 3))
    #   ))
    # })
    # 
    # 
    # 
    # m <- MapperRef$new(noisy_circle)$
    #   use_cover(filter_values = f_x, type = "fixed rectangular", number_intervals = 5L, percent_overlap = 0.35)$
    #   use_clustering_algorithm(cl = "single", num_bins = 10)$
    #   use_distance_measure(measure = "euclidean")$
    #   compute_vertices()$
    #   compute_edges()
    # 
    # m$cover$index_set
    # 
    # ## 1D case
    # ms_mapper$update_multi_cover(c(1L))
    # ms_mapper$update_multi_cover(c(0L))
    # 
    # 
    # plot_2d(c(1, 1))
    # ms_mapper$update_multi_cover(c(0L, -1L))
    # ms_mapper$update_multi_cover(c(1L, -1L))
    # 
    # ms_mapper$update_cover(0L, 0L)
    # ms_mapper$update_cover(1L, 0L)
    # ms_mapper$update_cover(2L, 0L)
    # ms_mapper$update_cover(3L, 0L)
    # 
    # ms_mapper$update_cover(0L, 0L)
    # 
    # 
    # 
    # 
    # mtree <- Mapper:::MultiSegmentTree$new(endpts = endpts)
    # mtree$insert_pts(filter_values)
    # mtree$query_dimension(0, 5, 0)
    # 
    # ms_mapper$compute_ls_idx(i = 1, d_i = 0)
    # 
    # for (i in (length(pt_idx[[1]]) - 1):0){
    #   plot_configuration(1L, i+1)
    #   ms_mapper$update_cover(i, 0L)
    #   print(format(ms_mapper$get_level_sets()))
    #   invisible(readline(prompt="Press [enter] to continue"))
    # }
    # 
    # lsfi_to_lsmi <- vector(mode = "list", length = 5*4)
    # for (x in 4:0){
    #   for (y in 3:0){
    #     lsfi_to_lsmi[[ms_mapper$query_ls_map(c(x, y))+1]] <- c(x, y) # paste0("(", x, " ", y, ")")
    #   }
    # }
    # 
    # worst_case <- rep(choose(prod(number_intervals), 2), nrow(target_idx)-1)
    # rec_based <- sapply(2:nrow(target_idx), function(i){
    #   c_swap <- sapply(1:filter_dim, function(d_i) swap_idx[[d_i]][i])
    #   nrow(ms_mapper$ls_to_change(c_swap))
    # })
    # optimal <- sapply(2:nrow(target_idx), function(i){
    #   nrow(ms_mapper$get_ls_that_change(target_idx[i-1,], target_idx[i,]))
    # })
    # suitable_g <- g[1:(length(g)-1)]
    # plot(x = suitable_g, y = log(cumsum(worst_case)), type = "l", ylim = c(0, max(log(cumsum(worst_case)))), 
    #      xlab = "Overlap percentage", ylab = "Cumulative # of LS intersections (log scale)", 
    #      main = "Complexity of building the 1-skeleton", 
    #      col = "red", lwd = 2)
    # mtext("As a function of g")
    # lines(x = suitable_g, y = log(cumsum(rec_based)), lwd = 2, col = "blue")
    # lines(x = suitable_g, y = log(cumsum(optimal)), lwd = 2, col = "green")
    # 
    # 
    # sapply()
    # ms_mapper$ls_to_change(c(1, 1))
    # 
    # animation::saveGIF({
    #   for (i in 2:200){
    #     plot_2d(idx = c(i, i), prev = c(i-1, i-1))
    #   }
    # }, "boxes_changing_2d.gif", interval = 0.2)
    # 
    # 
    # 
    # make_square <- function(x, y, ...){
    #   lines(x = c(x[1], x[2]), y = c(y[1], y[1]), ...)
    #   lines(x = c(x[2], x[2]), y = c(y[1], y[2]), ...)
    #   lines(x = c(x[2], x[1]), y = c(y[2], y[2]), ...)
    #   lines(x = c(x[1], x[1]), y = c(y[2], y[1]), ...)
    # }
    # get_range <- function(x){ range(x)+c(-1, 1)*diff(range(x))*0.10 }
    # 
    # lsmi <- as.matrix(do.call(expand.grid, lapply(number_intervals, seq)))-1L
    # 
    # plot_2d <- function(idx=NULL, g=NULL, prev=NULL, pt_ids = FALSE){
    #   if (missing(idx) && !missing(g)){
    #     ## Construct the level sets
    #     ls_endpts <- lapply(1L:filter_dim, function(d_i){
    #       ## Choose the interval (half) width to plot 
    #       if (g[d_i] == 0){
    #         eps <- (base_interval_length[[d_i]]/2)
    #       } else {
    #         eps <- (base_interval_length[[d_i]] + ((base_interval_length[[d_i]]*g[d_i])/(1 - g[d_i])))/2 # parameterized overlap
    #       }
    #       tmp <- as.vector(sapply(0L:(number_intervals[d_i] - 1L), function(idx){
    #         centroid <- filter_min[d_i] + (as.integer(idx)*base_interval_length[d_i]) + base_interval_length[d_i]/2.0
    #         c(centroid - eps, centroid + eps)
    #       }))
    #       matrix(tmp, ncol = 2, byrow = TRUE)
    #     })
    #   } else {
    #     if (length(idx) != length(dist_to_ls)){stop("nope")}
    #     ls_endpts <- lapply(1L:filter_dim, function(d_i){        ## Construct the level sets
    #       ## Choose the interval (half) width to plot 
    #       if (idx[d_i] == 0){
    #         eps <- (base_interval_length/2)
    #       } else {
    #         eps <- (base_interval_length/2) + dist_to_ls[[d_i]]$target_dist[dist_order[[d_i]]][idx[d_i]] ## parameterized overlap
    #       }
    #       tmp <- as.vector(sapply(0L:(number_intervals[d_i] - 1L), function(idx){
    #         centroid <- filter_min[d_i] + (as.integer(idx)*base_interval_length[d_i]) + base_interval_length[d_i]/2.0
    #         c(centroid - eps[d_i], centroid + eps[d_i])
    #       }))
    #       matrix(tmp, ncol = 2, byrow = TRUE)
    #     })
    #   }
    #   
    #   plot(filter_values, pch = 20, 
    #        xlim = get_range(filter_values[, 1L]), ylim = get_range(filter_values[, 2L]),
    #        xlab = "", ylab = "") # xaxt = "n", yaxt = "n"
    #   if (pt_ids){ text(filter_values, labels = 1:nrow(filter_values), pos = 3) }
    #   # 
    #   # abline(h = 0, col = "gray", lty = 3, lwd = 1.5)
    #   binned_color <- rev(rainbow(nrow(cart_prod), start = 0, end = 4/6))
    #   cc <- 0L
    #   for (i in 1:nrow(cart_prod)){
    #     { ii <- cart_prod[i,1]; jj <- cart_prod[i,2] }
    #     { cls1 <- ls_endpts[[1]]; cls2 <- ls_endpts[[2]] }
    #     col <- binned_color[ii*jj]
    #     make_square(cls1[ii,], cls2[jj,], col = col)
    #     points(x = mean(cls1[ii,]), y = mean(cls2[jj,]), col = col, pch = 3)
    #     text(x = mean(cls1[ii,]), y = mean(cls2[jj,]), labels = trimws(paste0(c(ii, jj) - 1L, collapse = ",")), cex = 0.50, pos = 3)
    #   }
    #   for (d_i in 1:filter_dim){
    #     # browser()
    #     if (idx[[d_i]] > 0){
    #       pt <- pt_idx[[d_i]][[idx[[d_i]]]]
    #       points(x = filter_values[pt,1], y = filter_values[pt,2], col = "purple")
    #     }
    #   }
    #   if (!missing(prev)){
    #     # browser()
    #     ls_changed <- ms_mapper$get_ls_that_change(prev, idx)+1
    #     if (nrow(ls_changed) > 0){
    #       apply(ls_changed, 1, function(x) {
    #         ls1 <- lsmi[x[1], ]
    #         ls2 <- lsmi[x[2], ]
    #         rect(xleft = ls_endpts[[1]][ls1[1],1], xright = ls_endpts[[1]][ls1[1],2], 
    #              ybottom = ls_endpts[[2]][ls1[2],1], ytop = ls_endpts[[2]][ls1[2],2], 
    #              col = rgb(0.8, 0, 0.8, alpha = 0.3))
    #         rect(xleft = ls_endpts[[1]][ls2[1],1], xright = ls_endpts[[1]][ls2[1],2], 
    #              ybottom = ls_endpts[[2]][ls2[2],1], ytop = ls_endpts[[2]][ls2[2],2], 
    #              col = rgb(0.8, 0, 0.8, alpha = 0.3))
    #       })
    #     }
    #   }
    # } ## plot_2d
    
    # plot_segments <- function(d_i, fix_y = TRUE, fixed = max(filter_values[, d_i])){
    #   cls <- ls_endpts[[d_i]]
    #   for (i in 1:nrow(cls)){
    #     ls <- cls[i, ]
    #     if (i %in% c(1, nrow(cls))){
    #       if (i == 1){
    #         if (fix_y){ params <- list() }
    #         lines(x = c(ls[1], cls[i+1L, 1L]), y = rep(0.51, 2))
    #         text(x = mean(c(ls[1], cls[i+1L, 1L])), y = 0.5, pos = 3, labels = as.character(cc))
    #         cc <- cc + 1
    #         lines(x = c(cls[i+1L, 1L], ls[2]), y = rep(0.51, 2))
    #         text(x = mean(c(cls[i+1L, 1L], ls[2])), y = 0.5, pos = 3, labels = as.character(cc))
    #         cc <- cc + 1
    #       }
    #       else if (i == nrow(cls)){
    #         lines(x = c(cls[i-1L, 2L], ls[2]), y = rep(0.51, 2))
    #         text(x = mean(c(cls[i-1L, 2L], ls[2])), y = 0.5, pos = 3, labels = as.character(cc))
    #         cc <- cc + 1
    #       }
    #     } else {
    #       lines(x = c(cls[i-1L, 2L], cls[i+1L, 1L]), y = rep(0.51, 2))
    #       text(x = mean(c(cls[i-1L, 2L], cls[i+1L, 1L])), y = 0.5, pos = 3, labels = as.character(cc))
    #       cc <- cc + 1
    #       lines(x = c(cls[i+1L, 1L], ls[2]), y = rep(0.51, 2))
    #       text(x = mean(c(cls[i+1L, 1L], ls[2])), y = 0.5, pos = 3, labels = as.character(cc))
    #       cc <- cc + 1
    #     }
    #   }
    # }
 
#       lines(x = c(ls[1], ls[2]), y = c(-0.5, -0.5), col = binned_color[i])
#       lines(x = c(ls[2], ls[2]), y = c(-0.5, 0.5), col = binned_color[i])
#       lines(x = c(ls[2], ls[1]), y = c(0.5, 0.5), col = binned_color[i])
#       lines(x = c(ls[1], ls[1]), y = c(0.5, -0.5), col = binned_color[i])
#       text(labels = as.character(i - 1L), x = mean(c(ls[1], ls[2])), y = -0.5, col = binned_color[i], pos = 1)
#       points(x = mean(c(ls[1], ls[2])), y = 0, pch = 3, col = binned_color[i])
#       i <- i + 1
# 
#     points(filter_values[, d_i], rep(0, 10), pch = 20)
#     points(filter_values[pt_idx[[d_i]][idx], d_i], 0, pch = 21, cex = 1.5, col = "purple")
#     text(filter_values[, d_i], 0, labels = 1:10, pos = 3)
#     
#     animation::saveGIF({
#       for (i in 1:length(pt_idx[[1]])){ plot_configuration(1L, i) }
#     }, movie.name = "expanding_boxes.gif", interval = 0.2)
#     
#     plot_configuration(1L, 0)
#     plot_configuration <- function(d_i, idx){
#       if (idx == 0){
#         eps <- (base_interval_length/2)
#       } else {
#         ## Choose the interval (half) width to plot 
#         eps <- (base_interval_length/2) + dist_to_ls[[d_i]]$target_dist[dist_order[[d_i]]][idx] ## parameterized overlap
#       }
#     
#       ## Construct the level sets
#       ls_endpts <- lapply(1L:filter_dim, function(d_i){
#         tmp <- as.vector(sapply(0L:(number_intervals[d_i] - 1L), function(idx){
#           centroid <- filter_min[d_i] + (as.integer(idx)*base_interval_length[d_i]) + base_interval_length[d_i]/2.0
#           c(centroid - eps[d_i], centroid + eps[d_i])
#         }))
#         matrix(tmp, ncol = 2, byrow = TRUE)
#       })
#       
#       plot(filter_values[, d_i], rep(0, 10), pch = 20, ylim = c(-0.75, 0.75), xlab = "", ylab = "", yaxt = "n")
#       text(filter_values[, d_i], 0, labels = 1:10, pos = 3)
#       abline(h = 0, col = "gray", lty = 3, lwd = 1.5)
#       binned_color <- rev(rainbow(nrow(ls_endpts[[d_i]]), start = 0, end = 4/6))
#       cc <- 0L
#       cls <- ls_endpts[[d_i]]
#       for (i in 1:nrow(cls)){
#         ls <- cls[i, ]
#         if (i %in% c(1, nrow(cls))){
#           if (i == 1){
#             lines(x = c(ls[1], cls[i+1L, 1L]), y = rep(0.51, 2))
#             text(x = mean(c(ls[1], cls[i+1L, 1L])), y = 0.5, pos = 3, labels = as.character(cc))
#             cc <- cc + 1
#             lines(x = c(cls[i+1L, 1L], ls[2]), y = rep(0.51, 2))
#             text(x = mean(c(cls[i+1L, 1L], ls[2])), y = 0.5, pos = 3, labels = as.character(cc))
#             cc <- cc + 1
#           }
#           else if (i == nrow(cls)){
#             lines(x = c(cls[i-1L, 2L], ls[2]), y = rep(0.51, 2))
#             text(x = mean(c(cls[i-1L, 2L], ls[2])), y = 0.5, pos = 3, labels = as.character(cc))
#             cc <- cc + 1
#           }
#         } else {
#           lines(x = c(cls[i-1L, 2L], cls[i+1L, 1L]), y = rep(0.51, 2))
#           text(x = mean(c(cls[i-1L, 2L], cls[i+1L, 1L])), y = 0.5, pos = 3, labels = as.character(cc))
#           cc <- cc + 1
#           lines(x = c(cls[i+1L, 1L], ls[2]), y = rep(0.51, 2))
#           text(x = mean(c(cls[i+1L, 1L], ls[2])), y = 0.5, pos = 3, labels = as.character(cc))
#           cc <- cc + 1
#         }
#         
#         lines(x = c(ls[1], ls[2]), y = c(-0.5, -0.5), col = binned_color[i])
#         lines(x = c(ls[2], ls[2]), y = c(-0.5, 0.5), col = binned_color[i])
#         lines(x = c(ls[2], ls[1]), y = c(0.5, 0.5), col = binned_color[i])
#         lines(x = c(ls[1], ls[1]), y = c(0.5, -0.5), col = binned_color[i])
#         text(labels = as.character(i - 1L), x = mean(c(ls[1], ls[2])), y = -0.5, col = binned_color[i], pos = 1)
#         points(x = mean(c(ls[1], ls[2])), y = 0, pch = 3, col = binned_color[i])
#         i <- i + 1
#       }
#       text(x = cls[nrow(cls), 2L], y = 0.5, pos = 3, labels = as.character(cc))
#       points(filter_values[, d_i], rep(0, nrow(filter_values)), pch = 20)
#       if (idx > 0){
#         points(filter_values[pt_idx[[d_i]][idx], d_i], 0, pch = 21, cex = 1.5, col = "purple")
#       }
#       text(filter_values[, d_i], 0, labels = 1:nrow(filter_values), pos = 3)
#     }
#    
#     # abline(v = (base_interval_length[1]*0:number_intervals[1]) + filter_min[1], col = "red")
#     
#     # x_minmax <- c(unique(ls_bnds[, 1]) - .Machine$double.eps, unique(ls_bnds[, 3]) + .Machine$double.eps)
#     # 
#     # level_set_query_indices <-  tapply(0:((number_intervals[1]*2L) - 1L), rep(1:number_intervals[1], each = 2), identity, simplify = FALSE)
#     # 
#     #     ## Swaps the upper integer at index i with the lower integer at index j
#     #     swap_ul <- function(i, j){
#     #       upper <- level_set_query_indices[[i]][2]
#     #       lower <- level_set_query_indices[[j]][1]
#     #       level_set_query_indices[[i]][2] <- lower
#     #       level_set_query_indices[[j]][1] <- upper
#     #       return(level_set_query_indices)
#     #     }
#     #     
#     #     number_intervals[1]
#     #     
#     
#     stree <- Mapper::segment_tree(sort(as.vector(ls_bnds)))
#     pts <- as.vector(filter_values)
#     stree$insert_points(pts)
#     # stree_ptr <- stree$as_XPtr()
#     
#     ms_mapper$query(c(0L), c(1L))
#     
#     test_cover_compute <- sapply(0:6, function(i){
#       ms_mapper$compute_ls_idx(i)
#     })
#     
#     
#     # ms_mapper$build_segment_trees()
#     #ms_mapper$set_ls_idx(as.vector(unlist(level_set_query_indices)))
#     #ms_mapper$build_segment_trees()
#     #ms_mapper$get_parameterization(30L)
#     ms_mapper$get_ls_idx()
#     ms_mapper$swap_ls(1)
#     ms_mapper$get_ls_idx()
#     ms_mapper$swap_ls(2)
#     ms_mapper$get_ls_idx()
#     ms_mapper$swap_ls(3)
#     ms_mapper$get_ls_idx()
#     ms_mapper$swap_ls(4)
#     ms_mapper$get_ls_idx()
#     
#     stree$n
#     
#     ms_mapper$get_parameterization(9L)
#     
#     make_ls_idx <- function(i, n){
#       #  res <- vector(mode = "integer", length = n)
#       res <- matrix(0, nrow = n/2L, ncol = n)
#       max_c <- (n - 2L)/2L
#       hc <- (n - 2L)/2L 
#       lc <- 0L
#       for (j in 0:max_c){
#         if (j %% 2 == 0){
#           if (lc == 0L){ res[, j] <- rep(j, n/2L) }
#           else {
#             res[, j+1] <- c(seq(j, j-lc, length.out=lc), rep(j, (n/2L) - lc))
#           }
#           lc <- lc + 1L
#         }
# 
#         hc <- hc - 1L
#       }
#       return(res)   
#     }
#     
#     x_minmax <- c(unique(ls_bnds[, 1]), unique(ls_bnds[, 3]))
#     stree <- Mapper::segment_tree(intervals = x_minmax)
#     stree$insert_points(filter_values[, 1])
#     
#     
#     res_idx <- lapply(level_set_query_indices, function(ls_idx){
#       stree$queryInterval(ls_idx[1], ls_idx[2]+1L)
#     })
#     idx <- vector(mode = "integer", length = 10)
#     
#     m <- MapperRef$new(noisy_circle)$
#       set_cover(filter_values = f_x, type = "fixed rectangular", number_intervals = 5L, percent_overlap = 0)
#     
#     
#       compute_vertices(num_bins = 10)$
#       compute_edges()
#     
#     1L:number_intervals[1]
#     dist_to_lower_ls[1,]
#     as.integer(A[1,]) - 1L
#     
#     ## Testing out swapping 
#     X <- c(x_1=0.2, x_2=4/3, x_3=1.6)
#     endpts <- c(0, 1, 1, 2)
#     stree <- Mapper::segment_tree(endpts)
#     stree$insert_points(X)
#     stree$queryInterval(4-4, 6-4)
#     stree$queryInterval(4-4, 7-4) 
#     
#     m$cover$number_intervals
#     set_clustering_algorithm(cl = "single")$
#     set_distance_measure(measure = "euclidean")$
#     compute_vertices(num_bins = 10)$
#     compute_edges()
#   }
# )
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

## More preprocessing...
rowmatch <- function(A,B) { 
  # Rows in A that match the rows in B
  f <- function(...) paste(..., sep=":")
  if(!is.matrix(B)) B <- matrix(B, 1, length(B))
  a <- do.call("f", as.data.frame(A))
  b <- do.call("f", as.data.frame(B))
  match(a, b)
}
