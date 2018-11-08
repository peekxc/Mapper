
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
      
      ## Collect statistics if requested 
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
        ls_pair_idx <- apply(ls_pairs_to_update, 2, function(idx){
          self$cover$index_set[idx+1]
        })
        self$compute_edges(ls_pair_idx)
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

## More preprocessing...
rowmatch <- function(A,B) { 
  # Rows in A that match the rows in B
  f <- function(...) paste(..., sep=":")
  if(!is.matrix(B)) B <- matrix(B, 1, length(B))
  a <- do.call("f", as.data.frame(A))
  b <- do.call("f", as.data.frame(B))
  match(a, b)
}

## Multiscale Module
Rcpp::loadModule("multiscale_module", TRUE)
