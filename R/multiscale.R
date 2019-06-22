#' Multiscale Mapper 
#' @description Computes the Multiscale Mapper construction. See details. 
#' @param m The initial mapper to start with.
#' @param max_dim The maximum dimension of the simplices to consider.  
#' @param max_overlap The maximum overlap to consider. The multiscale construction will build only up to this value. 
#' @param time The notion of 'time' to measure birth/death times with. See details.
#' @param f Function to evaluate at evry index. The mapper instance is passed as the first argument.
#' @param stats Collect statistics on number of operations required to compute the structure? Defaults to false. 
#' @details This method constructs a sequence of Mappers connected by simplicial maps, which can be readily 
#' passed to e.g. the \code{Simpers} library to produce persistence diagrams. 
#' Because constructing a sequence of Mappers is a computationally expensive operation, a special cover
#' indexing structure is used to express the computation in terms of elementary simplicial maps. This method requires 
#' a fixed interval cover with a constant number of open sets, and uses a specific clustering algorithm. \cr
#' \cr
#' The \code{time} parameter may be specified as "integer", "measure", or (average) "overlap". 
#' Due to the restrictions on the cover this method assumes, these are all monotonically related to each other, 
#' however they may produce visually different persistence diagrams. 
multiscale <- function(m, max_dim=1L, max_overlap = 50, time = c("integer", "measure", "overlap"), f = NULL, stats = FALSE){
  stopifnot(m$cover$typename == "Fixed Interval")## TODO: add maximum diameter condition
  stopifnot(!is.null(m$filter))
  m$cover$validate(m$filter)
  
  ## Get filter values
  fv <- m$filter()
  
  ## Get initial eps parameter to minimize bias on the clusters
  eps_vals <- sapply(m$cover$index_set, function(alpha){
    preimage <- m$cover$construct_cover(m$filter, alpha)
    if (length(preimage) < 3){ return(0) }
    hcl <- hclust(dist(m$X(preimage)), method = "single")
    Mapper::cutoff_first_threshold(hcl)
  })
  global_eps <- mean(Filter(function(x){ x != 0 }, eps_vals))
  
  ## Setup initial clustering 
  m$simplicial_complex$id_policy <- "compressed"
  m$use_clustering_algorithm(cl = stable_clustering(global_eps), run_internal = TRUE)
  m$construct_pullback()  
  initial_complex <- m$simplicial_complex$serialize()
  m$simplicial_complex$id_policy <- "unique"  ## All newly generated ids from eps -> infty should be unique across map
  environment(m$clustering_algorithm)[["setup"]] <- TRUE ## Allow simplicial maps to be captured
  
  ## Create indexed cover 
  indexed_cover <- construct_indexed_cover(m, m$cover, ensure_unique = TRUE)
  # ls1 <- lapply(seq(0, prod(m$cover$number_intervals)-1), function(i) indexed_cover$extract_level_set(i))
  # ls2 <- m$cover$level_sets
  
  ## Retrieve all distinct (interpolated) eps values
  distinct <- Mapper:::nondecreasing_seq(indexed_cover$eps) # + .Machine$double.eps
  
  ## Reduce to only indices that have overlap < max_overlap
  distinct <- local({
    base_eps <- diff(apply(fv, 2, range))/m$cover$number_intervals
    distinct_g <- matrix(apply(distinct$eps, 1, function(eps) 1 - (base_eps / eps)), nrow = ncol(fv)) 
    reduced_idx <- apply(distinct_g, 2, function(g) all((g*100) < max_overlap))
    within(distinct, {
      idx <- idx[reduced_idx,,drop=FALSE]
      eps <- eps[reduced_idx,,drop=FALSE]
    })
  })
  
  ## Fixed for multidimensional start
  if (ncol(distinct$eps) > 1){
    d <- ncol(distinct$eps)
    fix <- vector(mode = "list", length=d)
    c_idx <- rep(-1, d)
    c_eps <- rep(0, d)
    cc <- 1L
    for (ii in head(order(distinct$eps[1,]), d-1)){
      c_idx[ii] <- 0
      c_eps[ii] <- distinct$eps[1,ii]
      fix[[cc]] <- list(idx=c_idx, eps=c_eps)
      cc <- cc + 1L
    }
    distinct$idx <- rbind(do.call(rbind, lapply(seq(d-1), function(ii){ fix[[ii]]$idx })), distinct$idx)
    distinct$eps <- rbind(do.call(rbind, lapply(seq(d-1), function(ii){ fix[[ii]]$eps })), distinct$eps)
  }

  ## plot the cover
  # if (ncol(m$cover$filter_values) == 1){
  #   stripchart(m$cover$filter_values)
  #   text(x = m$cover$filter_values, y = 1, labels = seq(nrow(m$cover$filter_values))-1L, pos = 3)
  #   abline(v = m$cover$interval_bounds()[,1], col = "blue")
  #   abline(v = m$cover$interval_bounds()[,2], col = "red")
  # } else if (ncol(m$cover$filter_values) == 2){
  #   plot(m$cover$filter_values)
  #   text(x = m$cover$filter_values, labels = seq(nrow(m$cover$filter_values))-1L, pos = 3)
  #   abline(v = m$cover$interval_bounds()[,1], col = "blue")
  #   abline(h = m$cover$interval_bounds()[,2], col = "blue")
  #   abline(v = m$cover$interval_bounds()[,3], col = "red")
  #   abline(h = m$cover$interval_bounds()[,4], col = "red")
  # }
  if (stats){
    number_pullbacks <- 0L
    number_edges 
  }
  
  pb <- txtProgressBar(min = 0, max = nrow(distinct$idx), style = 3)
  for (i in seq(nrow(distinct$idx))){
    # browser()
    tidx <- distinct$idx[i,]
    changes <- indexed_cover$update_index(tidx, max_dim)
    
    ## To save adjacent pullback indices 
    # pid_comb_update <- vector(mode = "list", length(changes$indices_changed))
    pid_comb_update <- vector(mode="list", length = max_dim)
    cc <- 1L
    
    ## Update the individual pullback(s) that changed
    for (a_i in changes$indices_changed){
      pid_to_update <- m$cover$index_set[a_i + 1L]
      new_ls <- indexed_cover$extract_level_set(a_i)+1L
      new_idx <- changes$points_updated+1L
      if (length(new_idx) != 1){ stop("added more than 1 point.") }
      
      ## Get the pullback ids of sets intersecting the pullback that just changed. 
      ## This is equivalent to getting the open sets intersecting the cofaces of the vertices in the current pullback to update.
      adj_pids <- Mapper:::connected_pullbacks(pid_to_update, m$pullback, m$simplicial_complex$as_XPtr())
      if (length(adj_pids) > 0){
        for (k in 1L:max_dim){
          pid_combs <- append(list(pid_to_update), lapply(1L:k, function(i){ adj_pids[[pid_to_update]] }))
          pid_comb_update[[k]] <- rbind(pid_comb_update[[k]], as.matrix(do.call(expand.grid, pid_combs)))
        }
      }
      
      ## Run the clustering algorithm. This updates the points in the vertices, potentially changes complex,
      ## and (if the complex changes) records elementary simplicial maps
      # message(sprintf("Updating pullback id %s w/ point %d", pid_to_update, new_idx))
      m$clustering_algorithm(pid_to_update, new_ls, new_idx)
    }
    
    ## Update the nerve
    comb_len <- sapply(changes$pullback_indices_to_check, length)
    for (k in 1L:max_dim){
      ## Collect the pullback indices adjacent to sets that changed
      pid_labels <- lapply(which(comb_len == (k+1)), function(i){
        idx <- changes$pullback_indices_to_check[[i]]
        m$cover$index_set[idx+1L]
      })
      pid_labels_to_update <- unique(rbind(pid_comb_update[[k]], do.call(rbind, pid_labels)))
    
      # idx_pairs <- matrix(apply(changes$pullback_indices_to_check, 2, function(idx){ m$cover$index_set[idx+1L] }), ncol = 2)
      # idx_pairs <- unique(rbind(idx_pairs, do.call(rbind, pairs_to_update)))
      
      ## Update the nerve locally
      simplices <- m$construct_nerve(indices = pid_labels_to_update, k = k, modify = FALSE)
      c_env <- environment(m$clustering_algorithm)
      insert <- FALSE
      for (j in seq(nrow(simplices))){
        if (!m$simplicial_complex$find(simplices[j,])){
          m$simplicial_complex$insert(simplices[j,])
          c_env[["s_maps"]] <- c(c_env[["s_maps"]], sprintf("i %s", paste0(simplices[j,], collapse = " ")))
          insert <- TRUE
        }
      }
      if (insert){ c_env[["s_maps"]] <- c(c_env[["s_maps"]], sprintf("# %d", c_env[["c_time"]])) }
      }
    ## Optionally validate 
    # validate_inc(m, indexed_cover)
    
    ## Evaluate the user-function, if specified 
    if (!missing(f) && is.function(f)){ f(m) }
    
    ## Update progress 
    setTxtProgressBar(pb, value = i)
    c_env[["c_time"]] <- c_env[["c_time"]]+1L
  } # for (i in seq(nrow(distinct$idx)))
  close(pb)
  
  # What unit to consider as 'time' in the diagram? 
  if (missing(time) || time == "integer"){
    time <- as.integer(seq(nrow(distinct$eps)))
  } else if (time == "measure"){
    # all(order(apply(distinct$eps, 1, prod)) == seq(nrow(distinct$eps)))
    time <- as.numeric(apply(distinct$eps, 1, prod))
  } else if (time == "overlap"){
    base_eps <- diff(apply(fv, 2, range))/m$cover$number_intervals
    time <- as.numeric(colMeans(matrix(apply(distinct$eps, 1, function(eps) 1 - (base_eps / eps)), nrow = ncol(fv))))
  } else {
    stop(sprintf("Unknown unit name %s", time))
  }
  stopifnot(length(time) == nrow(distinct$idx))
  
  ## Map the timings to the specified unit
  s_maps <- environment(m$clustering_algorithm)[["s_maps"]]
  bd_idx <- grep(x = s_maps, pattern = "#\\s*(\\d+)", value = FALSE)
  int_times <- as.integer(gsub(x = s_maps[bd_idx], pattern = "[#]\\s*(\\d+)", replacement = "\\1"))
  s_maps[bd_idx] <- sprintf("# %f", time[int_times])
  
  ## Aggregate results
  initial <- simplextree::simplex_tree()
  initial$deserialize(initial_complex)
  result <- list(initial_mapper = initial, simplicial_maps = s_maps, cluster_f = basic_cluster(global_eps))
  return(result)
}

## Clustering algorithm that also creates a set of simplicial maps
## Need a 'stable' clustering algorithm, i.e. an algorithm which lends itself to creating elementary insertions 
## or collapses if a single point is inserted into a prior clustering
stable_clustering <- function(g_eps){
  requireNamespace("parallelDist")
  
  ## Simpers requires the vertex indices be consistent across all simplicial maps. This implies that if new 
  ## vertices are requested, they must not intersect the names of vertices that e.g. were removed by a collapse.
  setup <- FALSE
  s_maps <- c()
  c_time <- 1L
  c_eps <- list()
  
  ## Given a set of disjoint vertices and associated points, assigns each point into 
  ## a connected component based on which vertex it falls in
  cc <- function(vids, vertices, pt_ids){
    c_cc <- vector(mode="integer", length(pt_ids))
    c_check <- vector(mode="integer", length(pt_ids))
    for (vid in vids){
      idx <- match(intersect(vertices[[as.character(vid)]], pt_ids), pt_ids) 
      if (length(idx)){ 
        c_cc[idx] <- as.integer(vid) 
        c_check[idx] <- 1L
      }
    }
    ## Map points to vertex ids
    # c_cc <- sapply(pt_ids, function(pt_id) { 
    #   cvid <- sapply(vertices[vids], function(pts){ pt_id %in% pts })
    #   if (length(cvid) !- 1){ stop("Invalid vertex mapping.") }
    #   return(vids[cvid])
    # })
    if (any(c_check == 0)){ browser() }
    return(c_cc)
  }
  ## Retrieves the epsilon value for a given set of points 'x'
  get_eps <- function(x) { 
    return(g_eps)
    # if (nrow(x) <= 1){ return(0) }
    # if (nrow(x) == 2){ return(dist(matrix(x))) }
    # hcl <- hclust(parallelDist::parallelDist(x), method = "single")
    # return(Mapper::cutoff_first_threshold(hcl))
  }
  
  cluster <- function(pid, idx, new_idx=NULL){
    ## Setup chunk. This is run once per open set. 
    if (!setup){
      if (0 %in% idx){ stop("0-based indices given. Expects 1-based.") }
      if (is.null(idx) || length(idx) == 0){ return(NULL) }
      if (length(idx) <= 2L){ return(rep(1L, length(idx))); }
      base_hcl <- hclust(parallelDist::parallelDist(self$X(idx)), method = "single")
      c_eps[[pid]] <<- get_eps(self$X(idx))
      return(cutree(base_hcl, h = c_eps[[pid]]))
    } else {
      ## Run cutting heuristic to get new epsilon-value
      c_eps[[pid]] <<- max(c_eps[[pid]], get_eps(self$X(idx)))
      
      ## Eps-neighborhood check: where in the metric space does this new point fall? 
      p_idx <- setdiff(idx, new_idx)
      nn_eps <- (if (length(p_idx) == 0){ NULL } 
                 else { RANN::nn2(data=self$X(p_idx), query=self$X(new_idx), searchtype="radius", radius=c_eps[[pid]]) })
      nn_idx <- (if (length(p_idx) == 0){ numeric(0) } else { nn_eps$nn.idx[nn_eps$nn.idx != 0] })
      
      ## If the point is not in the eps-neighborhood of any CC, it is a new cluster, perform elementary inclusion
      if (length(nn_idx) == 0){
        ## BEGIN ELEMENTARY INCLUSION
        new_cl <- self$simplicial_complex$generate_ids(1)
        if (as.character(new_cl) %in% names(self$vertices)){ stop("Invalid vertix id generated.") }
        
        ## Modify the vertices, the pullback, and the simplicial complex
        self$vertices[[as.character(new_cl)]] <- new_idx
        self$pullback[[pid]] <- unique(c(self$pullback[[pid]], new_cl))
        self$simplicial_complex$insert(new_cl)
        
        ## Record the elementary inclusion
        s_maps <<- c(s_maps, sprintf("i %d", new_cl))
        s_maps <<- c(s_maps, sprintf("# %d", c_time))
        # message(tail(s_maps, 1))
        ## END ELEMENTARY INCLUSION
      } else {
        ## Vertex ids in the pullback
        pvids <- as.character(self$pullback[[pid]])
        pt_cc <- cc(vids = pvids, vertices = self$vertices, pt_ids = p_idx[nn_idx])
        
        ## BEGIN VERTEX INSERTION
        ## If the point is in the eps-neighborhood of only 1 CC, no map is needed. 
        ## For mapper, this reduces to a point insertion into one of the vertices
        if (length(unique(pt_cc)) == 1L){
          c_vid <- as.character(unique(pt_cc))
          self$vertices[[c_vid]] <- unname(c(self$vertices[[c_vid]], new_idx))
          return(NULL) ## Don't need to update pullback or simplicial complex
        }
        ## END VERTEX INSERTION
        
        ## BEGIN ELEMENTARY COLLAPSE
        ## Otherwise, the point is in the eps-neighborhood of 2 or more CCs. Collapse the 
        ## vertices to a single vertex. 
        intersecting_vids <- as.character(unique(pt_cc))
        target_vid <- head(intersecting_vids, 1)
        collapse_vids <- setdiff(intersecting_vids, target_vid)
        
        ## Remove collapsed vertices, update target vertex, including new point
        tmp_pts <- as.vector(unlist(self$vertices[collapse_vids]))
        self$vertices <- self$vertices[Filter(function(x) !x %in% collapse_vids, names(self$vertices))]
        self$vertices[[target_vid]] <- unique(c(self$vertices[[target_vid]], tmp_pts, new_idx))
        
        ## Remove collapsed vertices from pullback
        self$pullback <- lapply(self$pullback, function(pids){ pids[!pids %in% collapse_vids] })
        
        ## Perform the collapse on the complex
        for (c_id in collapse_vids){
          self$simplicial_complex$collapse(as.integer(c_id), as.integer(target_vid), as.integer(target_vid))
        }
        
        ## Record the elementary collapses
        s_maps <<- c(s_maps, sprintf("c %s t %s", paste0(collapse_vids, collapse = " "), target_vid))
        s_maps <<- c(s_maps, sprintf("# %d", c_time))
        # message(tail(s_maps, 1))
        ## END ELEMENTARY COLLAPSE
      } 
    }
  } # cluster function
  return(cluster)
}

# clustering function based on global cut threshold
basic_cluster <- function(g_eps){
  function(pid, idx){
    if (0 %in% idx){ stop("0-based indices given. Expects 1-based.") }
    if (is.null(idx) || length(idx) == 0){ return(NULL) }
    if (length(idx) <= 2L){ return(rep(1L, length(idx))); }
    base_hcl <- hclust(parallelDist::parallelDist(self$X(idx)), method = "single")
    cutree(base_hcl, h = g_eps)
  }
}


construct_indexed_cover <- function(m, cover, ensure_unique = FALSE){
  ## Locals
  fv <- m$filter()
  n <- nrow(fv)
  d <- ncol(fv)
  base_interval_length <- diff(apply(fv, 2, range))/cover$number_intervals
  
  ## R^(k x 2d) matrix of interval lower/upper bounds in form of (x_lb, y_lb, ..., x_ub, y_ub, ...)
  interval_bnds <- cover$interval_bounds(m$filter)
  
  ## Encode cover with multi-index
  key_to_ind <- function(idx_set) { lapply(strsplit(substr(idx_set, 2, nchar(idx_set)-1), " "), as.integer) }
  cover_multi_idx <- do.call(rbind, key_to_ind(cover$index_set))
  
  ## Encode points by both a flat and multi index. 
  pt_flat_idx <- Mapper:::constructLevelSetIndex(fv, interval_bnds)
  pt_multi_idx <- cover_multi_idx[pt_flat_idx,,drop=FALSE]
  
  ## Calculate point-to-set distances to non-intersecting open sets
  {
    ## Get upper and lower bounds of the set per point
    lb <- interval_bnds[pt_flat_idx, seq(d), drop=FALSE]
    ub <- interval_bnds[pt_flat_idx, seq(d+1, 2*d), drop=FALSE]
    dist_to_lb <- matrix(sapply(seq(n), function(i){ abs(fv[i,] - lb[i,]) }), nrow = d)
    dist_to_ub <- matrix(sapply(seq(n), function(i){ abs(fv[i,] - ub[i,]) }), nrow = d)
    
    ## Use these distances to compute the distances to 'target' level sets in each dimension
    interval_params <- lapply(seq(d), function(d_i){
      seq_int <- seq(cover$number_intervals[d_i])
      tmp <- mapply(function(i, relative_idx){
        target_pos <- setdiff(seq_int, relative_idx)
        offset <- (base_interval_length[d_i]*(abs(relative_idx - target_pos)-1L))
        eps <- base_interval_length[d_i] + 
          ifelse(relative_idx < target_pos, offset + dist_to_ub[d_i, i], offset + dist_to_lb[d_i, i])*2
        list(target_pos, eps)
      }, seq(n), pt_multi_idx[,d_i])
      list(target_pos = do.call(rbind, tmp[1,]), target_eps = do.call(rbind, tmp[2,]))
    })
    eps_order <- lapply(seq(d), function(d_i) order(interval_params[[d_i]]$target_eps))
    interval_sizes <- lapply(seq(d), function(d_i) sort(interval_params[[d_i]]$target_eps))
    
    ## If requested, make the interval sizes unique. This uses a custom comparator to determine equality. 
    m_eps <- sqrt(.Machine$double.eps)
    # dupped <- function(x, n = length(x)){ abs(x[1:(n-1)] - x[2:n]) < m_eps }
    if (ensure_unique){
      for (d_i in seq(d)){
        while(any(duplicated(interval_sizes[[d_i]]))){
          dup_idx <- which(duplicated(interval_sizes[[d_i]]))
          noise <- m_eps
          interval_sizes[[d_i]][dup_idx] <- interval_sizes[[d_i]][dup_idx] + noise
        }
        stopifnot(all(order(interval_sizes[[d_i]]) == seq(length(interval_sizes[[d_i]]))))
      }
    }
  } 
  
  ## Extract ordered paths 
  {
    ## Given matrices (x,y), attempts to find the index matching each row in x to a row in y
    rowmatch <- function(x, y){ match(data.frame(t(x)), data.frame(t(y))) }
    
    ## Generate the index paths
    idx_paths <- vector(mode = "list", length = d) ## contains the starting index the path take into unique paths
    uniq_idx_paths <- vector(mode = "list", length = d)
    for (d_i in seq(d)){
      ordered_path <- function(i) { 
        target_idx <- with(interval_params[[d_i]], { target_pos[i, order(target_eps[i,])] })
        c(pt_multi_idx[i, d_i], target_idx)
      }
      pt_idx_paths <- do.call(rbind, lapply(seq(n), ordered_path))
      uniq_idx_paths[[d_i]] <- unique(pt_idx_paths)
      idx_paths[[d_i]] <- rowmatch(pt_idx_paths, uniq_idx_paths[[d_i]]) ## use includes
    }
  }
  
  ## Extract canonical cut indices
  { 
    ## Partitions the interval sizes into which canonical cover each belongs too
    canonical_cuts <- lapply(1:d, function(d_i) {
      critical_dists <- ((base_interval_length[d_i]/2) * seq(2, cover$number_intervals[d_i] - 1L))*2
      canonical_idx <- findInterval(interval_sizes[[d_i]], vec = critical_dists)
      tmp <- rle(canonical_idx)[["lengths"]]
      head(c(0L, cumsum(tmp)), length(tmp))
    })
  }
  
  ## Create the indexing structure. The minimal information needed is as follows:
  ## 1) multi point index for each point
  ## 2) critical ordering (eps_order) + corresponding eps-parameters 
  ## 3) point path information + corresponding possible paths
  ## 4) canonical cover indexing information
  to_0_based <- function(x) { x-1L }
  initializer <- list(
    n = n, 
    number_intervals = as.vector(cover$number_intervals), 
    pt_multi_idx = pt_multi_idx-1L, 
    pt_path_idx = do.call(cbind, idx_paths)-1L, 
    target_order = lapply(eps_order, to_0_based), 
    eps_values = interval_sizes,
    unique_paths = lapply(uniq_idx_paths, to_0_based), 
    canonical_cuts = canonical_cuts ## already 0-based
  )
  ## Make the indexed cover
  indexed_cover <- Mapper:::MultiScale$new(initializer)
  return(indexed_cover)
}

## Multiscale Module
Rcpp::loadModule("multiscale_module", TRUE)
