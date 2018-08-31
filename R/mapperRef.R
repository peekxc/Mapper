#' @title Mapper Reference Class (R6) implementation
#' @docType class
#' @description Composes a set of classes to perform Mapper efficiently.
#' @author Matt Piekenbrock
#' @import methods
#' @export MapperRef
MapperRef <- R6Class("MapperRef", 
  private = list(.X=NA, .cover=NA, .clustering_algorithm=NA, .measure=NA, k_skeleton=NA, config=NA)
)

MapperRef$set("public", "initialize", function(X, cover){
  private$.X <- X
  self$cover <- cover
})

## The cover stores the filter values
MapperRef$set("active", "cover", 
  function(value){ #function(fv, type = c("restrained rectangular"), ...)
    if (missing(value)){
      private$.cover
    } else {
      stopifnot(inherits(value, "CoverRef"))
      private$.cover <- value
      invisible(self)
    }
  }
)
  # if (is.null(dim(fv)) && is.numeric(fv)){ fv <- matrix(fv, nrow = length(fv), ncol = 1) }
  # if (is.null(dim(fv)) || (!"matrix" %in% class(fv))) { stop("Filter values must be numeric and in matrix form") }
  # if (missing(type) || type == "restrained rectangular"){
  #   cover <<- r_rect_cover$new(filter_values = fv, type = "restrained rectangular", ...) # additional arguments to pass to the rectangular cover initialization method
  # } else if (type == "fixed rectangular"){
  #   cover <<- f_rect_cover$new(filter_values = fv, type = "fixed rectangular", ...) # additional arguments to pass to the rectangular cover initialization method
  # } else {
  #   stop(sprintf("%s cover not supported. Perhaps consider contributing it to the package w/ a pull request?", type))
  # }

MapperRef$set("active", "clustering_algorithm", 
  function(value){
    if (missing(value)){ private$.clustering_algorithm }
    else {
      stopifnot(is.function(value))
      private$.clustering_algorithm <- value
    }
  }
)

## Changes the clustering algorithm used by mapper.
## Must accept a 'dist' object and return a static or 'flat' clustering result
MapperRef$set("public", "set_clustering_algorithm", 
  function(cl = c("single", "ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid"), num_bins = 10L, ...){
    ## Use a linkage criterion + cutting rule
    if (missing(cl)) { cl <- "single" }
    if (class(cl) == "character"){
      hclust_opts <- c("single", "ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
      if (!cl %in% hclust_opts){ stop(sprint("Unknown linkage method passed. Please use one of (%s). See ?hclust for details.", paste0(hclust_opts, collapse = ", "))) }
      self$clustering_algorithm <- function(X, idx, num_bins = num_bins, ...){
        if (length(idx) <= 1){ return(1L); }
        dist_x <- parallelDist::parallelDist(X[idx,], method = "euclidean")
        hcl <- fastcluster::hclust(dist_x, method = cl)
        cutoff_first_bin(hcl, num_bins)
      }
    }
  }
)

## Sets the distance measure to use
## Supports any measure in proxy::pr_DB$get_entry_names()
MapperRef$set("public", "set_distance_measure", function(measure){
  available_measures <- toupper(proxy::pr_DB$get_entry_names())
  stopifnot(is.character(measure))
  stopifnot(toupper(measure) %in% available_measures)
  private$.measure <- measure
})


## Computes the distance matrix for each level set
MapperRef$set("public", "computeLevelSetDist", function(...){
  for (i in 1:length(self$cover$level_sets)){
    pt_idx <- self$cover$level_sets[[i]]$points_in_level_set
    if (length(pt_idx) <= 1){ self$cover$level_sets[[i]]$dist <- dist(0L) }
    else {
      if ("dist" %in% class(X)){
        self$cover$level_sets[[i]]$dist <- dist_subset(dist = self$X, idx = pt_idx)
      } else {
        if (is.character(measure)){
          self$cover$level_sets[[i]]$dist <- parallelDist::parallelDist(self$X[pt_idx,], method = measure)
        } else if (is.function(measure)){
          self$cover$level_sets[[i]]$dist <- measure(self$X[pt_idx,], ...)
        } else {
          stop("Unknown 'measure' argument. Must be either character string or function.")
        }
      }
    }
  }
  return(self)
})

## Computes the K-skeleton in a similar fashion as described by Singh et. al, section 3.2. Note that this
## procedure implicitly assumes the full abstract simplicial complex produced by Mapper can be treated as
## a flag complex.
MapperRef$set("public", "compute_k_skeleton", function(k, ...){
  stopifnot(k >= 0)
  if (k >= 0L){ self$computeNodes(...) }
  if (k >= 1L){ self$computeEdges(...) }
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
MapperRef$set("public", "compute_vertices", function(...){
  stopifnot(!is.na(self$cover$level_sets))
  
  ## Initialize the graph as an empty list
  self$k_skeleton <- list()

  ## Iterate through the 'dist' objects stored for each level set, perform the clustering
  cl_res <- lapply(self$cover$level_sets, function(ls) {
    if (length(ls$dist) > 0) clustering_algorithm(ls$dist, ...)
  })

  ## Precompute useful variables to know
  n_vertices <- sum(sapply(cl_res, function(cl) length(unique(cl))))
  vertice_idx <- unlist(mapply(function(cl, ls_i) if (length(cl) > 0) paste0(ls_i, ".", unique(cl)), cl_res, 1:length(cl_res)))
  n_lvlsets <- length(self$cover$level_sets)

  ## Agglomerate the nodes into a list. This matches up the original indexes of the filter values with the
  ## the clustering results, such that each node stores the original filter index values as well as the
  ## creating a correspondence between the node and it's corresponding level set flat index (lsfi)
  ## TODO: Cleanup and either vectorize or send down to C++
  self$k_skeleton$vertices <- vector(mode = "list", length = n_vertices)
  v_i <- 1L
  for (lsfi in 1:n_lvlsets){
    cl_i <- cl_res[[lsfi]]
    if (!is.null(cl_i)){
      ## Extract the vertex point indices for each cluster
      vertex_pt_idx <- lapply(unique(cl_i), function(cl_idx) self$cover$level_sets[[lsfi]]$points_in_level_set[which(cl_i == cl_idx)])
      for (vertex in vertex_pt_idx){
        attr(vertex, "level_set") <- lsfi
        if (any(is.na(vertex))){ browser() }
        self$k_skeleton$vertices[[v_i]] <- vertex
        v_i <- v_i + 1L
      }
    }
  }
})

## Computes the edges composing the topological graph (1-skeleton). 
## Assumes the nodes have been computed. 
MapperRef$set("public", "compute_edges", function(level_sets = NULL){

    ## Retrieve the level set flat indices (LSFI) for each corresponding node
    node_lsfi <- sapply(self$k_skeleton$vertices, function(vertex) attr(node, "level_set")) # which level set (by value) each node (by index) is in
  
    ## Create map from the level set flat index (by index) to the node indices the level set stores
    ## Note in this map empty level sets are NULL
    ls_node_map <- lapply(seq(length(self$cover$level_sets)), function(lvl_set_idx) {
      node_indices <- which(node_lsfi == lvl_set_idx)
      if (length(node_indices) == 0){ return(NULL) } else { return(node_indices) }
    })
    
    ## Retrieve the valid level set index pairs to compare. In the worst case, with no cover-specific optimization
    ## or 1-skeleton assumption, this may just be all pairwise combinations of LSFI's for the full simplicial complex.
    ## If the specific set of LSFI's were given
    if (missing(self$cover$level_sets) || is.null(self$cover$level_sets)){
      ## Let Rcpp handle the O(n^2) non-empty intersection checks
      ls_to_compare <- self$cover$valid_pairs()
      self$k_skeleton$adjacency <- adjacencyCpp(ls_pairs = ls_to_compare, nodes = self$k_skeleton$nodes, ls_node_map = ls_node_map);
    } else {
      self$k_skeleton$edgelist <- Mapper:::edgeList_int(ls_pairs = level_sets, nodes = self$k_skeleton$nodes, ls_node_map = ls_node_map)
    }
  }
)


## plotNetwork uses the 'network' package to plot the Mapper construction, w/ suitable defaults
## corresponding to what's commonly used in practice, all of which can be overridden.
# MapperRef$methods(plotNetwork = function(...){
# 
#   normalize <- function(x){ (x - min(x))/(max(x) - min(x)) }
#   ## Turn any given parameters into list
#   params <- list(...)
# 
#   ## If vertex color supplied, great! If not, use rainbow scale from blue --> red with values
#   ## corresponding to average filter values
#   if (is.null(params[["vertex.col"]])){
#     agg_pt_fv <- sapply(G$nodes, function(n_idx){ apply(matrix(cover$filter_values[n_idx,], nrow = length(n_idx)), 1, mean)})
#     agg_node_fv <- sapply(agg_pt_fv, mean)
#     rbw_pal <- rev(rainbow(100, start = 0, end = 4/6))
#     binned_idx <- cut(agg_node_fv, breaks = 100, labels = F)
#     params[["vertex.col"]] <- rbw_pal[binned_idx]
#   }
# 
#   ## If vertex size specified, great! If not, scale node size logarithmically with the number of points each contains
#   if (is.null(params[["vertex.cex"]])){
#     node_sizes <- sapply(G$nodes, function(n_idx){ length(n_idx) })
#     node_cex <- (6 - 1)*normalize(log(pmax(2, node_sizes)))+1 ## Scale node size logarithmically between [0.1, 4]
#     params[["vertex.cex"]] <- node_cex
#   }
# 
#   ## Other defaults
#   if (is.null(params[["label"]])){ params[["label"]] <- 1:length(G$nodes) }
#   if (is.null(params[["displaylabels"]])){ params[["displaylabels"]] <- TRUE }
# 
#   ## Construct the core network w/ network package
#   base_net <- network::network.initialize(nrow(G$adjacency), directed = F)
#   mapper_graph <- network::network.adjacency(x = G$adjacency, g = base_net)
#   params[["x"]] <- mapper_graph
# 
#   ## Evaluate in terms of passed in parameters
#   do.call(network::plot.network, params)
# })

## S3-like print override
MapperRef$set("public", "format", function(...){
  if ("dist" %in% class(self$X)){ n <- attr(self$X, "Size") } else { n <- nrow(self$X) }
  # if (!is.null(self$config$call))
  #   cat("\nCall:\n", deparse(self$config$call), "\n\n", sep = "")
  message <- sprintf("Mapper construction for %d objects", n)
  message <- append(message, format(self$cover))
  return(message)
})

only_combinations <- function(mat){
  mn <- pmin(mat[, 2], mat[, 1])
  mx <- pmax(mat[, 2], mat[, 1])
  int <- as.numeric(interaction(mn, mx))
  mat[match(unique(int), int),]
}

## Exports the internal mapper core structures to a TDAmapper output
# MapperRef$methods(exportTDAmapper = function(){
#   level_of_vertex <- sapply(self$k_skeleton$nodes, function(ni) attr(ni, "level_set"))
#   structure(
#     list(
#       adjacency = self$k_skeleton$adjacency,
#       num_vertices = length(self$k_skeleton$nodes),
#       level_of_vertex = level_of_vertex,
#       points_in_vertex = lapply(self$k_skeleton$nodes, as.vector),
#       points_in_level = unname(lapply(self$cover$level_sets, function(ls) ls$points_in_level_set)),
#       vertices_in_level = lapply(1:length(self$cover$level_sets), function(ls_idx) {
#         tmp <- which(level_of_vertex == ls_idx)
#         if (length(tmp) == 0){ return(-1) }
#         else return(tmp)
#       })),
#     class = "TDAmapper"
#   )
# })

## Exports the internal mapper core structures to a Mapper object, which is a list containing:
## 1. nodes member := list of integer vectors representing the indices of the original data the current node intersects with. Also
##                    contains attribute data storing the level set flat index of the level set the node is in.
## 2. adjacency := adjacency matrix of the resulting graph.
# MapperRef$methods(exportMapper = function(){
#   result <- .self$G
#   node_lsfi <- sapply(result$nodes, function(node) attr(node, "level_set"))
#   result$level_sets <- lapply(1:length(cover$level_sets), function(lsfi){ which(node_lsfi == lsfi) })
#   names(result$level_sets) <- apply(cover$index_set, 1, function(lsmi) paste0("(", paste(lsmi, collapse = ","), ")"))
#   cover_type <- paste0(toupper(substr(cover$type, start = 1, stop = 1)), tolower(substr(cover$type, start = 2, stop = nchar(cover$type))))
#   z_d <- ncol(cover$filter_values)
#   attr(result, ".summary") <- c(sprintf("Mapper object with filter function f: %s -> %s",
#                                         ifelse(is(X, "dist"), "dist(X)",
#                                                ifelse(ncol(X) > 1, sprintf("X^%d", ncol(X)), "X")),
#                                         ifelse(z_d > 1, sprintf("Z^%d", z_d), "Z")),
#                                 sprintf("Configured with a %s cover comprising %d open sets", cover$type, length(cover$level_sets)),
#                                 sprintf("The 1-skeleton contains %d nodes and %d edges", length(G$nodes), sum(G$adjacency == 1L)/2L))
#   class(result) <- "Mapper"
#   return(result)
# })



