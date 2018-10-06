#' @title Mapper Reference Class (R6) implementation
#' @docType class
#' @description Composes a set of classes to perform Mapper efficiently.
#' @author Matt Piekenbrock
#' @import methods
#' @export MapperRef
MapperRef <- R6Class("MapperRef", 
  private = list(.X=NA, .cover=NA, .clustering_algorithm=NA, .measure=NA, .simplicial_complex = NA, .vertices=list(), .cl_map=list(), .config=NA),
  lock_class = FALSE,  ## Feel free to add your own members
  lock_objects = FALSE ## Or change existing ones 
)

MapperRef$set("public", "initialize", function(X){
  stopifnot(is.numeric(X))
  if (is.null(dim(X)) && !is(X, "dist")){ X <- as.matrix(X) }
  private$.X <- X
  private$.simplicial_complex <- simplex_tree()
})

## The level_set -> vertex mapping is available to view, but read-only
MapperRef$set("active", "ls_vertex_map", 
  function(value){ 
    if (missing(value)){
      private$.cl_map
    } else {
      stop("'ls_vertex_map' is read-only.")
    }
  }
)

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

MapperRef$set("active", "vertices", 
  function(value){
    if(missing(value)){
      private$.vertices
    } else {
      stop("`$vertices` is read-only. To update the vertex membership, use the 'compute_vertices' function.", call. = FALSE)
    }
  }  
)

## Mapper stores the simplicial complex as a simplex tree.
MapperRef$set("active", "simplicial_complex", 
  function(value){
    if (missing(value)){ private$.simplicial_complex }
    else {
      stop("`$simplicial_complex` is read-only. To change the complex, use the objects (simplex tree) methods directly.", call. = FALSE)
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
MapperRef$set("public", "use_clustering_algorithm", 
  function(cl = c("single", "ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid"), num_bins = 10L){
    ## Use a linkage criterion + cutting rule
    if (missing(cl)) { cl <- "single" }
    if (class(cl) == "character"){
      hclust_opts <- c("single", "ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
      if (!cl %in% hclust_opts){ stop(sprint("Unknown linkage method passed. Please use one of (%s). See ?hclust for details.", paste0(hclust_opts, collapse = ", "))) }
      
      ## Closured representation to substitute default parameters. 
      create_cl <- function(cl, num_bins.default){
        function(X, idx, num_bins = num_bins.default){
          if (length(idx) <= 1){ return(1L); }
          dist_x <- parallelDist::parallelDist(X[idx,], method = "euclidean")
          hcl <- fastcluster::hclust(dist_x, method = cl)
          cutoff_first_bin(hcl, num_bins)
        }
      }
      self$clustering_algorithm <- create_cl(cl = cl, num_bins.default = force(num_bins))
    }
    invisible(self)
  }
)

## Sets the distance measure to use
## Supports any measure in proxy::pr_DB$get_entry_names()
MapperRef$set("public", "use_distance_measure", function(measure){
  available_measures <- toupper(proxy::pr_DB$get_entry_names())
  stopifnot(is.character(measure))
  stopifnot(toupper(measure) %in% available_measures)
  private$.measure <- measure
  invisible(self)
})

MapperRef$set("public", "use_cover", function(filter_values, type=c("fixed rectangular", "restrained rectangular"), ...){
  stopifnot(is.matrix(filter_values))
  stopifnot(nrow(filter_values) == nrow(private$.X))
  if (missing(type)){ type <- "fixed rectangular"}
  self$cover <- switch(type, 
    "fixed rectangular"=FixedRectangularCover$new(filter_values, ...)$construct_cover(), 
    "restrained rectangular"=RestrainedRectangularCover$new(filter_values, ...)$construct_cover(),
    stop(sprintf("Unknown cover type: %s", type))
  )
  invisible(self)
})


## Computes the distance matrix for each level set
# MapperRef$set("public", "computeLevelSetDist", function(...){
#   for (i in 1:length(self$cover$level_sets)){
#     pt_idx <- self$cover$level_sets[[i]]$points_in_level_set
#     if (length(pt_idx) <= 1){ self$cover$level_sets[[i]]$dist <- dist(0L) }
#     else {
#       if ("dist" %in% class(X)){
#         self$cover$level_sets[[i]]$dist <- dist_subset(dist = self$X, idx = pt_idx)
#       } else {
#         if (is.character(measure)){
#           self$cover$level_sets[[i]]$dist <- parallelDist::parallelDist(self$X[pt_idx,], method = measure)
#         } else if (is.function(measure)){
#           self$cover$level_sets[[i]]$dist <- measure(self$X[pt_idx,], ...)
#         } else {
#           stop("Unknown 'measure' argument. Must be either character string or function.")
#         }
#       }
#     }
#   }
#   return(self)
# })

## Computes the K-skeleton in a similar fashion as described by Singh et. al, section 3.2. Note that this
## procedure implicitly assumes the full abstract simplicial complex produced by Mapper can be treated as
## a flag complex.
MapperRef$set("public", "compute_k_skeleton", function(k, ...){
  stopifnot(k >= 0)
  if (k >= 0L){ self$compute_vertices(...) }
  if (k >= 1L){ self$compute_edges(...) }
  if (k >= 2L){
    stop("k-skeletons for k > 1 are an experimental feature!")
    # k_simplex <- vector(mode = "list", length = k - 2L)
    # for (k_i in 3L:k){
    #   pt_map <- lapply(1:length(m$nodes), function(lsfi) cbind(pt_id = m$nodes[[lsfi]], node_id = rep(lsfi)))
    #   pt_map <- data.table::data.table(do.call(rbind, pt_map))
    #   intersection_pts <- pt_map[, .(n_intersects = nrow(.SD)), by = pt_id][n_intersects == k_i]$pt_id
    #   simplex_pts <- pt_map[pt_id %in% intersection_pts, .SD, by = pt_id]
    #   index <- (seq(nrow(simplex_pts))-1) %% k_i
    #   k_simplex_mat <- unique(do.call(cbind, lapply(0:(k_i - 1L), function(idx){ simplex_pts$node_id[index == idx] })))
    #   k_simplex[[k_i - 2L]] <- lapply(1:nrow(k_simplex_mat), function(i) k_simplex_mat[i,])
    # }
  }
})

## Executes the clustering algorithm for the level sets indexed by the 'which_levels' parameter. If not given, 
## runs the clustering algorithm and computes the subsequent vertices for all the available level sets. Additional 
## parameters passed via the '...' are passed to the clustering algorithm. 
MapperRef$set("public", "compute_vertices", function(which_levels=NULL, ...){
  stopifnot(!is.na(self$cover$level_sets))
  stopifnot(is.function(private$.clustering_algorithm))
  if (!missing(which_levels) && !is.null(which_levels)){
    if (!all(which_levels %in% self$cover$index_set)){
      stop("If specified, 'which_levels' must be a vector of indexes in the covers index set.")
    }
  } else { which_levels <- self$cover$index_set }
  
  ## If not populated yet, create a default mapping from the covers index set to the vertex ids
  if (length(private$.cl_map) == 0){
    private$.cl_map <- lapply(self$cover$index_set, function(x) integer(0L))
    names(private$.cl_map) <- self$cover$index_set
  }
  
  ## Update the vertices for the given level sets. This requires updating both the simplex tree's 
  ## internal representation as well as the outward-facing vertex list. The new vertices are returned 
  ## to replace the current list. 
  which_level_idx <- match(which_levels, self$cover$index_set)-1L # 0-based
  stree_ptr <- private$.simplicial_complex$as_XPtr()
  private$.vertices <- Mapper:::build_0_skeleton(
    which_levels = which_level_idx, 
    X = private$.X, 
    f = self$clustering_algorithm, 
    level_sets = self$cover$level_sets, 
    vertices = private$.vertices, 
    ls_vertex_map = private$.cl_map, 
    stree = stree_ptr
  )
  
  ## Return self
  invisible(self)
})

## Computes the edges composing the topological graph (1-skeleton). 
## Assumes the nodes have been computed. 
MapperRef$set("public", "compute_edges", function(which_level_pairs = NULL){
  stopifnot(!is.na(self$cover$level_sets))
  if (!missing(which_level_pairs) && !is.null(which_level_pairs)){
    stopifnot(is.matrix(which_level_pairs))
    stopifnot(dim(which_level_pairs)[[2]] == 2)
    if (!all(apply(which_level_pairs, 2, function(x) x >= 1 && x <= length(self$cover$index_set)))){
      stop("If specified, 'which_level_pairs' must be an (n x 2) matrix of integers representing which level sets to compare.")
    }
  } else { which_level_pairs <- self$cover$level_sets_to_compare() }
  
  ## Retrieve the valid level set index pairs to compare. In the worst case, with no cover-specific optimization, 
  ## this may just be all pairwise combinations of LSFI's for the full simplicial complex.
  stree_ptr <- private$.simplicial_complex$as_XPtr()
  build_1_skeleton(ls_pairs = which_level_pairs-1L, vertices = private$.vertices, ls_vertex_map = private$.cl_map, stree = stree_ptr)

  ## Return self
  invisible(self)
})


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
  if ("dist" %in% class(private$.X)){ n <- attr(private$.X, "Size") } else { n <- nrow(private$.X) }
  # if (!is.null(self$config$call))
  #   cat("\nCall:\n", deparse(self$config$call), "\n\n", sep = "")
  message <- sprintf("Mapper construction for %d objects", n)
  message <- append(message, format(self$cover))
  return(message)
})


MapperRef$set("public", "plot_interactive", function(...){
  require("grapher")
  json_config <- list()
  rbw_pal <- rev(rainbow(100, start = 0, end = 4/6))
  
  ## Color nodes and edges by a default rainbow palette
  agg_node_val <- sapply(sapply(private$.vertices, function(v_idx){ 
    apply(as.matrix(self$cover$filter_values[v_idx,]), 1, mean)
  }), mean)
  binned_idx <- cut(agg_node_val, breaks = 100, labels = F)
  vertex_colors <- rbw_pal[binned_idx]
  
  ## Vertex sizes 
  vertex_sizes <- sapply(private$.vertices, length) 
  
  ## Make the igraph graph
  am <- private$.simplicial_complex$as_adjacency_matrix()
  G <- igraph::graph_from_adjacency_matrix(am, mode = "undirected", add.colnames = NA) 

  ## Normalize between 0-1, unless all the same
  normalize <- function(x) { 
    if (all(x == x[1])){ return(rep(1, length(x))) }
    else {  (x - min(x))/(max(x) - min(x)) }
  }
    
  ## Create the vertices. By default, color them on a rainbow palette according to their mean filter value.
  if (length(igraph::V(G)) > 0){
    vertex_radii <- (15L - 10L)*normalize(log(vertex_sizes)) + 10L
    vertex_xy <- apply(igraph::layout.auto(G), 2, normalize)
    json_config$nodes <- data.frame(x=vertex_xy[, 1], y=vertex_xy[, 2], r=vertex_radii,
                                    color=grapher::hex2rgba(vertex_colors),
                                    index = 0:(length(vertex_sizes)-1))
  }
    
  ## Create the edges w/ a similar coloring scheme.
  if (length(igraph::E(G)) > 0){
    el <- igraph::as_edgelist(G, names = FALSE)
    edge_binned_idx <- apply(el, 1, function(vertex_ids) { (binned_idx[vertex_ids[1]] + binned_idx[vertex_ids[2]])/2 })
    edge_links <- cbind(data.frame(matrix(apply(el, 2, as.integer) - 1L, ncol = 2)), substr(rbw_pal[edge_binned_idx], start = 0, stop = 7))
    json_config$links <- structure(edge_links, names = c("from", "to", "color"))
  }
    
  ## Return the grapher instance
  grapher::grapher(json_config)
})

# only_combinations <- function(mat){
#   mn <- pmin(mat[, 2], mat[, 1])
#   mx <- pmax(mat[, 2], mat[, 1])
#   int <- as.numeric(interaction(mn, mx))
#   mat[match(unique(int), int),]
# }

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
  # result <- self$
  # node_lsfi <- sapply(result$nodes, function(node) attr(node, "level_set"))
  # result$level_sets <- lapply(1:length(cover$level_sets), function(lsfi){ which(node_lsfi == lsfi) })
  # names(result$level_sets) <- apply(cover$index_set, 1, function(lsmi) paste0("(", paste(lsmi, collapse = ","), ")"))
  # cover_type <- paste0(toupper(substr(cover$type, start = 1, stop = 1)), tolower(substr(cover$type, start = 2, stop = nchar(cover$type))))
  # z_d <- ncol(cover$filter_values)
  # attr(result, ".summary") <- c(sprintf("Mapper object with filter function f: %s -> %s",
  #                                       ifelse(is(X, "dist"), "dist(X)",
  #                                              ifelse(ncol(X) > 1, sprintf("X^%d", ncol(X)), "X")),
  #                                       ifelse(z_d > 1, sprintf("Z^%d", z_d), "Z")),
  #                               sprintf("Configured with a %s cover comprising %d open sets", cover$type, length(cover$level_sets)),
  #                               sprintf("The 1-skeleton contains %d nodes and %d edges", length(G$nodes), sum(G$adjacency == 1L)/2L))
  # class(result) <- "Mapper"
  # return(result)
# })

## Load the exported Simplex Tree class into the package namespace
Rcpp::loadModule("simplex_tree_module", TRUE)

