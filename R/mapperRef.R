#' Mapper Reference Class (R6) implementation
#' @docType class
#' @description Composes a set of classes to perform Mapper efficiently.
#' @format An \code{\link{R6Class}} generator object
#' @keywords keyword
#' @return Instance object of the \code{\link{MapperRef}} class with methods for building the mapper.
#' 
#' @field X The data matrix.
#' @field cover The cover. 
#' @field clustering_algorithm The clustering algorithm to use. 
#' @field measure String value of the measure to use to compute distances in ambient space. Read-only. See \code{use_distance_measure} for more details.  
#' @field vertices The vertices. 
#' @field ls_vertex_map Map where the level set (by index) are keys and the vertices (by id) are the values.  
#' @field simplicial_complex A simplex tree object.
#' 
#' @section Methods:
#' \describe{
#'   \item{Documentation}{ For full documentation see.}
#'   \item{\code{new(X)}}{This method uses \code{parameter_1} to ...}
#'   \item{\code{new(X)}}{This method uses \code{parameter_1} to ...}
#' }
#' 
#' @import methods
#' @importFrom Rcpp sourceCpp
#'
#' @author Matt Piekenbrock, \email{matt.piekenbrock@@gmail.com}
#' @encoding UTF-8
#' @references Singh, Gurjeet, Facundo Memoli, and Gunnar E. Carlsson. "Topological methods for the analysis of high dimensional data sets and 3d object recognition." SPBG. 2007.
#' @useDynLib Mapper
#' @export
MapperRef <- R6::R6Class("MapperRef", 
  private = list(.X=NA, .cover=NA, .clustering_algorithm=NA, .measure=NA, .simplicial_complex = NA, .vertices=list(), .cl_map=list(), .config=NA),
  lock_class = FALSE,  ## Feel free to add your own members
  lock_objects = FALSE ## Or change existing ones 
)

MapperRef$set("public", "initialize", function(X){
  stopifnot(is.numeric(X))
  if (is.null(dim(X)) && !is(X, "dist")){ X <- as.matrix(X) }
  private$.X <- X
  private$.simplicial_complex <- simplex_tree()
  return(self)
})

## To add a public memebr function 
MapperRef$set("public", "add_function", function(name, FUN) {
  self[[name]] <- FUN
  environment(self[[name]]) <- environment(self$add_function)
})


## The level_set -> vertex mapping is available to view, but read-only
MapperRef$set("active", "ls_vertex_map", 
  function(value){ 
    if (missing(value)){ private$.cl_map } 
    else { stop("'ls_vertex_map' is read-only.") }
  }
)

## The cover stores the filter values
MapperRef$set("active", "cover", 
  function(value){ #function(fv, type = c("restrained rectangular"), ...)
    if (missing(value)){ private$.cover } 
    else {
      stopifnot(inherits(value, "CoverRef"))
      private$.cover <- value
      invisible(self)
    }
  }
)

## Mapper stores the vertices as a list
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

## Clustering algorithm must be a function
MapperRef$set("active", "clustering_algorithm", 
  function(value){
    if (missing(value)){ private$.clustering_algorithm }
    else {
      stopifnot(is.function(value))
      private$.clustering_algorithm <- value
    }
  }
)

## The data should be held fixed
MapperRef$set("active", "X", 
  function(value){
    if (missing(value)){ return(private$.X) }
    else {
      stop("`$X` is read-only. The data points 'X' are specific to a MapperRef object.")
    }
  }
)

## Active binding for the distance measure
MapperRef$set("active", "measure",
    function(value){
      if (missing(value)){ private$.measure }
      else {
        has_proxy <- requireNamespace("proxy", quietly = TRUE)
        available_measures <- if (has_proxy) proxy::pr_DB$get_entry_names() else c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
        stopifnot(is.character(value), toupper(value) %in% toupper(available_measures))
        dist_f <- ifelse(has_proxy, parallelDist::parallelDist, stats::dist)
        private$.measure <- value
      }
    }
)

#' @name use_clustering_algorithm 
#' @title Sets the clustering algorithm 
#' @description Sets a default hierarchical clustering algorithm.
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
          
          ## Use parallelDist package if available
          has_pd <- requireNamespace("parallelDist", quietly = TRUE)
          dist_f <- ifelse(has_pd, parallelDist::parallelDist, stats::dist)
          dist_x <- do.call(dist_f, list(x=X[idx,], method=self$measure))
          
          ## Use fastcluster if available 
          has_fc <- requireNamespace("fastcluster", quietly = TRUE)
          cl_f <- ifelse(has_fc, fastcluster::hclust, stats::hclust)
          hcl <- do.call(cl_f, list(d=dist_x, method=cl))
        
          ## Use the histogram binning heuristic to produce the partitioning
          cutoff_first_bin(hcl, diam = max(dist_x), num_bins)
        }
      }
      self$clustering_algorithm <- create_cl(cl = cl, num_bins.default = force(num_bins))
    } else if (is.function(cl)){
      self$clustering_algorithm <- cl
      environment(self$clustering_algorithm) <- environment(self$initialize)
    } else { stop("Invalid parameter type 'cl'") }
    invisible(self)
  }
)

## Sets the distance measure to use
## Supports any measure in proxy::pr_DB$get_entry_names()
MapperRef$set("public", "use_distance_measure", function(measure){
  self$measure <- measure
  invisible(self)
})

#' @name use_cover 
#' @title Selects one of the available covering methods to use.  
#' @description 
#' Convenience method to select, parameterize, and construct a cover and associate it with the calling objects \code{cover} field.   
#' @param filter_values (n x d) numeric matrix of values giving the results of the map.  
#' @param typename The name of the cover to use. Accepts any one of the typenames printed in the \code{available_covers}.  
#' @param ... Additional parameter values to pass to the covering generators initialize method. 
#' @details The cover is automatically constructed before being assigned to the calling Mappers \code{$cover} field.  
#' If this is undesirable, create the cover using the appropriate R6 class generators and assign to \code{$cover} field directly. 
MapperRef$set("public", "use_cover", function(filter_values, typename="fixed rectangular", ...){
  stopifnot(is.matrix(filter_values))
  stopifnot(nrow(filter_values) == nrow(private$.X))
  if (missing(typename)){ typename <- "fixed rectangular"}
  self$cover <- switch(typename, 
    "fixed rectangular"=FixedRectangularCover$new(filter_values, ...)$construct_cover(), 
    "restrained rectangular"=RestrainedRectangularCover$new(filter_values, ...)$construct_cover(),
    "adaptive"=AdaptiveCover$new(filter_values, ...)$construct_cover(),
    "ball"=BallCover$new(filter_values, ...)$construct_cover(),
    stop(sprintf("Unknown cover type: %s", typename))
  )
  invisible(self)
})

#' @name compute_k_skeleton 
#' @title Computes the K-skeleton 
#' @description For the details on how this is computed, see Singh et. al, section 3.2. 
MapperRef$set("public", "compute_k_skeleton", function(k, ...){
  stopifnot(k >= 0)
  if (k >= 0L){ self$compute_vertices(...) }
  if (k >= 1L){ self$compute_edges(...) }
  if (k >= 2L){
    stop("k-skeletons for k > 1 are an experimental feature!")
    ## TODO: offer option to use flag complex k-expansion algorithm w/ Simplex tree or 
    ## more general brute-force intersection method. 
    
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
  invisible(self)
})

#' @name compute_vertices 
#' @title Computes the vertices of mapper. 
#' @description Executes the clustering algorithm for the level sets indexed by the 'which_levels' parameter. If not given, 
#' runs the clustering algorithm and computes the subsequent vertices for all the available level sets. Additional 
#' parameters passed via the '...' are passed to the clustering algorithm. 
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

#' @name compute_edges 
#' @title Computes the edges of the mapper. 
#' @description Computes (or updates) the edges composing the topological graph (1-skeleton). 
#' @param which_level_pairs (n x 2) matrix of indices of the covers index set to check.
#' @param min_weight minimum intersection size to consider as an edge. Defaults to 1.
MapperRef$set("public", "compute_edges", function(which_level_pairs = NULL, min_weight=1L){
  stopifnot(!is.na(self$cover$level_sets))
  if (!missing(which_level_pairs) && !is.null(which_level_pairs)){
    stopifnot(is.matrix(which_level_pairs))
    stopifnot(dim(which_level_pairs)[[2]] == 2)
    stopifnot(all(unlist(which_level_pairs) %in% self$cover$index_set))
  } else { which_level_pairs <- self$cover$level_sets_to_compare() }
  
  ## Retrieve the valid level set index pairs to compare. In the worst case, with no cover-specific optimization, 
  ## this may just be all pairwise combinations of LSFI's for the full simplicial complex.
  stree_ptr <- private$.simplicial_complex$as_XPtr()
  ls_pairs <- apply(which_level_pairs, 2, function(x){ match(x, self$cover$index_set)-1L }) # 0-based
  ls_pairs <- matrix(ls_pairs, ncol = 2)
  build_1_skeleton(ls_pairs = ls_pairs, min_sz = min_weight, vertices = self$vertices, ls_vertex_map = private$.cl_map, stree = stree_ptr)

  ## Return self
  invisible(self)
})

## S3-like print override
MapperRef$set("public", "format", function(...){
  if ("dist" %in% class(private$.X)){ n <- attr(private$.X, "Size") } else { n <- nrow(private$.X) }
  message <- sprintf("Mapper construction for %d objects", n)
  if (is(self$cover, "CoverRef")){ message <- append(message, format(self$cover)) }
  return(message)
})

#' @name as_igraph 
#' @title Exports the 1-skeleton as an igraph object.
#' @description uses the igraph library. 
MapperRef$set("public", "as_igraph", function(){
  requireNamespace("igraph", quietly = TRUE)
  am <- private$.simplicial_complex$as_adjacency_matrix()
  G <- igraph::graph_from_adjacency_matrix(am, mode = "undirected", add.colnames = NA) 
  
  ## Color nodes and edges by a default rainbow palette
  rbw_pal <- rev(rainbow(100, start = 0, end = 4/6))
  agg_node_val <- sapply(sapply(private$.vertices, function(v_idx){ 
    apply(as.matrix(self$cover$filter_values[v_idx,]), 1, mean)
  }), mean)
  igraph::vertex_attr(G, name = "color") <- rbw_pal[cut(agg_node_val, breaks = 100, labels = F)]
  
  ## Normalize between 0-1, unless all the same
  normalize <- function(x) { 
    if (all(x == x[1])){ return(rep(1, length(x))) }
    else {  (x - min(x))/(max(x) - min(x)) }
  }
  vertex_sizes <- sapply(private$.vertices, length) 
  igraph::vertex_attr(G, "size") <- (15L - 10L)*normalize(log(vertex_sizes)) + 10L

  ## Fill in labels with id:size
  igraph::vertex_attr(G, "label") <- apply(cbind(names(private$.vertices), vertex_sizes), 1, function(x){
    paste0(x, collapse = ":")
  })
  
  return(G)
})

#' @name as_grapher
#' @title Exports the 1-skeleton as a grapher object.
#' @param construct_widget whether to construct the htmlwidget or just the grapher configuration.
#' @description Uses the grapher library. 
MapperRef$set("public", "as_grapher", function(construct_widget=TRUE, ...){
  requireNamespace("grapher", quietly = TRUE)
  
  ## Make the igraph 
  am <- private$.simplicial_complex$as_adjacency_matrix()
  G <- igraph::graph_from_adjacency_matrix(am, mode = "undirected", add.colnames = NA) 
  json_config <- grapher::getDefaultJsonConfig(network=G)

  ## Color nodes and edges by a default rainbow palette
  rbw_pal <- rev(rainbow(100, start = 0, end = 4/6))
  agg_node_val <- sapply(sapply(private$.vertices, function(v_idx){ 
    apply(as.matrix(self$cover$filter_values[v_idx,]), 1, mean)
  }), mean)
  binned_idx <- cut(agg_node_val, breaks = 100, labels = F)
  vertex_colors <- rbw_pal[binned_idx]
  vertex_sizes <- sapply(private$.vertices, length) 
  
  ## Normalize between 0-1, unless all the same
  normalize <- function(x) { 
    if (all(x == x[1])){ return(rep(1, length(x))) }
    else {  (x - min(x))/(max(x) - min(x)) }
  }
    
  ## Create the vertices. By default, color them on a rainbow palette according to their mean filter value.
  if (length(igraph::V(G)) > 0){
    vertex_radii <- (15L - 10L)*normalize(log(vertex_sizes)) + 10L
    vertex_xy <- apply(igraph::layout.auto(G), 2, normalize)
    json_config$nodes$x <- vertex_xy[, 1]
    json_config$nodes$y <- vertex_xy[, 2]
    json_config$nodes$r <- vertex_radii
    json_config$nodes$color <- grapher::hex2rgba(vertex_colors)
    # index = 0:(length(vertex_sizes)-1))
  } else {
    json_config$nodes <- integer(length = 0L)
  }
    
  ## Create the edges w/ a similar coloring scheme.
  if (length(igraph::E(G)) > 0){
    el <- igraph::as_edgelist(G, names = FALSE)
    edge_binned_idx <- apply(el, 1, function(vertex_ids) { (binned_idx[vertex_ids[1]] + binned_idx[vertex_ids[2]])/2 })
    edge_links <- matrix(apply(el, 2, as.integer) - 1L, ncol = 2)
    json_config$links$from <- edge_links[,1]
    json_config$links$to <- edge_links[,2]
    json_config$links$color <- substr(rbw_pal[edge_binned_idx], start = 0, stop = 7) 
  } else {
    json_config$links <- integer(length = 0L)
  }
    
  ## Return the grapher instance or just the json config if requested 
  if (construct_widget){
    grapher::grapher(json_config) 
  } else {
    return(json_config)
  }
})


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

#' @name exportMapper
#' @title Exports the essential mapper graph information 
#' @param graph_type export preference on the structure the graph output.
#' @return list with the following members: 
#' \itemize{
#' \item vertices list of the indices of the original data the current vertex intersects with.
#' \item graph some adjacency representation of the mapper graph.
#' \item level_sets map connecting the which vertices belong to the preimage of the sets in the cover. 
#' }
MapperRef$set("public", "exportMapper", function(graph_type=c("adjacency_matrix", "adjacency_list", "edgelist")){
  result <- list(level_sets = private$.cl_map, vertices = private$.vertices)
  if (missing(graph_type)){ graph_type <- "adjacency_matrix" }
  result$graph <- switch(graph_type, 
                         "adjacency_matrix"=self$simplicial_complex$as_adjacency_matrix(),
                         "adjacency_list"=self$simplicial_complex$as_adjacency_list(), 
                         "edgelist"=self$simplicial_complex$as_edge_list(), 
                         stop("'graph_type' must be one of: 'adjacency_matrix', 'adjacency_list', 'edgelist'"))
  n_simplices <- self$simplicial_complex$n_simplexes
  z_d <- ncol(self$cover$filter_values)
  attr(result, ".summary") <- c(sprintf("Mapper with filter f: %s -> %s",
                                        ifelse(is(self$X, "dist"), "dist(X)",
                                               ifelse(ncol(self$X) > 1, sprintf("X^%d", ncol(self$X)), "X")),
                                        ifelse(z_d > 1, sprintf("Z^%d", z_d), "Z")),
                                sprintf("Configured with a %s cover comprising %d open sets", self$cover$typename, length(self$cover$level_sets)),
                                sprintf("The graph contains %d vertices and %d edges", n_simplices[[1]], n_simplices[[2]]))
  class(result) <- "Mapper"
  return(result)
})

