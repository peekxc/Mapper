#' Mapper Reference Class (R6) implementation
#' @docType class
#' @description R6 utility class that enables computing mappers efficiently.
#' @format An \code{\link{R6Class}} generator object
#' @keywords keyword
#' @return Instance object of the \code{\link{MapperRef}} class with methods for building the mapper.
#' 
#' @field X The data matrix.
#' @field cover The cover. 
#' @field clustering_algorithm The clustering algorithm to use in the pullback. 
#' @field measure String value of the distance measure to use to compute distances in ambient space. Read-only. See \code{use_distance_measure} for more details.  
#' @field pullback Mapping between the sets in the cover (by index) and the vertices (by id).  
#' @field vertices The mapper vertices. 
#' @field simplicial_complex A \code{\link[simplextree:simplextree]{simplex tree}} object.
#' 
#' @section Methods:
#' \describe{
#'   \item{Documentation}{Full documentation available \href{https://peekxc.github.io/Mapper}{online}.}
#'   \item{\code{new(X)}}{This method uses \code{parameter_1} to ...}
#'   \item{\code{new(X)}}{This method uses \code{parameter_1} to ...}
#' }
#' 
#' @import methods
#' @importFrom Rcpp sourceCpp
#'
#' @author Matt Piekenbrock, \email{matt.piekenbrock@@gmail.com}
#' @encoding UTF-8
#' @references Gurjeet Singh, Facundo Mémoli, and Gunnar Carlsson. "Topological methods for the analysis of high dimensional data sets and 3d object recognition." SPBG. 2007.
#' @useDynLib Mapper
#' @export
MapperRef <- R6::R6Class("MapperRef", 
  private = list(.X=NA, .cover=NA, .clustering_algorithm=NA, .measure="euclidean", .simplicial_complex = NA, .vertices=list(), .pullback=list(), .config=NA),
  lock_class = FALSE,  ## Feel free to add your own members
  lock_objects = FALSE ## Or change existing ones 
)

## initialize ----
MapperRef$set("public", "initialize", function(X){
  # 
  if (!is.function(X)){
    stopifnot(is.numeric(X))
    if (is.null(dim(X)) && !is(X, "dist")){ X <- as.matrix(X) }
    Xf <- function(data_matrix){
      function(idx){ return(data_matrix[idx,,drop=FALSE]) }
    }
    private$.X <- Xf(X)
  } else { private$.X <- X }
  private$.simplicial_complex <- simplextree::simplex_tree()
  self$use_distance_measure(measure = "euclidean")
  self$use_clustering_algorithm(cl = "single", cutoff_method = "continuous")
  return(self)
})

## To add a public member function 
## add function ----
MapperRef$set("public", "add_function", function(name, FUN) {
  self[[name]] <- FUN
  environment(self[[name]]) <- environment(self$add_function)
})

## The set index -> (vertex) decomposition mapping. Read-only
## pullback ----
MapperRef$set("active", "pullback", 
  function(value){ 
    if (missing(value)){ private$.pullback } 
    else { private$.pullback <- value }
  }
)

## cover ----
#' @title Mapper cover
#' @name cover
#' @description Every \code{\link{MapperRef}} object requires a \code{\link{CoverRef}} object as 
#' is \code{cover} member field. In the context of Mapper, a cover is used to discretize the filter 
#' space into a partition, which is then used via a \emph{pullback} operation to construct the vertices. \cr 
#' \cr 
#' The \code{\link{MapperRef}} class makes no restrictions on the cover that is used; only that it fulfills the 
#' requirements of being a valid \code{\link{CoverRef}} instance.
#' @seealso \code{\link{CoverRef}} 
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
## vertices ----
MapperRef$set("active", "vertices", 
  function(value){
    if(missing(value)){
      private$.vertices
    } else {
      private$.vertices <- value
      # stop("`$vertices` is read-only. To update the vertex membership, use the 'compute_vertices' function.", call. = FALSE)
    }
  }  
)

## simplicial_complex ----
#' @title Simplicial Complex 
#' @name simplicial_complex 
#' @description The relational information of the Mapper construction. 
#' @details The primary output of the Mapper method is a simplicial complex. With \code{\link{MapperRef}} objects, 
#' the simplicial complex is stored as a \code{\link[simplextree]{simplextree}}. 
#' \cr 
#' The underlying complex is completely maintained by \code{\link{MapperRef}} methods 
#' (e.g. \code{\link{compute_vertices}}, etc.). The complex may also be modified directly
#' via \code{\link[simplextree]{simplextree}} methods. 
MapperRef$set("active", "simplicial_complex", 
  function(value){
    if (missing(value)){ private$.simplicial_complex }
    else {
      private$.simplicial_complex <- value 
      # stop("`$simplicial_complex` is read-only. To change the complex, use the objects (simplex tree) methods directly.", call. = FALSE)
    }
  }
)

## clustering_algorithm ----
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

## X ----
## The data should be held fixed
MapperRef$set("active", "X", 
  function(value){
    if (missing(value)){ return(private$.X) }
    else {
      stop("`$X` is read-only. The data points 'X' are specific to a MapperRef object.")
    }
  }
)

## measure ----
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

## use_clustering_algorithm ----
#' @name use_clustering_algorithm 
#' @title Sets the clustering algorithm 
#' @description Sets the clustering algorithm used to construct the connected components in the pullback cover.
#' @param cl Either one of the link criteria used in the \code{\link[stats]{hclust}} function (string) or a function. See details.
#' @param cutoff_method Type of heuristic to determine cut value. See details. Ignored is \code{cl} is a function.
#' @param run_internal Whether to run the clustering within the \code{MapperRef} instance environment.
#' @param ... Additional parameters passed as defaults to the cutoff method. See details. 
#' @details If \code{cl} is a linkage criterion, a standard hierarchical clustering algorithm is used with suitable default parameters. 
#' If it is a function, it must have a signature compatible with \code{function(X, idx, ...)} where \code{X} is the data matrix and \code{idx}
#' is a vector of integer indices giving which subset of \code{X} to cluster over. \cr 
#' \cr
#' The cutoff method may be one of either c("histogram", "continuous"). See \code{\link{cutoff_first_bin}} and \code{\link{cutoff_first_threshold}}
#' for what they each correspond to, respectively. Additional named parameters passed via \code{...} act as defaults to the cutting method chosen. 
#' If not chosen, reasonable defaults are used. \cr 
#' \cr 
#' If \code{run_internal} is set to TRUE, the \code{MapperRef} instance environment is imported, such that field members
#' may be used inside \code{cl}s definition. Defaults to FALSE. 
#' \cr
#' Additional named parameters passed via the dots in \code{$clustering_algorithm} 
#' will taken precedence and override all previously set default settings. 
MapperRef$set("public", "use_clustering_algorithm", 
  function(cl = c("single", "ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid"), 
           cutoff_method = c("continuous", "histogram"),
           run_internal = FALSE,
           ...){
    ## Use a linkage criterion + cutting rule
    default_cutting_params <- list(...)
    if (missing(cl)) { cl <- "single" }
    if (missing(cutoff_method)){ cutoff_method <- "histogram" } 
    if (class(cl) == "character"){
      hclust_opts <- c("single", "ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
      if (!cl %in% hclust_opts){ stop(sprint("Unknown linkage method passed. Please use one of (%s). See ?hclust for details.", paste0(hclust_opts, collapse = ", "))) }
      stopifnot(is.character(self$measure))
      ## Closured representation to substitute default parameters. 
      create_cl <- function(cl, cutoff_method, cut_defaults){
        cutoff_f <- switch(cutoff_method, "histogram"=cutoff_first_bin, "continuous"=cutoff_first_threshold)
        function(pid, idx, ...){
          if (is.null(idx) || length(idx) == 0){ return(integer(0L)) }
          if (length(idx) <= 2L){ return(rep(1L, length(idx))); }
          override_params <- list(...)
          
          ## Use parallelDist package if available
          has_pd <- requireNamespace("parallelDist", quietly = TRUE)
          dist_f <- ifelse(has_pd, parallelDist::parallelDist, stats::dist)
          dist_x <- do.call(dist_f, list(x=self$X(idx), method=self$measure))
          
          ## Use fastcluster if available 
          has_fc <- requireNamespace("fastcluster", quietly = TRUE)
          cl_f <- ifelse(has_fc, fastcluster::hclust, stats::hclust)
          hcl <- do.call(cl_f, list(d=dist_x, method=cl))
          
          cutoff_params <- if (cutoff_method == "histogram"){ 
            list(hcl = hcl, num_bins = 10L) 
          } else { list(hcl = hcl, threshold = 0.0) }
          cutoff_params <- modifyList(modifyList(cutoff_params, cut_defaults), override_params)
        
          ## Use heuristic cutting value to produce the partitioning
          eps <- do.call(cutoff_f, cutoff_params)
          return(cutree(hcl, h = eps))
        }
      }
      cl_f <- create_cl(cl, cutoff_method, default_cutting_params)
      parent.env(environment(cl_f)) <- environment(self$initialize)
      self$clustering_algorithm <- cl_f
      # environment(self$clustering_algorithm) <- environment(self$initialize)
    } else if (is.function(cl)){
      self$clustering_algorithm <- cl
      ## If internal objects are requested
      if (run_internal){
        parent.env(environment(self$clustering_algorithm)) <- environment(self$initialize)
        # environment(self$clustering_algorithm) <- environment(self$initialize)
      }
    } else { stop("Invalid parameter type 'cl'") }
    invisible(self)
  }
)

## use_distance_measure ----
## Sets the distance measure to use
## Supports any measure in proxy::pr_DB$get_entry_names()
MapperRef$set("public", "use_distance_measure", function(measure){
  self$measure <- measure
  invisible(self)
})

## use_cover ----
#' @name use_cover 
#' @title Selects one of the available covering methods to use.  
#' @description 
#' Convenience method to select, parameterize, and construct a cover and associate it with the calling objects \code{cover} field.   
#' @param filter_values (n x d) numeric matrix of values giving the results of the map.  
#' @param typename The name of the cover to use. Accepts any one of the typenames printed in the \code{available_covers}.  
#' @param ... Additional parameter values to pass to the covering generators initialize method. 
#' @details The cover is automatically constructed before being assigned to the calling Mappers \code{$cover} field.  
#' If this is undesirable, create the cover using the appropriate R6 class generators and assign to \code{$cover} field directly. 
MapperRef$set("public", "use_cover", function(filter_values, typename="fixed interval", ...){
  stopifnot(is.matrix(filter_values))
  stopifnot(nrow(filter_values) == nrow(private$.X))
  if (missing(typename)){ typename <- "fixed interval"}
  self$cover <- switch(typename, 
    "fixed interval"=FixedIntervalCover$new(filter_values, ...)$construct_cover(), 
    "restrained interval"=RestrainedIntervalCover$new(filter_values, ...)$construct_cover(),
    "adaptive"=AdaptiveCover$new(filter_values, ...)$construct_cover(),
    "ball"=BallCover$new(filter_values, ...)$construct_cover(),
    stop(sprintf("Unknown cover type: %s", typename))
  )
  invisible(self)
})

## compute_k_skeleton ----
#' @name compute_k_skeleton 
#' @title Compute the K-skeleton 
#' @description Computes the k-skeleton of the mapper by computing the nerve of the 
#' the pullback cover induced by the cover set by the \code{cover} member. \cr
#' \cr
#' This function computes the k-skeleton inductively, e.g. by first computing the vertices, 
#' then the edges, etc. Computing higher dimensions is supported, however may be significantly 
#' slower than computing the 1-skeleton. 
#' @details For the details on how this is computed, see Singh et. al, section 3.2. 
MapperRef$set("public", "construct_k_skeleton", function(k=1L, ...){
  stopifnot(k >= 0)
  self$simplicial_complex$clear()
  if (k >= 0L){ self$construct_pullback(...) }
  if (k >= 1L){
    for (k_i in seq(1L, k)){ self$construct_nerve(k = k_i) }
  }
  invisible(self)
})

## construct_pullback ----
#' @name construct_pullback 
#' @title Computes the vertices of the mapper at a given index. 
#' @description Executes the clustering algorithm for sets indexed by \code{pullback_ids}. 
#' @details If ids are not supplied, the pullback is computed for all the sets in the cover.
#' parameters passed via the '...' are passed to the clustering algorithm. \cr
#' \cr
#' Note that this method removes the vertices associated with \code{pullback_ids} from the simplicial 
#' complex, which subsequently also any of their cofaces. New vertices computed with \code{clustering_algorithm}
#' have their ids generated by \code{generate_ids}.
MapperRef$set("public", "construct_pullback", function(pullback_ids=NULL, ...){
  stopifnot(!is.na(self$cover$level_sets))
  stopifnot(is.function(private$.clustering_algorithm))
  pids_supplied <- (!missing(pullback_ids) && !is.null(pullback_ids))
  if (pids_supplied){ stopifnot(all(pullback_ids %in% self$cover$index_set)) }
  pullback_ids <- if (pids_supplied){ pullback_ids } else { self$cover$index_set }
  
  ## If not populated yet, create a default pullback mapping from the covers index set to the vertex ids
  if (length(private$.pullback) == 0){
    private$.pullback <- lapply(self$cover$index_set, function(x) integer(0L))
    names(private$.pullback) <- self$cover$index_set
  }
  
  ## Update the vertices for the given level sets. This requires updating both the simplex tree's 
  ## internal representation as well as the outward-facing vertex list. The new vertices are returned 
  ## to replace the current list. The pullback map is also updated. 
  private$.vertices <- Mapper:::build_0_skeleton(
    pullback_ids = as.character(pullback_ids),
    cluster_f = self$clustering_algorithm, 
    level_set_f = self$cover$construct_cover, 
    vertices = private$.vertices, 
    pullback = private$.pullback, 
    stree = private$.simplicial_complex$as_XPtr()
  )
  
  ## Return self
  invisible(self)
})

## construct_nerve ----
#' @name construct_nerve 
#' @title Compute the nerve of the cover. 
#' @description Computes (or updates) the k-simplices composing the Mapper, where k >= 1. 
#' @param k The order of the simplices to construct. See details.  
#' @param indices (n x k) matrix of indices of the covers index set to update.
#' @param min_weight minimum intersection size to consider as a simplex. Defaults to 1.
#' @details Compared to \code{construct_k_skeleton}, this method \emph{only} intersections 
#' between (k-1) simplices in the complex.
MapperRef$set("public", "construct_nerve", function(k = 1, indices = NULL, min_weight=1L, modify=TRUE){
  stopifnot(length(private$.vertices) > 0)
  stopifnot(k >= 1L)
  idx_specified <- (!missing(indices) && !is.null(indices))
  if (idx_specified){
    stopifnot(is.matrix(indices))
    stopifnot(all(unlist(indices) %in% self$cover$index_set))
    k <- ncol(indices)-1L ## ensure k is equal to the number of columns of the indices - 1
  }
  indices <- if (idx_specified){ indices } else { self$cover$neighborhood(k) }
  
  
  ## Retrieve the valid level set index pairs to compare. In the worst case, with no cover-specific optimization, 
  ## this may just be all pairwise combinations of LSFI's for the full simplicial complex.
  stree_ptr <- private$.simplicial_complex$as_XPtr()
  params <- list(pullback_ids = indices, vertices = self$vertices, pullback = private$.pullback, stree = stree_ptr, modify = modify)
  if (k == 1){ params <- modifyList(params, list(min_sz = min_weight)) } 
  else { params <- modifyList(params, list(k = k)) }
  sk_f <- ifelse(k == 1, build_1_skeleton, build_k_skeleton)
  if (params$modify){
    do.call(sk_f, params)
    return(invisible(self))
  } else {
    return( do.call(rbind, do.call(sk_f, params)) )
  }
})

## format ----
## S3-like print override
MapperRef$set("public", "format", function(...){
  #if ("dist" %in% class(private$.X)){ n <- attr(private$.X, "Size") } else { n <- nrow(private$.X) }
  max_k <- length(private$.simplicial_complex)
  if (max_k == 0){ message <- "Mapper construction (empty)" }
  else {
    simplex_info <- sprintf("(%s) (%s)-simplices\n", paste0(private$.simplicial_complex$n_simplices, collapse = ", "), paste0(0L:(max_k-1L), collapse = ", "))
    message <- paste0("Mapper construction with ", simplex_info)
  }
  if (is(self$cover, "CoverRef")){ message <- append(message, format(self$cover)) }
  return(message)
})

## as_igraph ----
#' @name as_igraph 
#' @title Exports Mapper as an igraph object.
#' @description Exports the 1-skeleton to a graph using the igraph library.
#' @param vertex_scale scaling function for the vertex sizes. 
#' @param vertex_min minimum vertex size. 
#' @param vertex_min maximum vertex size. 
#' @param col_pal color palette to color the vertices by. 
#' @details By default, the vertices of the output graph are given "color", "size", and "label" 
#' attributes. The vertex colors are colored according on a blue to red rainbow according 
#' to their mean filter value (see \code{\link{bin_color}}). The vertex sizes are scaled according 
#' to the number of points they contain, scaled by \code{vertex_scale}, and bounded between 
#' (\code{vertex_min}, \code{vertex_max}). The vertex labels are in the format "<id>:<size>".
MapperRef$set("public", "as_igraph", function(vertex_scale=c("linear", "log"), vertex_min=10L, vertex_max=15L, col_pal="rainbow"){
  requireNamespace("igraph", quietly = TRUE)
  am <- private$.simplicial_complex$as_adjacency_matrix()
  G <- igraph::graph_from_adjacency_matrix(am, mode = "undirected", add.colnames = NA) 
  
  ## Color nodes and edges by a default rainbow palette
  agg_node_val <- sapply(private$.vertices, function(v_idx){ 
    mean(rowMeans(self$cover$filter_values[v_idx,,drop=FALSE]))
  })
  igraph::vertex_attr(G, name = "color") <- bin_color(agg_node_val, col_pal)
  
  ## Normalize between 0-1, unless all the same
  normalize <- function(x) { 
    if (all(x == x[1])){ return(rep(1, length(x))) }
    else {  (x - min(x))/(max(x) - min(x)) }
  }
  if (missing(vertex_scale)){ vertex_scale <- "linear"}
  vertex_scale <- switch(vertex_scale, "linear"=identity, "log"=log)
  vertex_sizes <- sapply(private$.vertices, length) 
  igraph::vertex_attr(G, "size") <- (vertex_max - vertex_min)*normalize(vertex_scale(vertex_sizes)) + vertex_min

  ## Fill in labels with id:size
  igraph::vertex_attr(G, "label") <- apply(cbind(names(private$.vertices), vertex_sizes), 1, function(x){
    paste0(x, collapse = ":")
  })
  return(G)
})

## as_grapher ----
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
    vertex_radii <- (7L - 2L)*normalize(vertex_sizes) + 2L
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

## as_pixiplex ----
#' @name as_pixiplex
#' @title Exports the complex as a pixiplex object.
#' @description Uses the pixiplex library. 
MapperRef$set("public", "as_pixiplex", function(...){
  requireNamespace("pixiplex", quietly = TRUE)
  pp <- pixiplex::pixiplex(private$.simplicial_complex)
  
  ## Default styling
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
  pp$nodes$color <- vertex_colors
  pp$nodes$radius <- (15L - 5L)*normalize(vertex_sizes) + 5L
  return(pp)
})

## exportMapper ----
#' @name exportMapper
#' @title Exports mapper information 
#' @param graph_type export preference on the structure the graph output.
#' @return list with the following members: 
#' \itemize{
#' \item \emph{vertices} list of the indices of the original data the current vertex intersects with.
#' \item \emph{graph} some adjacency representation of the mapper graph.
#' \item \emph{level_sets} map connecting the which vertices belong to the preimage of the sets in the cover. 
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

