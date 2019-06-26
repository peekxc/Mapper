#' Mapper Reference Class (R6)
#' @docType class
#' @description R6 class built for parameterizing and computing mappers efficiently.
#' @format An \code{\link{R6Class}} generator object
#' @return Instance object of the \code{\link{MapperRef}} class with methods for building the mapper.
#' 
#' @field X The data, as a function. Evaluation returns the data as a matrix.  
#' @field filter The filter, as a function. Evaluation returns the filtered data as a matrix.  
#' @field cover The cover (a \code{\link{CoverRef}} derived object). 
#' @field clustering_algorithm The clustering algorithm to use in the pullback. 
#' @field measure Distance measure to use to compute distances in ambient space. See \code{use_distance_measure} for more details.  
#' @field pullback Mapping between the sets in the cover (by index) and the vertices (by id).  
#' @field vertices The mapper vertices. 
#' @field simplicial_complex A \code{\link[simplextree:simplextree]{simplex tree}} object.
#' 
#' @details To creates a new \code{MapperRef} instance object, instantiate the class with the \code{\link[R6:R6Class]{R6 new}} operator. 
#' Instantiation of a \code{MapperRef} objects requires a data matrix, or a function returning one.
#' @section Methods:
#' \itemize{
#'   \item{\code{\link{use_filter}}: Specifies the filter.}
#'   \item{\code{\link{use_distance_measure}}: Specifies the distance measure.}
#'   \item{\code{\link{use_cover}}: Specifies the cover. Must be a \code{\link{CoverRef}} object.}
#'   \item{\code{\link{use_clustering_algorithm}}: Specifies the algorithm to decompose the pullback.}
#'   \item{\code{\link{construct_pullback}}: Decomposes the preimages into connected components.}
#'   \item{\code{\link{construct_nerve}}: Constructs simplices at a given dimension.}
#'   \item{\code{\link{construct_k_skeleton}}: Constructs the simplicial complex up to a given dimension.}
#'   \item{\code{\link{as_igraph}}: Converts the 1-skeleton to an \code{\link[igraph:igraph]{igraph object}}.}
#'   \item{\code{\link{as_pixiplex}}: Converts the simplicial complex to a \code{\link[pixiplex:pixiplex]{pixiplex object}}.}
#'   \item{\code{\link{exportMapper}}: Exports the core information of the mapper construction.}
#' }
#' @section More information:
#' Full documentation available \href{https://peekxc.github.io/Mapper}{online}.
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
  private = list(
    .X = NA, 
    .cover=NA, 
    .clustering_algorithm=NA, 
    .measure="euclidean", 
    .vertices=list(),
    .simplicial_complex = NA, 
    .pullback=list(), 
    .filter = NULL, 
    .config = NA
  ),
  lock_class = FALSE,  ## Feel free to add your own members
  lock_objects = FALSE ## Or change existing ones 
)

## To add a public member function 
## add function ----
MapperRef$set("public", "add_function", function(name, FUN) {
  self[[name]] <- FUN
  environment(self[[name]]) <- environment(self$add_function)
})

## The set index -> (vertex) decomposition mapping.
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
#' The underlying complex does not need to be modified by the user, i.e. is completely maintained by \code{\link{MapperRef}} methods 
#' (e.g. \code{\link{construct_nerve}}, \code{\link{construct_k_skeleton}}, etc.).
MapperRef$set("active", "simplicial_complex", 
  function(value){
    if (missing(value)){ private$.simplicial_complex }
    else {
      stopifnot(is(value, "Rcpp_SimplexTree"))
      private$.simplicial_complex <- value 
      # stop("`$simplicial_complex` is read-only. To change the complex, use the objects (simplex tree) methods directly.", call. = FALSE)
    }
  }
)

## clustering_algorithm ----
## Clustering algorithm must be a function that takes as arguments 'pid' and 'idx' 
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
      if (is.matrix(value)){
        X_acc <- function(data_matrix){
          function(idx=NULL){ 
            if (missing(idx)){ return(data_matrix) }
            return(data_matrix[idx,,drop=FALSE]) 
          }
        }
        private$.X <- X_acc(value)
      } else if (is.function(value)){
        private$.X <- value
      } else {
        stop("X must be either a matrix or a function that returns a matrix of coordinates.")
      }
      # stop("`$X` is read-only. The data points 'X' are specific to a MapperRef object.")
    }
  }
)

## filter ----
## The filter function associated with the Mapper instance
MapperRef$set("active", "filter", 
  function(value){
    if (missing(value)){ return(private$.filter) }
    else {
      # browser()
      if (is.matrix(value)){
        filter_acc <- function(filter_matrix){
          function(idx=NULL){ 
            if (missing(idx) || is.null(idx)){ return(filter_matrix) }
            return(filter_matrix[idx,,drop=FALSE]) 
          }
        }
        private$.filter <- filter_acc(value)
      } else if (is.function(value)){
        private$.filter <- value
      } else {
        stop("Filter must be a matrix of coordinate values or a function that returns a matrix of coordinate values.")
      }
      return(self)
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
        private$.measure <- value
      }
    }
)

## initialize ----
## Initialization method relies on active binding 
MapperRef$set("public", "initialize", function(X){
  self$X <- X
  private$.simplicial_complex <- simplextree::simplex_tree()
  self$use_distance_measure(measure = "euclidean")
  self$use_clustering_algorithm(cl = "single", cutoff_method = "continuous")
  return(self)
})

## use_filter ----
#' @name use_filter
#' @title Sets the filter
#' @description Sets the map, or \emph{filter}, to associate with the instance
#' @param filter either the filter name to use, a matrix, or a function. See details. 
#' @param ... additional parameters to pass to the filter function.
#' @details \code{filter} must be either a matrix of coordinate values, a function that returns a matrix of 
#' coordinate values, or one of the following predefined filters: 
#' \itemize{
#'   \item{\strong{PC}}{ Principle components (with \code{\link[stats:prcomp]{prcomp}})}
#'   \item{\strong{IC}}{ Independent components (with \code{\link[fastICA:fastICA]{fastICA}})}
#'   \item{\strong{ECC}}{ Eccentricity (internal, change the norm by passing one of \eqn{p} = [1, 2, Inf])}
#'   \item{\strong{KDE}}{ Kernel Density Estimate (with \code{\link[ks:ks]{kde}})}
#'   \item{\strong{DTM}}{ Distance to measure (with \code{\link[TDA:dtm]{dtm}})}
#'   \item{\strong{MDS}}{ Classic (Metric) Multidimensional Scaling (with \code{\link[stats:cmdscale]{cmdscale}})}
#'   \item{\strong{ISOMAP}}{ Isometric feature mapping (with \code{\link[vegan:isomap]{isomap}})}
#'   \item{\strong{LE}}{ Laplacian Eigenmaps (with \code{\link[geigen:geigen]{geigen}})}
#'   \item{\strong{UMAP}}{ Uniform Manifold Approximation and Projection (with \code{\link[umap:umap]{umap}})}
#' }
#' Nearly all the pre-configured filters essentially call functions in other packages with
#' somewhat reasonable default parameters to perform the mapping. Any parameters supplied to \code{...} 
#' are passed to the corresponding package function linked above, which override any default parameters. 
#' If the package needed to compute the filter is not installed, a prompt is given asking the user 
#' whether they would like to install it. \cr
#' \cr
#' \strong{NOTE:} The predefined filters are meant to be used for exploratory or illustrative purposes only---this function 
#' is \emph{not} meant to act as a comprensive interface to the functions each filter corresponds too.
MapperRef$set("public", "use_filter", function(filter=c("PC", "IC", "ECC", "KDE", "DTM", "MDS", "ISOMAP", "LE", "UMAP"), ...){
  if (is.function(filter)){ private$.filter <- filter } 
  else if (is.matrix(filter)){
    self$filter <- filter
  } else if (is.character(filter)){
    filter <- toupper(filter)
    filter_types <- c("PC", "IC", "ECC", "KDE", "DTM", "MDS", "ISOMAP", "LE", "UMAP")
    if (!filter %in% filter_types){
      stop(sprintf("Filter type %s not recognized. Must be one of: %s", filter, paste0(filter_types, collapse=", ")))
    }
    given_params <- list(...)
    require_but_ask <- function(pkg){
      pkg_installed <- requireNamespace("fastICA", quietly = TRUE)
      if (!pkg_installed){
        messsage(sprintf("Using this filter requires the package '%s' to be installed.", pkg))
        response <- readline(prompt = "Would you like to install it? y/n: ")
        if (toupper(substr(response, 1, 1)) == "Y"){ install.packages(pkg) }
      }
    }
    ## Basic filters
    if (filter == "PC"){
      default_params <- list(x = self$X(), scale. = TRUE, center = TRUE, rank. = 2L)
      params <- modifyList(default_params, given_params)
      res <- do.call(stats::prcomp, params)
      self$filter <- matrix(res$x, ncol = params[["rank."]])  
      # c("mapper_filter", "function")
    }
    else if (filter == "IC"){
      require_but_ask("fastICA")
      default_params <- list(X=self$X(), n.comp=2, method="C", alg.typ="parallel", fun="logcosh")
      params <- modifyList(default_params, given_params)
      res <- do.call(fastICA::fastICA, params)
      self$filter <- matrix(res$S, ncol = params[["n.comp"]]) ## S stores independent components
    } 
    else if (filter == "ECC"){ ## eccentricity
      has_pd <- requireNamespace("parallelDist", quietly = TRUE)
      dist_f <- ifelse(has_pd, parallelDist::parallelDist, stats::dist)
      params <- modifyList(list(p=1), given_params)
      p_str <- c("manhattan", "euclidean", "maximum")[match(params$p, list(1, 2, Inf))]
      if (params$p != Inf) { 
        self$filter <- matrix(colMeans(as.matrix(dist_f(self$X())^(1/params$p))), ncol = 1)
      } else { self$filter <- matrix(apply(as.matrix(dist_f(self$X())), 2, max), ncol = 1) } 
    }
    else if (filter == "KDE"){
      require_but_ask("ks")
      X <- self$X()
      default_params <- list(x = X, eval.points = X, verbose = FALSE)
      if (is.null(given_params[["H"]])){
        H <- if (ncol(X) <= 4L){ ks::Hpi(X) } else { diag(apply(X, 2, stats::bw.nrd0)) }
        default_params[["H"]] <- H
      }
      params <- modifyList(default_params, given_params)
      self$filter <- matrix(do.call(ks::kde, params)$estimate, ncol = 1L)
    } else if (filter == "DTM"){
      require_but_ask("TDA")
      X <- self$X()
      params <- modifyList(list(X=X, Grid=X, m0=0.20, r=2), given_params)
      self$filter <- matrix(do.call(TDA::dtm, params), ncol = 1L)
    } else if (filter == "MDS"){
      if (is.null(given_params[["d"]])){
        has_pd <- requireNamespace("parallelDist", quietly = TRUE)
        dist_f <- ifelse(has_pd, parallelDist::parallelDist, stats::dist)
        given_params[["d"]] <- dist_f(self$X())
      }
      params <- modifyList(list(k=2), given_params)
      self$filter <- matrix(do.call(stats::cmdscale, params), ncol = params[["k"]])
    } else if (filter == "ISOMAP"){
      require_but_ask("vegan")
      if (is.null(given_params[["dist"]])){
        has_pd <- requireNamespace("parallelDist", quietly = TRUE)
        dist_f <- ifelse(has_pd, parallelDist::parallelDist, stats::dist)
        given_params[["dist"]] <- dist_f(self$X())
      }
      ## Get the smallest epsilon to make the graph connected
      if (is.null(given_params[["epsilon"]]) && is.null(given_params[["k"]])){
        eps <- max(vegan::spantree(given_params[["dist"]])$dist)
        eps_add <- min(given_params[["dist"]])
        given_params[["epsilon"]] <- eps + eps_add
      }
      params <- modifyList(list(ndim=2), given_params)
      self$filter <- matrix(do.call(vegan::isomap, params)$points, ncol = params[["ndim"]])
    } else if (filter == "LE"){
      require_but_ask("geigen")
      if (is.null(given_params[["dist"]])){
        has_pd <- requireNamespace("parallelDist", quietly = TRUE)
        dist_f <- ifelse(has_pd, parallelDist::parallelDist, stats::dist)
        given_params[["dist"]] <- dist_f(self$X())
      }
      params <- modifyList(list(k=2L, sigma=mean(apply(self$X(), 2, stats::bw.nrd0))), given_params)
      W <- as.matrix(exp(-(params[["dist"]]/(params[["sigma"]]))))
      D <- diag(colSums(W))
      L <- D - W
      res <- geigen::geigen(A = L, B = D, symmetric = TRUE, only.values = FALSE)
      self$filter <- res$vector[,seq(2, 2L+(params[["k"]]-1)),drop=FALSE]
    } else if (filter == "UMAP"){
      require_but_ask("umap")
      self$filter <- do.call(umap::umap, list(d=self$X()))$layout
    } else {
      stop(sprintf("Unknown filter: %s", filter))
    }
  }
  # print.mapper_filter <- function(x){
  #   
  # }
  invisible(self)
})

## use_distance_measure ----
#' @name  use_distance_measure
#' @title Sets the distance measure to use
#' @param measure The distance measure to use (string).
#' @details Specifies the distance measure to use to compute distances when decomposing the preimages. 
#' \code{measure} must be either one of ["euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"], or, 
#' if the \pkg{proxy} package is installed, any measure in \code{proxy::pr_DB$get_entry_names()}. 
MapperRef$set("public", "use_distance_measure", function(measure){
  self$measure <- measure
  invisible(self)
})

## use_cover ----
#' @name use_cover 
#' @title Selects one of the available covering methods to use.  
#' @description 
#' Convenience method to select, parameterize, and construct a cover and associate it with the calling objects \code{cover} field.   
#' @param cover Either a pre-configured cover, or name of the cover to use. 
#' @param ... Additional parameter values to pass to the covering generators initialize method. 
#' @details The \code{cover} parameter can be any \code{\link{CoverRef}}-derived instance (e.g. a \code{\link{FixedIntervalCover}}), 
#' or one of the cover typenames listed in the \code{\link{covers_available}()}. 
#' If a typename is given, the cover is automatically constructed before being assigned to the \code{$cover} field.  
MapperRef$set("public", "use_cover", function(cover="fixed interval", ...){
  stopifnot(!is.null(self$filter))
  if (missing(cover)){ cover <- "fixed interval"}
  self$cover <- switch(cover, 
    "fixed interval"=FixedIntervalCover$new(...)$construct_cover(self$filter), 
    "restrained interval"=RestrainedIntervalCover$new(...)$construct_cover(self$filter),
    # "adaptive"=AdaptiveCover$new(...)$construct_cover(self$filter),
    "ball"=BallCover$new(...)$construct_cover(self$filter),
    stop(sprintf("Unknown cover type: %s, please specify a cover typename listed in `covers_available()`", cover))
  )
  invisible(self)
})
#@param filter_values (n x d) numeric matrix of values giving the results of the map. 

## use_clustering_algorithm ----
#' @name use_clustering_algorithm 
#' @title Sets the clustering algorithm 
#' @description Sets the clustering algorithm used to construct the connected components in the pullback cover.
#' @param cl Either one of the link criteria used in the \code{\link[stats]{hclust}} function (string) or a function. See details.
#' @param cutoff_method Type of heuristic to determine cut value. See details. Ignored is \code{cl} is a function.
#' @param ... Additional parameters passed as defaults to the cutoff method. See details. 
#' @details If \code{cl} is a linkage criterion, a standard hierarchical clustering algorithm is used with suitable default parameters. 
#' If it is a function, it must have a signature compatible with \code{function(X, idx, ...)} where \code{X} is the data matrix and \code{idx}
#' is a vector of integer indices giving which subset of \code{X} to cluster over. The enclosing environment of \code{cl}
#' is set to  \cr 
#' \cr
#' The cutoff method may be one of either c("histogram", "continuous"). See \code{\link{cutoff_first_bin}} and \code{\link{cutoff_first_threshold}}
#' for what they each correspond to, respectively. Additional named parameters passed via \code{...} act as defaults to the cutting method chosen. 
#' If not chosen, reasonable defaults are used. \cr 
#' \cr
#' Additional named parameters passed via the dots in \code{$clustering_algorithm} 
#' will taken precedence and override all previously set default settings. 
MapperRef$set("public", "use_clustering_algorithm", 
  function(cl = c("single", "ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid"), 
           cutoff_method = c("continuous", "histogram"),
           ...){
    ## Use a linkage criterion + cutting rule
    default_cutting_params <- list(...)
    if (missing(cl)) { cl <- "single" }
    if (missing(cutoff_method)){ cutoff_method <- "continuous" } 
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
      # cl_f <- create_cl(cl, cutoff_method, default_cutting_params)
      # parent.env(environment(cl_f)) <- environment(self$initialize)
      self$clustering_algorithm <- create_cl(cl, cutoff_method, default_cutting_params)
      # environment(self$clustering_algorithm) <- environment(self$initialize)
    } else if (is.function(cl)){
      ## Run given functions under the instance environment
      # parent.env(environment(cl)) <- environment(self$initialize)
      # if (environmentName(environment(cl)) == "R_GlobalEnv"){
      ## Unencloses the passed function parameters; based on pryr's 'unenclose' method 
      env <- environment(cl)
      sub <- function(x, env){
        call <- substitute(substitute(x, env), list(x = x))
        eval(call)
      }
      cl_sub <- self$clustering_algorithm <- eval(call("function", as.pairlist(formals(cl)), sub(body(cl), env)), parent.env(env))
      self$add_function("clustering_algorithm", cl_sub)
    } else { stop("Invalid parameter type 'cl'") }
    invisible(self)
  }
)

## construct_pullback ----
#' @name construct_pullback 
#' @title Constructs the (decomposed) pullback cover. 
#' @description Executes the clustering algorithm on subsets of \code{X}, 
#' @param pullback_ids indices of the \code{\link[Mapper:CoverRef]{covers}} \code{index_set}, or \code{NULL}.
#' @param ... additional parameters to pass to the \code{clustering_algorithm}.
#' @details This methods uses the function given by \code{clustering_algorithm} field to 
#' decompose the preimages returned by the \code{cover} member into connected components, which are 
#' stored as \code{vertices}. Indices may be passed to limit which sets are decomposed, otherwise 
#' the sets in the cover all considered.
#' \cr
#' Note that this method removes the \code{vertices} associated with \code{pullback_ids}, but does not 
#' modify the simplicial complex. 
#' @seealso \code{\link{construct_k_skeleton}} \code{\link{construct_nerve}} \code{\link{use_clustering_algorithm}}
MapperRef$set("public", "construct_pullback", function(pullback_ids=NULL, ...){
  stopifnot(!is.na(self$cover$level_sets))
  stopifnot(is.function(private$.clustering_algorithm))
  pids_supplied <- (!missing(pullback_ids) && !is.null(pullback_ids))
  if (pids_supplied){ stopifnot(all(pullback_ids %in% self$cover$index_set)) }
  pullback_ids <- if (pids_supplied){ pullback_ids } else { self$cover$index_set }
  
  ## If specific pullback ids not given, then resist the pullback mapping
  if (!pids_supplied || length(private$.pullback) == 0){
    n_sets <- length(self$cover$index_set)
    self$simplicial_complex$clear()
    self$vertices <- list()
    private$.pullback <- structure(replicate(n_sets, integer(0)), names = self$cover$index_set)
  }
  
  ## Update the vertices for the given level sets. This requires updating both the simplex tree's 
  ## internal representation as well as the outward-facing vertex list. The new vertices are returned 
  ## to replace the current list. The pullback map is also updated. 
  calc_preimage <- function(){
    function(index){ self$cover$construct_cover(self$filter, index) }
  }
  partial_cluster <- function(...){
    extra <- list(...)
    function(pid, idx){ do.call(self$clustering_algorithm, append(list(pid, idx), extra)) }
  }
  # id_gen <- function(n){
  #   ## make cpp function smallest_not_in
  #   private$.vertices
  #   as.integer(names(self$vertices))
  #   private$.simplicial_complex
  #   ids <- private$.simplicial_complex$generate_ids(n)
  #   private$.simplicial_complex$insert(as.list(ids))
  #   return(ids)
  # }
  
  ## Do the decomposition
  private$.vertices <- Mapper:::decompose_preimages(
    pullback_ids = as.character(pullback_ids),
    cluster_f = partial_cluster(...), 
    level_set_f = calc_preimage(), 
    vertices = private$.vertices,
    # id_generator = id_gen,
    pullback = private$.pullback
  )
  
  ## Return self
  invisible(self)
})

## construct_nerve ----
#' @name construct_nerve 
#' @title Compute the nerve of the cover. 
#' @description Computes (or updates) the k-simplices composing the Mapper, where k >= 0. 
#' @param k The order of the simplices to construct. See details.  
#' @param indices (n x k) matrix of indices of the covers index set to update.
#' @param min_weight minimum intersection size to consider as a simplex. Defaults to 1.
#' @details Compared to \code{construct_k_skeleton}, this method \emph{only} intersections 
#' between (k-1) simplices in the complex.
MapperRef$set("public", "construct_nerve", function(k, indices = NULL, min_weight=1L, modify=TRUE){
  stopifnot(length(private$.vertices) > 0)
  stopifnot(k == trunc(k), k >= 0)
  idx_specified <- (!missing(indices) && !is.null(indices))
  if (idx_specified){
    stopifnot(is.matrix(indices))
    stopifnot(all(unlist(indices) %in% self$cover$index_set))
    stopifnot(k >= 1)
  }
  indices <- if (idx_specified){ indices } else { self$cover$neighborhood(self$filter, k) }
  ## Retrieve the valid level set index pairs to compare. In the worst case, with no cover-specific optimization, 
  ## this may just be all pairwise combinations of LSFI's for the full simplicial complex.
  stree_ptr <- private$.simplicial_complex$as_XPtr()
  params <- list(pullback_ids = indices, vertices = self$vertices, pullback = private$.pullback, stree = stree_ptr, modify = modify)
  if (k == 0){ 
    if (!params$modify){ return(as.integer(names(self$vertices))) } 
    build_0_skeleton(as.integer(names(self$vertices)), stree_ptr)
    return(invisible(self))
  }
  else if (k == 1){ 
    params <- modifyList(params, list(min_sz = min_weight)) 
    if (!params$modify){ return(do.call(rbind, do.call(build_1_skeleton, params)) ) } 
    do.call(build_1_skeleton, params)
    return(invisible(self))
  } 
  else { 
    params <- modifyList(params, list(k = k)) 
    if (!params$modify){ return(do.call(rbind, do.call(build_k_skeleton, params)) ) } 
    do.call(build_k_skeleton, params)
    return(invisible(self))
  }  
})

## construct_k_skeleton ----
#' @name construct_k_skeleton 
#' @title Constructs the k-skeleton
#' @param k the maximal dimension to consider.
#' @description Computes the k-skeleton of the mapper by computing the nerve of the 
#' the pullback cover induced by the \code{cover} member. \cr
#' \cr
#' This function computes the k-skeleton inductively, e.g. by first computing the vertices, 
#' then the edges, etc. A check is performed to ensure the pullback has been decomposed, and 
#' if not, then \code{\link{construct_pullback}} is called.
#' @details For the details on how this is computed, see Singh et. al, section 3.2. 
MapperRef$set("public", "construct_k_skeleton", function(k=1L){
  stopifnot(k >= 0)
  if (length(private$.vertices) == 0) { self$construct_pullback() }
  for (k_i in seq(0L, k)){ self$construct_nerve(k = k_i) }
  invisible(self)
})

## format ----
## S3-like print override
MapperRef$set("public", "format", function(...){
  #if ("dist" %in% class(private$.X)){ n <- attr(private$.X, "Size") } else { n <- nrow(private$.X) }
  max_k <- length(private$.simplicial_complex$n_simplices)
  if (max_k == 0){ message <- "(empty) Mapper construction" }
  else {
    simplex_info <- sprintf("(%s) (%s)-simplices", paste0(private$.simplicial_complex$n_simplices, collapse = ", "), paste0(0L:(max_k-1L), collapse = ", "))
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
    mean(rowMeans(self$filter(v_idx)))
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

# MapperRef$set("public", "as_grapher", function(construct_widget=TRUE, ...){
#   requireNamespace("grapher", quietly = TRUE)
#   
#   ## Make the igraph 
#   am <- private$.simplicial_complex$as_adjacency_matrix()
#   G <- igraph::graph_from_adjacency_matrix(am, mode = "undirected", add.colnames = NA) 
#   json_config <- grapher::getDefaultJsonConfig(network=G)
# 
#   ## Color nodes and edges by a default rainbow palette
#   rbw_pal <- rev(rainbow(100, start = 0, end = 4/6))
#   agg_node_val <- sapply(sapply(private$.vertices, function(v_idx){ 
#     apply(as.matrix(self$cover$filter_values[v_idx,]), 1, mean)
#   }), mean)
#   binned_idx <- cut(agg_node_val, breaks = 100, labels = F)
#   vertex_colors <- rbw_pal[binned_idx]
#   vertex_sizes <- sapply(private$.vertices, length) 
#   
#   ## Normalize between 0-1, unless all the same
#   normalize <- function(x) { 
#     if (all(x == x[1])){ return(rep(1, length(x))) }
#     else {  (x - min(x))/(max(x) - min(x)) }
#   }
#     
#   ## Create the vertices. By default, color them on a rainbow palette according to their mean filter value.
#   if (length(igraph::V(G)) > 0){
#     vertex_radii <- (7L - 2L)*normalize(vertex_sizes) + 2L
#     vertex_xy <- apply(igraph::layout.auto(G), 2, normalize)
#     json_config$nodes$x <- vertex_xy[, 1]
#     json_config$nodes$y <- vertex_xy[, 2]
#     json_config$nodes$r <- vertex_radii
#     json_config$nodes$color <- grapher::hex2rgba(vertex_colors)
#     # index = 0:(length(vertex_sizes)-1))
#   } else {
#     json_config$nodes <- integer(length = 0L)
#   }
#     
#   ## Create the edges w/ a similar coloring scheme.
#   if (length(igraph::E(G)) > 0){
#     el <- igraph::as_edgelist(G, names = FALSE)
#     edge_binned_idx <- apply(el, 1, function(vertex_ids) { (binned_idx[vertex_ids[1]] + binned_idx[vertex_ids[2]])/2 })
#     edge_links <- matrix(apply(el, 2, as.integer) - 1L, ncol = 2)
#     json_config$links$from <- edge_links[,1]
#     json_config$links$to <- edge_links[,2]
#     json_config$links$color <- substr(rbw_pal[edge_binned_idx], start = 0, stop = 7) 
#   } else {
#     json_config$links <- integer(length = 0L)
#   }
#     
#   ## Return the grapher instance or just the json config if requested 
#   if (construct_widget){
#     grapher::grapher(json_config) 
#   } else {
#     return(json_config)
#   }
# })

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
    rowMeans(self$filter(v_idx))
  }), mean)
  binned_idx <- cut(agg_node_val, breaks = 100, labels = F)
  vertex_colors <- rbw_pal[binned_idx]
  vertex_sizes <- sapply(private$.vertices, length) 
  
  ## Normalize between 0-1, unless all the same
  normalize <- function(x) { 
    if (all(x == x[1])){ return(rep(1, length(x))) }
    else {  (x - min(x))/(max(x) - min(x)) }
  }
  idx_match <- match(pp$nodes$id, as.integer(names(private$.vertices)))
  pp$nodes$color <- vertex_colors[idx_match]
  pp$nodes$radius <- ((15L - 5L)*normalize(vertex_sizes) + 5L)[idx_match]
  return(pp)
})

## exportMapper ----
#' @name exportMapper
#' @title Exports minimal information about the Mapper
#' @description This function exports a few core characteristics about the mapper complex as a list.
#' @param graph_type export preference on the structure the graph output.
#' @return List with the following members: 
#' \itemize{
#'   \item \emph{vertices} list of the indices of the original data the current vertex intersects with.
#'   \item \emph{graph} some adjacency representation of the mapper graph.
#'   \item \emph{pullback} map connecting the index set of the cover and the vertices.
#' }
#' @details \code{graph_type} must be one of 'adjacency_matrix', 'adjacency_list', or 'edgelist'. 
MapperRef$set("public", "exportMapper", function(graph_type=c("adjacency_matrix", "adjacency_list", "edgelist")){
  result <- list(pullback = private$.pullback, vertices = private$.vertices)
  if (missing(graph_type)){ graph_type <- "adjacency_matrix" }
  result$graph <- switch(graph_type, 
                         "adjacency_matrix"=self$simplicial_complex$as_adjacency_matrix(),
                         "adjacency_list"=self$simplicial_complex$as_adjacency_list(), 
                         "edgelist"=self$simplicial_complex$as_edge_list(), 
                         stop("'graph_type' must be one of: 'adjacency_matrix', 'adjacency_list', 'edgelist'"))
  n_simplices <- self$simplicial_complex$n_simplices
  x_dim <- ncol(self$X())
  z_dim <- ncol(self$filter())
  attr(result, ".summary") <- c(sprintf("Mapper with filter f: %s -> %s",
                                        ifelse(x_dim > 1, sprintf("X^%d", x_dim), "X"),
                                        ifelse(z_dim > 1, sprintf("Z^%d", z_dim), "Z")),
                                format(self$cover))
  class(result) <- "Mapper"
  return(result)
})

