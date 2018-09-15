#' @title mapper
#'
#' @description Computes the topological graph of the Mapper construction. Mapper is a tool for summarizing topological information from datasets and
#' maps defined on them. It takes as input a set of 'point cloud' data and a (possibly low dimensional) coordinatization of the data,
#' and produces a topological summary of the data expressed through a cover of the codomain of the map.
#'
#' @param X Either an \eqn{n x m} data matrix (preferred) or a \code{\link[stats]{dist}} object.
#' @param filter_values An \eqn{n x d} matrix (or vector) of values produced .
#' @param cover_type Name of the covering method to use.
#' @param return_reference Boolean whether or not to return the reference class used to construct Mapper. See the return value section.
#' @param ... Additional arguments passed to other components in the Mapper framework. See details.
#' @details
#' \code{mapper} is a generic function, whose parameterization allows complete control on the resulting construction.
#' Additional named parameters needed by any of these components may be passed through the three-dot ellipsis \code{...}.
#' See the examples section for usage cases.
#'
#' \strong{Additional arguments}{
#' \describe{
#' \item{\code{number_intervals}}{ Either a positive integer or a vector of \emph{m} positive integers indicating the number of intervals to distribute along the filter space.}
#' \item{\code{overlap}}{ Either a single real-valued number or m-length vector of values between 0 and 1 specifying how much adjacent intervals should overlap.}
#' \item{\code{measure}}{ If \code{X} is a point cloud, string-value indicating the measure to compute on \code{X}. Ignored if \code{X} is a \code{dist} object. }
#' \item{\code{cl}}{ String-value of the linkage criterion to using for the hierarchical clustering step, or a clustering function which returns an integer vector indicating cluster membership. See details. }
#' \item{\code{...}}{ All named parameters are dispatched to other functional components of Mapper. See examples. }
#' }
#' }
#' \code{X} may be a \code{dist} object or point cloud data in a matrix-coercible format (recommended).
#'
#' The \code{cl} argument controls the clustering algorithm to use. If a string is given, any named linkage criterion in \code{\link[stats]{hclust}} may be used.
#' Alternatively, a function may be given which accepts as input a \code{\link[stats]{dist}} object and returns an integer vector partitioning of the data.
#' By default, if a hierarchical algorithm is used, a histogram-based heuristic is used to find the cutting threshold. The \code{num_bins} parameter may be
#' given to use an alternative binning for the histogram (default is 10). See ?\code{\link[Mapper]{cutoff_first_bin}} for details.
#'
#' The \code{cover_type} can be either "fixed rectangular" (default) or "restrained rectangular." Both covers rely on the same parameters, \code{number_intervals} and
#' \code{overlap}.
#'
#' If \code{X} is a point cloud and hierarchical clustering is used, distances must be computed on subsets of \code{X}. The \code{measure} parameter controls which
#' distance measure to use. Supports any distance measure in \code{\link[parallelDist]{parallelDist}} package. Alternatively, \code{measure} may be given as a function
#' which accepts as input the subset of \code{X} to compute distances on.
#'
#' All named parameters in \code{...} are dispatched to their appropriate methods, and thus may be used to configure custom functions as well. For example, if a
#' custom clustering function is supplied to the \code{cl} argument, named parameters in \code{...} will be passed to the clustering algorithm as needed.
#'
#'
#' @section Return Value:
#' By default, returns a \code{mapper} object, which is a named list with members \code{adjacency} (adjacency matrix for the edges),
#' \code{nodes} (list with point membership indices), and list \code{level_sets} (list with node indices).
#' Each node in \code{nodes} also records the flat index of the level set that node intersects, recorded as an attribute.
#' The \code{level_sets} are ordered by their flat index, and named according to their multi-index.
#'
#' If \code{return_reference} is TRUE, the \code{MapperRef} object is returned. This object can be used to fine-tune parts of the mapper construction.
#'
#' @seealso \code{\link[stats]{dist}} \code{\link[stats]{hclust}} \code{\link[Mapper]{cutoff_first_bin}}
#' @useDynLib Mapper
#' @import methods parallelDist htmlwidgets
#' @importFrom Rcpp sourceCpp
#'
#' @author Matt Piekenbrock, \email{matt.piekenbrock@@gmail.com}
#' @encoding UTF-8
#' @references Singh, Gurjeet, Facundo Memoli, and Gunnar E. Carlsson. "Topological methods for the analysis of high dimensional data sets and 3d object recognition." SPBG. 2007.
#' @examples
#'
#' \dontrun{
#' #install.packages("igraph")
#' library(igraph)
#' g1 <- graph.adjacency(m1$adjacency, mode="undirected")
#' plot(g1, layout = layout.auto(g1) )
#' }
#' @export
mapper <- function(X, filter_values, cover_type = "fixed rectangular", return_reference = FALSE, ...) {

  ## Extract the given parameters as a named list.
  extra <- list(...)
  getParam <- function(param, default){  ifelse(is.null(extra[[param]]), default, extra[[param]]) }

  ## Setup 
  if (!is.null(dim(X)) && !methods::is(X, "dist")){ X <- as.matrix(X) } ## convert to matrix
  if (!class(X) %in% c("matrix", "dist")){ stop("Mapper expects 'X' to be either a matrix-coercible data type or a 'dist' object.") }
  m <- MapperRef$new(X = X)
  
  ## Configure mapper 
  m$use_cover(filter_values = filter_values, type = cover_type, getParam("number_intervals", 5L), getParam("percent_overlap", 0.35))$
    use_clustering_algorithm(cl = getParam("cl", "single"))$
    use_distance_measure(measure = getParam("measure", "euclidean"))$
    compute_vertices(num_bins = getParam("num_bins", 10L))$
    compute_edges()

  ## Convert to 'Mapper' object or, if the reference class is wanted, return that. The .summary attribute
  ## stores a useful string used for default printing characteristics of the Mapper.
  # mapperoutput <- if (return_reference) m else m$exportMapper()
  # if (class(mapperoutput) == "Mapper"){
  #   attr(mapperoutput, ".summary") <- c(attr(mapperoutput, ".summary"), paste0("\nCall:  ", paste0(trimws(deparse(match.call())), collapse = " "), sep = ""))
  # }
  return(m)
}

#' S3 method for default printing
#' @param x Mapper object.
#' @param ... unused. 
#' @export
print.Mapper <- function(x, ...){
  writeLines(attr(x, ".summary"))
}

#' S3 method for default plotting
#' @param x Mapper object.
#' @param ... unused.
#' @export
plot.Mapper <- function(x, ...){
  node_sizes <- sapply(x$nodes, length)
  g <- igraph::graph_from_adjacency_matrix(x$adjacency, mode = "undirected")
  xy_coords <- local({ set.seed(1234); igraph::layout.fruchterman.reingold(g) })
  edges <- apply(igraph::as_edgelist(g), 1, function(e){ list(xy_coords[e[1],], xy_coords[e[2],]) })
  edges <- do.call(rbind, lapply(edges, function(lst) do.call(rbind, lst)))
  params <- list(...)
  default_params <- list(xlab=NA, ylab=NA, xaxt='n', yaxt='n', cex=(1 + (node_sizes/sum(node_sizes))), main = attr(x, ".summary")[[1]], pch = 20)
  args <- names(default_params)
  params[args] <- ifelse(args %in% names(params), params[args], default_params[args])
  { params[["x"]] <- xy_coords; do.call(plot, params) }
  { params[["x"]] <- edges; do.call(lines, params) }
}


