#' @title mapper
#'
#' @description Computes the mapper graph. Mapper is a tool for summarizing topological information from datasets and
#' maps defined on them. It takes as input a set of 'point cloud' data, a (possibly lower dimensional) map defined on the data,
#' and produces a topological summary of the data expressed through a cover equipped to the codomain of the map.
#'
#' @param X Either an \eqn{n x D} data matrix.
#' @param filter An \eqn{n x d} data matrix, or a function.
#' @param cover_params Named list of cover parameters. See details. 
#' @param measure String indicating the measure in the data space. Accepts any in \link[proxy]{pr_DB}. 
#' @param cluster_params Named list of clustering parameters. See details. 
#' @param return_reference Boolean whether or not to return the reference class used to construct Mapper. See \link[Mapper]{MapperRef}.
#' @details
#' \code{mapper} is a generic function that concisely parameterizes the Mapper framework into a single function definition.
#' This function serves as a convenience wrapper around the \link[Mapper]{MapperRef} R6 generator for users used to functional programming. 
#' For finer control over the mapper construction, use \link[Mapper]{MapperRef} instead. If \code{return_reference} is TRUE, the 
#' \link[Mapper]{MapperRef} instance used by this function is returned. 
#' 
#' The \code{cover_params} must be a named list of all of the parameters needed by the cover. 
#' It must have a \code{typename} member matching one of the typenames listed in \link[Mapper]{covers_available}, along with all of the 
#' corresponding parameters for that cover type.   
#' 
#' The \code{cluster_params} must be a named list containing members \code{cl} (string) and \code{num_bins} (integer), where 
#' \code{cl} is one of the linkage criterion used by \link[stats]{hclust}, and \code{num_bins} corresponds the binning argument in 
#' the cutting heuristic \link[Mapper]{cutoff_first_bin}.  
#'
#' @return 
#' By default, returns a \code{mapper} object, which is a named list with members \code{adjacency} (adjacency matrix for the edges),
#' \code{vertices} (list with point membership indices), and list \code{level_sets} (list with vertex indices).
#' Each vertex in \code{vertices} also records the flat index of the level set that vertex intersects, recorded as an attribute.
#'
#' If \code{return_reference} is TRUE, the \code{MapperRef} object is returned. This object can be used to fine-tune parts of the mapper construction.
#'
#' @seealso \code{\link[stats]{hclust}} \code{\link[Mapper]{cutoff_first_bin}}
#' @useDynLib Mapper
#' @import methods
#' @importFrom Rcpp sourceCpp
#'
#' @author Matt Piekenbrock, \email{matt.piekenbrock@@gmail.com}
#' @encoding UTF-8
#' @references Singh, Gurjeet, Facundo Memoli, and Gunnar E. Carlsson. "Topological methods for the analysis of high dimensional data sets and 3d object recognition." SPBG. 2007.
#' @examples
#' data("noisy_circle", package="Mapper")
#' left_pt <- noisy_circle[which.min(noisy_circle[, 1]),]
#' f_x <- matrix(apply(noisy_circle, 1, function(pt) (pt - left_pt)[1]))
#' 
#' m <- mapper(X = noisy_circle, filter_values = f_x, 
#'             cover_params = list(typename="fixed interval", number_intervals=10L, percent_overlap=50),
#'             measure = "euclidean", 
#'             cluster_params = list(cl="single", threshold = 0.0))
#' @export
mapper <- function(X, filter_values, 
                   cover_params = c(typename="fixed interval", number_intervals=10L, percent_overlap=35), 
                   measure = "euclidean",
                   cluster_params = c(cl="single"),
                   return_reference = FALSE) {
  ## Extract the given parameters as a named list.
  # extra <- list(...)
  # getParam <- function(param, default){  ifelse(is.null(extra[[param]]), default, extra[[param]]) }

  ## Setup 
  if (!is.null(dim(X))){ X <- as.matrix(X) } ## convert to matrix
  if (!class(X) %in% c("matrix")){ stop("Mapper expects 'X' to be either a matrix-coercible data type.") }
  m <- MapperRef$new(X = X)
  
  ## Configure mapper 
  m$use_filter(filter = filter_values)
  do.call(m$use_cover, cover_params)
  do.call(m$use_distance_measure, list(measure=measure))
  do.call(m$use_clustering_algorithm, cluster_params)
  m$construct_k_skeleton(k = 1L)

  ## Convert to 'Mapper' object or, if the reference class is wanted, return that. The .summary attribute
  ## stores a useful string used for default printing characteristics of the Mapper.
  if (return_reference){
    return(m)
  } else {
    mapperoutput <- m$exportMapper()
    # if (class(mapperoutput) == "Mapper"){
    #   attr(mapperoutput, ".summary") <- c(attr(mapperoutput, ".summary"), paste0("\nCall:  ", paste0(trimws(deparse(match.call())), collapse = " "), sep = ""))
    # }
    rm(m)
    return(mapperoutput)
  }
}

#' S3 method for default printing
#' @param x Mapper object.
#' @param ... unused. 
#' @export
print.Mapper <- function(x, ...){
  writeLines(attr(x, ".summary"))
}

# Empty environment to allow passing parameters to dashboard
.dash_env <- new.env(parent = emptyenv())


