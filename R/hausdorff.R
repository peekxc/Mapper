#' Hausdorff Distance
#' @description Computes a smooth approximation to the Hausdorff distance between the vertices of a mapper.
#' Implements the equation in Section 5 of the original Mapper paper. 
#' \deqn{
#'  d_h(X_i, X_j) = \max{( \frac{\sum_y \min_x(d(x, y)) }{n_j}, \frac{\sum_x \min_y(d(x, y)) }{n_i})}
#' }{
#'  DH(Xi, Xj) = max\{ \sum min(d(x, y)) / nj, \sum min(d(x, y)) / ni\}
#' }
#' where \eqn{x, y} are elements of the clusters \eqn{C_i, C_j}{Ci, Cj} respectively, and \eqn{X_i, X_j}{Xi, Xj} are the clusters corresponding 
#' vertices in the mapper. See the reference below for more details. \cr
#' \cr
#' \strong{WARNING:} This function may be very computationally expensive. 
#' @param m A \code{MapperRef} object. 
#' @return A \code{\link[stats]{dist}} object with \code{Size} equal to the number of vertices. 
#' @details Currently, this requires the \pkg{RANN} to be installed, and only considers euclidean distance.
#' @references Singh, Gurjeet, Facundo Memoli, and Gunnar E. Carlsson. "Topological methods for the analysis of high dimensional data sets and 3d object recognition." SPBG. 2007.
#' @export
hausdorff_distance <- function(m){
  ## Smooth approximation of hausdorff distance
  n_vertex <- length(m$vertices)
  stopifnot(is.numeric(n_vertex))
  stopifnot(n_vertex > 0)
  idx <- utils::combn(n_vertex, 2)
  d_h.vec <- apply(idx, 2, function(ii){
    x1 <- m$X(m$vertices[[ii[1]]])
    x2 <- m$X(m$vertices[[ii[2]]])
    xy <- mean(RANN::nn2(data = x1, query = x2, k = 1)$nn.dist)
    yx <- mean(RANN::nn2(data = x2, query = x1, k = 1)$nn.dist)
    max(xy, yx)
  })
  d_h.matrix <- matrix(0, nrow=n_vertex, ncol=n_vertex)
  d_h.matrix[t(idx)] <- d_h.vec
  d_h <- stats::as.dist(t(d_h.matrix))
  attr(d_h, "method") <- "hausdorff"
  return(d_h)
}
