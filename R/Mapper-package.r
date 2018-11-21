#' Mapper Package
#' @description Provides an R/Rcpp implementation of Mapper and its related tools. 
#' The main documentation for the package can be found here: \url{https://peekxc.github.io/mapper/}.  
#' @author Matt Piekenbrock
#' @name Mapper
#' @import Rcpp
#' @importFrom Rcpp evalCpp Module cpp_object_initializer
#' @importFrom methods new
#' @useDynLib Mapper, .registration = TRUE
#' @docType package
NULL


#' bin_color
#' @param x A numeric vector whose magnitudes should be binned onto the color palette. 
#' @param col_pal Color palette to bin numeric values into. See details. 
#' @param output_format Whether the output should be RGBa (hex9) or RGB (hex7).
#' @description The \code{col_pal} must be a string vector of hexadecimal RGB or RGBa color codes.  
#' @details Given a numeric vector \code{x}, bins the values from low to high on a given color gradient. 
#' Defaults to the reversed rainbow gradient, where blue == low, red == high.
#' @export
bin_color <- function(x, col_pal = "rainbow", 
                     output_format = c("hex9", "hex7")){
  if (missing(col_pal) || col_pal == "rainbow"){ col_pal <- rev(grDevices::rainbow(100L, start = 0, end = 4/6)) }
  col_res <- length(col_pal)
  binned_idx <- cut(x, breaks = col_res, labels = FALSE)
  binned_colors <- col_pal[binned_idx]
  if (missing(output_format) || output_format == "hex9"){ return(binned_colors) }
  else if (output_format == "hex7"){ return(substr(binned_colors, start = 0L, stop = 7L)) }
  else { stop("'output_format' must be one of 'hex9' or 'hex7'.") }
}