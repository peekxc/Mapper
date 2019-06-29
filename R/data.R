#' 1,500 noisy points sampled on the perimeter of a circle.
#' 
#' The circle has radius ~ 2.1 and is based on example 3.2 of [1].
#' @format A two dimensional of points.
#' @source Generated.
#' @references 
#' Singh, Gurjeet, Facundo Mémoli, and Gunnar E. Carlsson. "Topological methods for the analysis of high dimensional data sets and 3d object recognition." SPBG. 2007.
#' @examples 
#' ## Example code to generate noisy points around the perimeter of a circle 
#' n <- 1500
#' t <- 2*pi*runif(n)
#' r <- runif(n, min = 2, max = 2.1)
#' noisy_circle <- cbind(r*cos(t), r*sin(t))
"noisy_circle"


#' Subset of World Values Survey
#' 
#' This data set contains a small subset of Wave 6 of the World Values Survey. The subset is 
#' specifically only responses collected in the U.S., and contains a fraction of the dimensions, 
#' which were renamed to more semantically clearer names.
#' 
#' The dimensions correspond the following question codes (in order):
#' V4, V5, V7, V8, V9, V10, V11, V23, V24, V25, V29, V30, V32, V55, V56, V57, V58, V59, V60, 
#' V61, V67, V68, V69, V84, V95, V97, V98, V101
#' 
#' Please refer to the WV6 codebook for more details on the questions and responses.
#' 
#' @format A 28-dimensional data.frame of questionnaire responses, unnormalized. 
#' @source \url{http://www.worldvaluessurvey.org/WVSDocumentationWV6.jsp}
"wvs_us_wave6"