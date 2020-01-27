#' Arbitrary Cover
#' @docType class
#' @description This class enables the user to use any fixed arbitrary cover. 
#' To facilitate this, rather than providing a parameterized family of covers as the others
#' do (see e.g. \code{\link{FixedIntervalCover}}), the user provides the indexed set of 
#' open sets themselves.
#' 
#' @author Matt Piekenbrock
#' @family cover
#' @export
ArbitraryCover <- R6::R6Class(
  classname = "ArbitraryCover",
  inherit = CoverRef,
  public = list(method="euclidean")
)

## initialize ------
ArbitraryCover$set("public", "initialize", function(...){
  super$initialize(typename="arbitrary")
  params <- list(...)
  if ("level_sets" %in% names(params)){ self$level_sets <- params[["level_sets"]] }
  if ("method" %in% names(params)){ self$method <- params[["method"]] }
})

## validate ------
ArbitraryCover$set("public", "validate", function(filter){
  stopifnot(!is.null(self$epsilon))
  stopifnot(is.character(self$method) || is.function(self$method))
})

## format ----
ArbitraryCover$set("public", "format", function(...){
  titlecase <- function(x){
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }
  sprintf("%s Cover: (epsilon = %.2f)", titlecase(private$.typename), self$epsilon)
})

## construct_cover ------
ArbitraryCover$set("public", "construct_cover", function(filter, index=NULL){
  self$validate()
  
  ## Return specific level sets if requested
  if (!missing(index)){ return(self$level_sets[[index]]) }
  
  ## Always return self 
  invisible(self)
})

ArbitraryCover$set("public", "plot", function(filter, add=FALSE, index=NULL, color_pal=NULL, ...){
  
  ## Get filter values 
  fv <- filter()
  f_dim <- ncol(fv)
  f_size <- nrow(fv)
  
  if (missing(add) || !add){
    set_rng <- range(fv)
    graphics::plot.new()
    graphics::plot.window(xlim=set_rng[1], ylim=c(0,n_overlap+1))
  }
  
  
  
})
