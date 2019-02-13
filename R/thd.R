# THD 
#' @export
thd <- function(filter_f, mapper_params, stopping_criterion=list(min_size=10L)){
  
  
  
  ## Get the connected components
  ds <- Mapper::union_find(length(m$vertices))
  el <- m$simplicial_complex$as_edge_list()
  invisible(apply(el, 1, function(e){ ds$union(e[1], e[2]) }))
  cc <- ds$connected_components()
  
  for (cc_id in unique(cc)){
    node_idx <- which(cc == cc_id)
    cc_pt_idx <- unname(unlist(m$vertices[names(m$vertices)[node_idx]]))
  }
  
  shape_mappers[[1]]
  Mapper::union_find()
}


thd_recursive <- function(){
  node$ids <- ids
  
}