## helper_mapper.R
## Helper functions to test validity of mappers

## Short command for getting the adjacency matrix
am <- function(m){ m$simplicial_complex$as_adjacency_matrix() }

## Make sure the two mappers have the same number of simplices, and that their neighborhood statistics are the same
check_neighborhoods <- function(m1, m2){
  require("igraph")
  g1 <- igraph::graph_from_adjacency_matrix(am(m1))
  g2 <- igraph::graph_from_adjacency_matrix(am(m2))
  initial_check <- all(c(
    names(m1$ls_vertex_map) == names(m2$ls_vertex_map),
    m1$simplicial_complex$n_simplexes == m2$simplicial_complex$n_simplexes,
    igraph::vcount(g1) == igraph::vcount(g2),
    igraph::ecount(g1) == igraph::ecount(g2),
    sort(igraph::local_scan(g1, k = 1)) == sort(igraph::local_scan(g2, k = 1)),
    sort(igraph::local_scan(g1, k = 2)) == sort(igraph::local_scan(g2, k = 2)),
    sort(igraph::local_scan(g1, k = 3)) == sort(igraph::local_scan(g2, k = 3))
  ))
}

## Comprehensive vertex check: Ensure that every vertex has at least one exact correspondence 
## between the two mappers
check_vertices <- function(m1, m2){
  f_check <- Vectorize(function(i, j){  
    if (length(m1$vertices[[as.character(i)]]) != length(m2$vertices[[as.character(j)]])){ return(FALSE) }
    all(m1$vertices[[as.character(i)]] %in% m2$vertices[[as.character(j)]]) 
  })
  exact_check <- sapply(names(m1$ls_vertex_map), function(ls_key){
    n_vertices <- length(m1$ls_vertex_map[[ls_key]])
    if (n_vertices == 0){ return(n_vertices == length(m2$ls_vertex_map[[ls_key]])) }
    else {
      sum(outer(m1$ls_vertex_map[[ls_key]], m2$ls_vertex_map[[ls_key]], f_check)) >= n_vertices
    }
  })
  return(exact_check)
}

## Comprehensive edge check: Ensure every connected and non-connected edge has a non-empty and empty 
## intersection, respectively
check_edges <- function(m){
  rowmatch <- function(A,B) {# Rows in A that match the rows in B
    f <- function(...) paste(..., sep=":")
    if(!is.matrix(B)) B <- matrix(B, 1, length(B))
    a <- do.call("f", as.data.frame(A))
    b <- do.call("f", as.data.frame(B))
    match(a, b)
  }
  connected_edges <- m$simplicial_complex$as_edge_list()
  all_edges <- apply(t(combn(length(m$vertices), 2)), 2, function(idx) names(m$vertices)[idx])
  conn_idx <- rowmatch(apply(connected_edges, 2, as.character), all_edges)
  
  ## Ensures edges recorded as connected have a non-empty intersection
  connected_check <- all(apply(connected_edges, 1, function(x){
    length(intersect(m$vertices[[as.character(x[1])]], m$vertices[[as.character(x[2])]])) > 0
  }))
  
  ## Ensures edges recorded as unconnected have an empty intersection
  unconnected_check <- all(apply(matrix(all_edges[-conn_idx,], ncol = 2), 1, function(x){
    length(intersect(m$vertices[[as.character(x[1])]], m$vertices[[as.character(x[2])]])) == 0
  }))
  
  ## Return 
  return(connected_check && unconnected_check)
}