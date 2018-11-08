## Sources.R
## Attaches the default sources to the Mapper dashboard app.
## Author: Matt Piekenbrock

## Normalize between 0-1, unless all the same
normalize <- function(x) { 
  if (all(x == x[1])){ return(rep(1, length(x))) }
  else {  (x - min(x))/(max(x) - min(x)) }
}

mapper_graph <- reactive({
  am <- M$.simplicial_complex$as_adjacency_matrix()
  igraph::graph_from_adjacency_matrix(am, mode = "undirected", add.colnames = NA) 
})

## Makes the default configuration depend on M
json_config <- reactive({
  MG <- mapper_graph()
  grapher::getDefaultJsonConfig(network=MG)
})

# update_color <- reactive({
#   ## Color nodes and edges by a default rainbow palette
#   rbw_pal <- rev(rainbow(100, start = 0, end = 4/6))
#   agg_node_val <- sapply(sapply(M$vertices, function(v_idx){ 
#     apply(as.matrix(M$cover$filter_values[v_idx,]), 1, mean)
#   }), mean)
#   binned_idx <- cut(agg_node_val, breaks = 100, labels = F)
#   vertex_colors <- rbw_pal[binned_idx]
#   json_config()$nodes$color <- grapher::hex2rgba(vertex_colors)
# })
# 
# vertex_sizes <- sapply(M$vertices, length) 

# 
#   
#   ## Create the vertices. By default, color them on a rainbow palette according to their mean filter value.
#   if (length(igraph::V(G)) > 0){
#     vertex_radii <- (15L - 10L)*normalize(log(vertex_sizes)) + 10L
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
# })



