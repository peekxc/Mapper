## color_functions.R
## Contains a set of functions that may be used to color nodes
## Each function must accept as input a 'MapperRef' object and return
## as output a vector of numbers of length equal to the number of nodes

## The current environment passed in
c_env <- environment()

## Computes the mean filter value for each node
c_env[["Mean filter value"]] <- function(M, ...){
  agg_pt_fv <- sapply(M$vertices, function(n_idx){ apply(as.matrix(M$cover$filter_values[n_idx,]), 1, mean)})
  agg_node_val <- sapply(agg_pt_fv, mean)
  return(agg_node_val)
}

## Computes the mean data value for each node
c_env[["Mean data value"]] <- function(M, ...){
  agg_pt_x <- sapply(M$vertices, function(n_idx){ apply(as.matrix(M$X[n_idx,]), 1, mean) })
  agg_node_x <- sapply(agg_pt_x, mean)
  return(agg_node_x)
}

## Density of nodes
c_env[["Density"]] <- function(M, ...){
  return(sapply(M$vertices, length))
}

rm(c_env)
