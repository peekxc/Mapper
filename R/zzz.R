## zzz.R 
## Load the exported modules. 
## The file is named this way to ensure the collation order 

## Simplex Tree module 
Rcpp::loadModule("simplex_tree_module", TRUE)

## Segment Tree module
Rcpp::loadModule("segment_tree_module", TRUE)

## Multiscale Module
Rcpp::loadModule("multiscale_module", TRUE)