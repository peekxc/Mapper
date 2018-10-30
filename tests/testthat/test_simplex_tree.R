## test_simplex_tree.R
library("Mapper")
library("testthat")
testthat::context("Testing Simplex Tree")

test_that("Can construct and deconstruct a Simplex Tree object", {
  stree <- Mapper::simplex_tree()
  expect_is(stree, "Rcpp_SimplexTree")
  expect_equal(stree$n_simplexes, numeric(0))
  expect_silent({ rm(stree); invisible(gc(verbose = FALSE)) })
})

test_that("Can add and remove vertices", {
  stree <- Mapper::simplex_tree()
  stree$add_vertices(5L)
  expect_equal(head(stree$n_simplexes, 1), 5L)
  stree$remove_vertices(seq(5))
  expect_equal(head(stree$n_simplexes, 1), 0L)
  rm(stree)
})


test_that("Can add and remove edges", {
  stree <- Mapper::simplex_tree()
  n_vertices <- sample(1:25, size = 1)
  stree$add_vertices(n_vertices)
  edges <- t(combn(n_vertices, 2L))
  
  expect_silent(invisible(apply(edges, 1, function(e){
    stree$insert_simplex(as.integer(e))
  })))
  expect_equal(stree$n_simplexes, c(n_vertices, choose(n_vertices, 2L)))
  
  ## Remove edges, check each time
  cc <- choose(n_vertices, 2L)
  expect_silent(invisible(apply(edges, 1, function(e){
    stree$remove_edge(as.integer(e))
  })))
  expect_equal(stree$n_simplexes, c(n_vertices, 0L))
  rm(stree)
})

test_that("Export types work", {
  stree <- Mapper::simplex_tree()
  n_vertices <- sample(2:25, size = 1)
  stree$add_vertices(n_vertices)
  edges <- t(combn(n_vertices, 2L))
  invisible(apply(edges, 1, function(e){
    stree$insert_simplex(as.integer(e))
  }))
  
  ## Test can export to adjacency matrix 
  expect_is(stree$as_adjacency_matrix(), class = "matrix")
  am <- stree$as_adjacency_matrix()
  expect_equal(nrow(am), n_vertices)
  expect_equal(sum(am == 1), choose(n_vertices, 2)*2) 
  
  ## Test can export to adjacency list 
  expect_is(stree$as_adjacency_list(), class = "list")
  al <- stree$as_adjacency_list()
  expect_equal(length(al), n_vertices)
  expect_true(all(sapply(names(al), function(key) length(al[[key]]) == n_vertices-1)))
  
  ## Test can export to an edgelist
  expect_is(stree$as_edge_list(), class = "matrix")
  el <- stree$as_edge_list()
  expect_equal(dim(el), c(choose(n_vertices, 2), 2L))
})






