## test_uf.R
library("Mapper")
library("testthat")
testthat::context("Testing union find data structure")

## Tests union find data structure
testthat::test_that("UnionFind is constructible", {
  uf <- Mapper::union_find(10)
  testthat::expect_is(uf, "Rcpp_UnionFind")
})


## Tests union find data structure
testthat::test_that("Can search connected components", {
  uf <- Mapper::union_find(10)
  testthat::expect_equal(uf$connected_components(), 0L:9L)
  testthat::expect_equal(uf$find(0), 0L)
  testthat::expect_equal(uf$find_all(0L:9L), 0L:9L)
})

testthat::test_that("Can union components", {
  uf <- Mapper::union_find(10)
  testthat::expect_silent(uf$union(0L, 1L))
  testthat::expect_equal(length(unique(uf$connected_components())), 9L)
  testthat::expect_silent(uf$union_all(0L:4L))
  testthat::expect_silent(uf$union_all(5L:9L))
  testthat::expect_equal(length(unique(uf$connected_components())), 2L)
})

testthat::test_that("Can add components", {
  uf <- Mapper::union_find(10)
  testthat::expect_equal(length(unique(uf$connected_components())), 10L)
  testthat::expect_silent(uf$add_sets(5L))
  testthat::expect_equal(length(unique(uf$connected_components())), 15L)
})

