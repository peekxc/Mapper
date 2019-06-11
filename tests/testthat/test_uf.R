## test_uf.R
library("Mapper")
library("testthat")
testthat::context("Testing mapper")

## Tests union find data structure
testthat::test_that("UnionFind is constructible", {
  uf <- Mapper::union_find(10)
  testthat::expect_is(uf, "Rcpp_UnionFind")
})