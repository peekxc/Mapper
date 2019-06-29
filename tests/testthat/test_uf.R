## test_uf.R
library("Mapper")
library("testthat")
testthat::context("Testing union find data structure")

## Tests union find data structure
testthat::test_that("UnionFind is constructible", {
  uf <- Mapper::union_find(10)
  testthat::expect_is(uf, "Rcpp_UnionFind")
})