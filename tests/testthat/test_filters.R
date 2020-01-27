## test_covers.R
library("Mapper")
library("testthat")
testthat::context("Testing filters")

## test data 
data("wvs_us_wave6", package = "Mapper")
X <- unique(scale(wvs_us_wave6))
X <- X[sample(seq(nrow(X)), size = 250, replace = FALSE),]

## Test assignment of filter
test_that("Can use filters", {
  m <- MapperRef$new(data = X)
  expect_is(m, "MapperRef")
  expect_silent(m$use_filter(filter = X))
  expect_equal(m$filter(), X)
})

## Start with global mapper object
m <- MapperRef$new(data = X)

## Test principle component filter
test_that("Principle Components filter works", {
  expect_silent(m$use_filter(filter = "PC"))
  expect_true(all(dim(m$filter()) == c(nrow(X), 2L)))
})

## Test eccentricity filter
test_that("Eccentricity filter works", {
  expect_silent(m$use_filter(filter = "ECC"))
  f_x <- matrix(colMeans(as.matrix(dist(X))), ncol = 1)
  expect_true(all(abs(m$filter() - f_x) < sqrt(.Machine$double.eps)))
})

## Removed from Suggests: fastICA, ks, TDA, vegan, geigen, umap,

# ## Test independent component filter
# test_that("Independent Components filter works", {
#   expect_silent(m$use_filter(filter = "IC"))
#   expect_true(all(dim(m$filter()) == c(nrow(X), 2L)))
# })
# 
# ## Test kernel density estimate filter
# test_that("KDE filter works", {
#   expect_silent(m$use_filter(filter = "KDE"))
#   f_x <- ks::kde(X, H = diag(apply(X, 2, stats::bw.nrd0)), eval.points = X, verbose = FALSE)$estimate
#   expect_true(all(abs(m$filter() - f_x) < sqrt(.Machine$double.eps)))
# })
# 
# ## Test distance to measure filter
# test_that("DTM filter works", {
#   expect_silent(m$use_filter(filter = "DTM"))
#   f_x <- matrix(do.call(TDA::dtm, list(X=X, Grid=X, m0=0.20, r=2)), ncol = 1L)
#   expect_true(all(abs(m$filter() - f_x) < sqrt(.Machine$double.eps)))
# })
# 
# ## Test distance to measure filter
# test_that("MDS filter works", {
#   expect_silent(m$use_filter(filter = "MDS"))
#   expect_true(all(dim(m$filter()) == c(nrow(X), 2L)))
# })
# 
# ## Test distance to measure filter
# test_that("ISOMAP filter works", {
#   expect_silent(m$use_filter(filter = "ISOMAP"))
#   expect_true(all(dim(m$filter()) == c(nrow(X), 2L)))
# })
# 
# ## Test lacplacian eigenmaps filter
# test_that("LE filter works", {
#   expect_silent(m$use_filter(filter = "LE"))
#   expect_true(all(dim(m$filter()) == c(nrow(X), 2L)))
# })
# 
# ## Test umap filter
# test_that("UMAP filter works", {
#   expect_silent(m$use_filter(filter = "UMAP"))
#   expect_true(all(dim(m$filter()) == c(nrow(X), 2L)))
# })