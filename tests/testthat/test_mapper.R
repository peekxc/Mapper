## test_mapper.R
library("Mapper")
library("testthat")
testthat::context("Testing mapper")
# skip_on_cran()

## Load noisy circle data 
data("noisy_circle", package = "Mapper")

## --------- Begin Tests --------- 
## Helper imported: am, check_neighborhoods, check_vertices, check_edges

test_that("Can construct MapperRef object", {
  m <- MapperRef$new(noisy_circle)
  expect_is(m, "MapperRef")
  expect_is(m$simplicial_complex, "Rcpp_SimplexTree")
})

## Let mapper object be global now
m <- MapperRef$new(X = noisy_circle)

## Define filter values equal to the distance from each point to the left-most point in the circle 
test_that("Can assign filter values", {
  left_pt <- noisy_circle[which.min(noisy_circle[, 1]),]
  f_x <- matrix(apply(noisy_circle, 1, function(pt) (pt - left_pt)[1]))
  expect_silent(m$use_filter(filter = f_x))
  expect_true(is.function(m$filter))
  expect_is(m$filter(), "matrix")
})

## Test cover construction
test_that("Can construct CoverRef object", {
  cover_params <- list(cover="fixed interval", number_intervals=5L, percent_overlap=20)
  expect_silent(do.call(m$use_cover, cover_params))
  expect_is(m$cover, "CoverRef")
})

test_that("Can create clustering algorithm", {
  cl_params <- list(cl = "single", threshold = 0.0)
  expect_silent(do.call(m$use_clustering_algorithm, cl_params))
  expect_is(m$clustering_algorithm, "function") 
})

test_that("Can assign distance measure", {
  expect_silent(do.call(m$use_distance_measure, list(measure = "euclidean")))
  expect_is(m$measure, "character")
})
  
test_that("Can construct pullback", {
  expect_silent(m$construct_pullback())
  expect_is(m$vertices, "list")
  expect_is(m$pullback, "list")
  for (pid in names(m$pullback)){
    pt_ids <- unname(unlist(m$vertices[ as.character(m$pullback[[pid]]) ]))
    expect_equal(sort(pt_ids), sort(m$cover$construct_cover(filter = m$filter, index = pid)))
  }
})

test_that("Can construct nerve", {
  expect_silent(m$construct_nerve(1L))
  expect_equal(m$simplicial_complex$n_simplices, c(8, 8))
})

## Stricter check that tests every connected and non-connected edge
test_that("Mapper is valid", {
  expect_true(check_edges(m))
})


