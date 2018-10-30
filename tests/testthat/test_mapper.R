## test_mapper.R
library("Mapper")
library("testthat")
testthat::context("Testing mapper")
# skip_on_cran()

## Load noisy circle data 
system.file(file.path("data", "noisy_circle.rdata"), package = "Mapper")

## Define filter values equal to the distance from each point to the left-most point in the circle 
left_pt <- noisy_circle[which.min(noisy_circle[, 1]),]
f_x <- matrix(apply(noisy_circle, 1, function(pt) (pt - left_pt)[1]))

## --------- Begin Tests --------- 
## Helper imported: am, check_neighborhoods, check_vertices, check_edges

test_that("Can construct MapperRef object", {
  m <- MapperRef$new(noisy_circle)
  expect_is(m, "MapperRef")
  expect_equal(m$cover, NA)
  expect_is(m$simplicial_complex, "Rcpp_SimplexTree")
})

## Let mapper object be global now
m <- MapperRef$new(noisy_circle)

test_that("Can construct CoverRef object", {
  cover_params <- list(filter_values = matrix(f_x, ncol = 1), type="fixed rectangular", number_intervals=5L, percent_overlap=20)
  expect_silent(do.call(m$use_cover, cover_params))
  expect_is(m$cover, "CoverRef")
})

test_that("Can create clustering algorithm", {
  cl_params <- list(cl = "single", num_bins = 10L)
  expect_silent(do.call(m$use_clustering_algorithm, cl_params))
  expect_is(m$clustering_algorithm, "function") 
})

test_that("Can assign distance measure", {
  expect_silent(do.call(m$use_distance_measure, list(measure = "euclidean")))
  expect_is(m$measure, "character")
})
  
test_that("Can compute vertices", {
  expect_silent(m$compute_vertices())
  expect_is(m$vertices, "list")
  expect_is(m$ls_vertex_map, "list")
})

test_that("Can compute edge", {
  expect_silent(m$compute_edges())
  expect_equal(m$simplicial_complex$n_simplexes, c(8, 8))
})

## Stricter check that tests every connected and non-connected edge
test_that("Mapper is valid", {
  expect_true(check_edges(m))
})


