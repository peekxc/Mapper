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
m <- MapperRef$new(noisy_circle)

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
  expect_is(m$measure, "proxy_registry_entry")
})

test_that("Can construct cover", {
  expect_silent(m$construct_cover())
  expect_true(length(m$cover$sets) > 0)
})
  
test_that("Can construct pullback", {
  expect_silent(m$construct_pullback())
  expect_is(m$vertices, "list")
  expect_is(m$pullback, "list")
  for (pid in names(m$pullback)){
    pt_ids <- unname(unlist(m$vertices[ as.character(m$pullback[[pid]]) ]))
    expect_equal(
      sort(pt_ids), 
      sort(as.vector(unlist(m$cover$construct(filter = m$filter, index = pid)), mode="integer"))
    )
  }
})

test_that("Can construct nerve", {
  expect_silent(m$construct_nerve(k = 0L))
  expect_equal(m$simplicial_complex$n_simplices, c(8))
  expect_silent(m$construct_nerve(k = 1L))
  expect_equal(m$simplicial_complex$n_simplices, c(8, 8))
})

## Stricter check that tests every connected and non-connected edge
test_that("Mapper is valid", {
  expect_true(check_edges(m))
})

## Check distance matrix works
test_that("Mapper can be constructed with a distance matrix", {
  left_pt <- noisy_circle[which.min(noisy_circle[, 1]),]
  f_x <- matrix(apply(noisy_circle, 1, function(pt) (pt - left_pt)[1]))
  m1 <- MapperRef$new()$
    use_data(dist(noisy_circle))$
    use_filter(f_x)$
    use_cover(cover="fixed interval", number_intervals=5L, percent_overlap=50)$
    use_clustering_algorithm(cl="single", cutoff_method = "continuous", threshold = 0.0)$
    construct_k_skeleton(k = 1L)
  expect_true(check_edges(m1))
  m2 <- MapperRef$new()$
    use_data(noisy_circle)$
    use_filter(f_x)$
    use_cover(cover="fixed interval", number_intervals=5L, percent_overlap=50)$
    use_distance_measure("euclidean")$
    use_clustering_algorithm(cl="single", cutoff_method = "continuous", threshold = 0.0)$
    construct_k_skeleton(k = 1L)
  testthat::expect_equivalent(m1$exportMapper(), m2$exportMapper())
})

## Check wrapper works 
test_that("Mapper wrapper matches MapperRef object", {
  left_pt <- noisy_circle[which.min(noisy_circle[, 1]),]
  f_x <- matrix(apply(noisy_circle, 1, function(pt) (pt - left_pt)[1]))
  m1 <- mapper(data = noisy_circle, 
               filter = f_x,
               cover = list(cover="fixed interval", number_intervals=5L, percent_overlap=50),
               distance_measure = "euclidean",
               clustering_algorithm = list(cl="single", cutoff_method = "continuous", threshold = 0.0), return_reference = TRUE)
  m2 <- MapperRef$new()$
    use_data(noisy_circle)$
    use_filter(f_x)$
    use_cover(cover="fixed interval", number_intervals=5L, percent_overlap=50)$
    use_distance_measure("euclidean")$
    use_clustering_algorithm(cl="single", cutoff_method = "continuous", threshold = 0.0)$
    construct_k_skeleton(k = 1L)
  testthat::expect_equal(m1, m2)
})




