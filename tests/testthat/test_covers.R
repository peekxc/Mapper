## test_covers.R
library("Mapper")
library("testthat")
testthat::context("Testing covers")

## test data 
set.seed(1234)
f_1 <- function(){ replicate(1, runif(150)) }
f_2 <- function(){ replicate(2, runif(150)) }

## Fixed Interval cover testing
test_that("Can construct fixed interval covers", {
  expect_is(FixedIntervalCover$new(), "CoverRef")
  expect_is(FixedIntervalCover$new(), "FixedIntervalCover")
  cover <- FixedIntervalCover$new(number_intervals = 5, percent_overlap = 25)
  expect_is(cover, "FixedIntervalCover")
  
  ## Test 1D 
  expect_silent(cover$construct_cover(filter = f_1))
  expect_equal(cover$index_set, names(cover$level_sets))
  expect_equal(length(unname(unique(unlist(cover$level_sets)))), 150L)
  
  ## Test 2D 
  expect_silent(cover$construct_cover(filter = f_2))
  expect_equal(cover$index_set, names(cover$level_sets))
  expect_equal(length(unname(unique(unlist(cover$level_sets)))), 150L)
})
  

## Restrained Interval cover testing
test_that("Can construct restrained interval covers", {
  expect_is(RestrainedIntervalCover$new(), "CoverRef")
  expect_is(RestrainedIntervalCover$new(), "RestrainedIntervalCover")
  cover <- RestrainedIntervalCover$new(number_intervals = 5, percent_overlap = 25)
  expect_is(cover, "RestrainedIntervalCover")
  
  ## Test 1D 
  expect_silent(cover$construct_cover(filter = f_1))
  expect_equal(cover$index_set, names(cover$level_sets))
  expect_equal(length(unname(unique(unlist(cover$level_sets)))), 150L)
  
  ## Test 2D 
  expect_silent(cover$construct_cover(filter = f_2))
  expect_equal(cover$index_set, names(cover$level_sets))
  expect_equal(length(unname(unique(unlist(cover$level_sets)))), 150L)
})

## Restrained Interval cover testing
test_that("Can construct ball covers", {
  expect_is(BallCover$new(), "CoverRef")
  expect_is(BallCover$new(), "BallCover")
  cover <- BallCover$new(epsilon = 0.10)
  expect_is(cover, "BallCover")
  
  ## Test 1D 
  cover$epsilon <- 0.15*max(dist(f_1()))
  expect_silent(cover$construct_cover(filter = f_1))
  expect_equal(cover$index_set, names(cover$level_sets))
  expect_equal(length(unname(unique(unlist(cover$level_sets)))), 150L)
  
  ## Test 2D 
  cover$epsilon <- 0.15*max(dist(f_2()))
  expect_silent(cover$construct_cover(filter = f_2))
  expect_equal(cover$index_set, names(cover$level_sets))
  expect_equal(length(unname(unique(unlist(cover$level_sets)))), 150L)
})


test_that("Can construct mappers with all covers", {
  m <- MapperRef$new()
  m$use_data(f_2())
  m$use_filter(f_1())
  
  ## Fixed interval cover
  m$use_cover("fixed interval", number_intervals=10, percent_overlap=20)
  expect_silent(m$construct_k_skeleton(k=1))
  expect_true(sum(m$simplicial_complex$n_simplices) > 0)
  
  ## Restrained interval cover
  m$use_cover("restrained interval", number_intervals=10, percent_overlap=20)
  expect_silent(m$construct_k_skeleton(k=1))
  expect_true(sum(m$simplicial_complex$n_simplices) > 0)
  
  ## Ball cover
  m$use_filter(f_2())
  m$use_cover("ball", epsilon = max(dist(m$filter()))/20)
  expect_silent(m$construct_k_skeleton(k=1))
  expect_true(sum(m$simplicial_complex$n_simplices) > 0)
})



