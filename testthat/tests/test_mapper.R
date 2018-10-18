## test_mapper.R
context("Testing mapper")
# skip_on_cran()

## Load noisy circle data 
system.file(file.path("data", "noisy_circle.rdata"), package = "Mapper")

## Define filter values equal to the distance from each point to the left-most point in the circle 
left_pt <- noisy_circle[which.min(noisy_circle[, 1]),]
f_x <- matrix(apply(noisy_circle, 1, function(pt) (pt - left_pt)[1]))

## --------- Begin Tests --------- 
## Helper imported: am, check_neighborhoods, check_vertices, check_edges

m <- MapperRef$new(noisy_circle)

## Test can construct mapper 
expect_is(m, "MapperRef")
expect_equal(m$cover, NA)
expect_is(m$simplicial_complex, "Rcpp_SimplexTree")

## Test can construct cover  
cover_params <- list(filter_values = matrix(f_x, ncol = 1), type="fixed rectangular", number_intervals=5L, percent_overlap=0.20)
expect_silent(do.call(m$use_cover, cover_params))
expect_is(m$cover, "CoverRef")

## Test can create clustering algorithm 
cl_params <- list(cl = "single", num_bins = 10L)
expect_silent(do.call(m$use_clustering_algorithm, cl_params))
expect_is(m$clustering_algorithm, "function")

## Test can assign distance measure 
expect_silent(do.call(m$use_distance_measure, list(measure = "euclidean")))
expect_is(m$measure, "character")
  
## Test can compute vertices 
expect_silent(m$compute_vertices())
expect_is(m$vertices, "list")
expect_is(m$ls_vertex_map, "list")

## Test can compute edges 
expect_silent(m$compute_edges())
expect_equal(m$simplicial_complex$n_simplexes, c(8, 8))

## Stricter checks
check_vertices(m)
check_edges(m)


