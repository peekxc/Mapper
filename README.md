
# Mapper 
This package provides a set of tools written in R/Rcpp for computing the _mapper_ construction, and other related algorithms. Mapper was originally introduced in the article below: 

> Singh, Gurjeet, Facundo Mémoli, and Gunnar E. Carlsson. "Topological methods for the analysis of high dimensional data sets and 3d object recognition." SPBG. 2007.

## Installation 

The current development version can be installed with the [devtools](https://github.com/r-lib/devtools) package: 
```R
require("devtools")
devtools::install_gitub("peekxc/mapper")
```

A stable CRAN release is planned for the future. 

## Usage

Given a data set, define the filter function. Here is an example using the noisy points sampled from the perimeter of a circle, similar to the example given by Example 3.2 in the original paper.   
```R
data("noisy_circle")

## Define filter values equal to the distance from each point to the left-most point in the circle 
left_pt <- noisy_circle[which.min(noisy_circle[, 1]),]
f_x <- matrix(apply(noisy_circle, 1, function(pt) (pt - left_pt)[1]))
```

You can construct a __mapper__ with [R6 method chaining](https://adv-r.hadley.nz/r6.html#method-chaining)
```R
## Define the main via chaining R6 methods
m <- MapperRef$new(noisy_circle)$
  use_cover(filter_values = matrix(f_x), type="fixed rectangular", number_intervals=5L, percent_overlap=20)$
  use_clustering_algorithm(cl = "single", num_bins = 10)$
  use_distance_measure(measure = "euclidean")$
  compute_k_skeleton(k=1L)
```

You can export to your favorite graph-based representation. 
```R
## Overview of the simplicial complex 
print(m$simplicial_complex)

## More detailed summary 
m$simplicial_complex$print_tree()

## Export 1-skeleton
am <- m$simplicial_complex$as_adjacency_matrix()
```

The vertices of the __mapper__ are stored as a simple list 
```R
View(m$vertices)
```

To view the graph interactively, consider the [grapher](https://github.com/peekxc/grapher) library. 
```R
library("grapher")
m$as_grapher() 
```

## Additional Information 

More comprehensive documentation is available [here](https://peekxc.github.io/mapper/)
See the more extensive [vignette on using the package](https://peekxc.github.io/mapper/articles/UsingMapper.html)

