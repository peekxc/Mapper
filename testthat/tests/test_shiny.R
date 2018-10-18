## test_shiny.R

## Load the libraries 
library("shiny")
library("keras")
library("htmlwidgets")
library("jsonlite")
library("grapher")
library("Mapper")

## Generate noisy points around the perimeter of a circle 
n <- 150
t <- 2*pi*runif(n)
r <- runif(n, min = 2, max = 2.1)
circ <- cbind(r*cos(t), r*sin(t))

## Define filter values equal to the distance from each point to the left-most point in the circle 
left_pt <- circ[which.min(circ[, 1]),]
f_x <- sqrt(colSums(apply(circ, 1, function(pt) (pt - left_pt)^2)))

## Make the mapper
m <- Mapper::mapper(X = circ, filter_values = f_x, number_intervals = 5, overlap = 0.20)

## Make sure it works 
library("R.utils")
R.utils::evalWithTimeout(shinyApp(ui = ui, server = server), timeout = 1, onTimeout = "error")

library("shinytest")
shinytest::recordTest(shinyApp(ui = ui, server = server))