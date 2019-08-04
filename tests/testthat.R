library(testthat)
library(Mapper)

## To test on rhub: 
# rhub::check(path = ".", platform = unlist(rhub:::check_shortcut_platforms[c("linux", "windows", "macos")]), email = "matt.piekenbrock@gmail.com", check_args = "--no-manual")

## Using testthat
test_check("Mapper")
