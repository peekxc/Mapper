#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp14)]] 
#include "../discreture/include/discreture.hpp"
using namespace discreture;

// [[Rcpp::export]]
NumericVector test_discreture(NumericVector x) {
  //using index_t = discreture::combination
  for (auto& x : combinations(6, 2))
  {
    std::cout << x << std::endl; 
  }
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
test_discreture(42)
*/
