#include "Rcpp.h"
using namespace Rcpp;

#include <limits>
#include <algorithm>

// [[Rcpp::export]]
NumericMatrix nondecreasing_matrix(const NumericMatrix& mat){
  const std::size_t n = mat.nrow(), d = mat.ncol();
  NumericMatrix res = no_init_matrix(n*d, d);
  NumericVector c_row = NumericVector(d, 0);
  std::vector< std::size_t > idx = std::vector< std::size_t >(d, 0);
  for (std::size_t i = 0; i < n*d; ++i){
    
    // Extract min dimension + value
    double min_val = std::numeric_limits<double>::infinity(); 
    std::size_t min_idx = 0, d_i = 0; 
    std::for_each(idx.begin(), idx.end(), [&](const std::size_t ii){
      double tmp = mat(i, ii);
      if (tmp < min_val && idx.at(d_i) <= n){ 
        min_idx = ii;
        min_val = tmp; 
      }
      d_i++;
    });
    
    // Assign 
    c_row.at(min_idx) = min_val;
    res(i, _) = c_row; 
    
    // Increment 
    idx.at(min_idx) = idx.at(min_idx) + 1;
  }
  return(res);
}
