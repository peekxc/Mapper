#include <Rcpp.h>
using namespace Rcpp;


#include "cartesian_product.h"

// template<typename T> 
// struct index_struct {
//   std::vector< std::vector<T> > index_SoA;
//   index_struct(const IntegerVector& idx){
//     
//   }
//   
// };

// [[Rcpp::export]]
SEXP index_test(IntegerVector lsmi) {
  return(self_match(lsmi));
}


/*** R
x <- sample(1L:5L, size = 15, replace = TRUE)
index_test(x)
*/
