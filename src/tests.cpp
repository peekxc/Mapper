#include <Rcpp.h>
using namespace Rcpp;

#include "utility_rcpp.h"


// [[Rcpp::export]]
void run_rcpp_tests() {
  using idx_t = std::vector< std::size_t >;
  auto vec = std::vector< idx_t >();
  std::size_t seq_idx[] = { 0, 1, 2, 3 };
  idx_t t1(seq_idx, seq_idx + 2);
  idx_t t2(seq_idx, seq_idx + 3);
  idx_t t3(seq_idx, seq_idx + 2);
  idx_t t4(seq_idx, seq_idx + 3);
  vec.push_back(t1);
  vec.push_back(t2);
  vec.push_back(t3);
  vec.push_back(t4);
  
  // Regular cartesian product 
  Rcout << "Normal Cartesian product" << std::endl; 
  CartesianProduct(vec, [](const idx_t elem){
    for (std::size_t i = 0; i < elem.size(); ++i){
      Rcout << elem.at(i) << ", ";
    }
    Rcout << std::endl;
  });
  
  
  // Cartesian product with dimensions ordered
  Rcout << "Ordered Cartesian product" << std::endl; 
  idx_t idx_order = { 3, 2, 1, 0 };
  CartesianProductOrdered(vec, [](const idx_t elem){
    for (std::size_t i = 0; i < elem.size(); ++i){
      Rcout << elem.at(i) << ", ";
    }
    Rcout << std::endl;
  }, idx_order);
}




/*** R
## Testing cartesian product 
run_rcpp_tests()

x <- as.matrix(dist(replicate(2, rnorm(10))))
y <- as.matrix(dist(replicate(2, rnorm(10))))
Mapper:::gh_make_Q(x, y)
*/
