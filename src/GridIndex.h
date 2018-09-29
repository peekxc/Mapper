#include <Rcpp.h>
using namespace Rcpp;

#include "cartesian_product.h"

// = GridIndex = 
// GridIndex is a simple, integral-templated indexing structure that allows a quick mapping between array multi-indexes to their flat equivalents, 
// and vice-versa. The grid is built using a fast stack-based iterative cartesian product to populate the array structure.
template<typename T>
struct GridIndex {
  static_assert(std::is_integral<T>::value, "Integral-type required as an index.");
  const std::size_t d;                    // Dimension of the multi-index
  std::size_t n_multi_indices;            // Total number of multi-indexes indices 
  std::vector< T > index_AoS;             // Implicitly represents an AoS layout
  std::vector< std::size_t > dim_sizes;   // Sizes of each dimension
  std::vector< std::size_t > cum_prod;    // Precomputed cumulative product to save on indexing 
  
  GridIndex(const IntegerVector idx);
  
  T flat_from_multi(std::vector<T> index);
  std::vector<T> multi_from_flat(std::size_t index);
    
  // Wrappers for interfacing from R for testing
  int flat_from_multi_rcpp(IntegerVector index);
  IntegerVector multi_from_flat_rcpp(std::size_t index);
  IntegerMatrix multi_matrix();
};