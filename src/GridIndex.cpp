
#include "GridIndex.h"

template <typename T>
GridIndex<T>::GridIndex(const IntegerVector idx) : d(idx.size()) {
  std::size_t d_i = 0;
  n_multi_indices = std::accumulate(idx.begin(), idx.end(), 1, std::multiplies<std::size_t>());
  
  // Cumulative product (backwards) to precompute offsets
  cum_prod = std::vector< std::size_t >(d, 0);
  for (d_i = d; d_i > 0; --d_i){
    cum_prod.at(d_i-1) = std::accumulate(idx.begin(), idx.begin()+d_i-1, 1, std::multiplies<std::size_t>());
    // Rprintf("d: %d := %d\n", d_i-1, cum_prod.at(d_i-1));
  }
  
  // Expand the given dimensions as sequences starting at 0
  std::vector< std::vector<T> > index_rngs = std::vector< std::vector<T> >(idx.size());
  std::for_each(idx.begin(), idx.end(), [&](const int max_i){
    std::vector<T> rng = std::vector<T>(max_i);
    std::iota(rng.begin(), rng.end(), 0);
    dim_sizes.push_back(max_i);
    index_rngs.at(d_i++) = rng;
  });
  
  // Store the cartesian product of these ranges
  index_AoS.reserve(std::accumulate(dim_sizes.begin(), dim_sizes.end(), 1, std::multiplies<std::size_t>()));
  CartesianProduct(index_rngs, [&](const std::vector<T> c_index){
    index_AoS.insert(index_AoS.end(), c_index.begin(), c_index.end());
  });
}

template <typename T>
std::string GridIndex<T>::multi_to_string(std::vector<T> index){
  std::stringstream result;
  std::copy(index.begin(), index.end(), std::ostream_iterator< int >(result, " "));
  std::string key = result.str();
  return(key.substr(0, key.length() - 1));
}

template <typename T>
T GridIndex<T>::flat_from_multi(std::vector<T> index){
  if (index.size() != d){ return(T(-1)); }
  std::size_t flat_index = 0, d_i = d;
  std::for_each(index.rbegin(), index.rend() - 1, [&](T dim_idx){
    flat_index += cum_prod.at(--d_i)*dim_idx;
  });
  flat_index += index.front();
  return(flat_index);
}

template <typename T>
int GridIndex<T>::flat_from_multi_rcpp(IntegerVector index){
  std::vector<T> tmp(index.begin(), index.end());
  return(static_cast<int>(flat_from_multi(tmp)));
}

template <typename T>
std::vector<T> GridIndex<T>::multi_from_flat(std::size_t index){
  if (index < 0 || index >= n_multi_indices) { stop("Flat index out of range."); }
  std::vector<T> res(d);
  std::size_t offset = (index*d);
  std::copy(index_AoS.begin() + offset, index_AoS.begin() + offset + d, res.begin());
  return(res);
}

template <typename T>
IntegerVector GridIndex<T>::multi_from_flat_rcpp(std::size_t index){
  std::vector<T> multi_index = multi_from_flat(index);
  IntegerVector res = IntegerVector(multi_index.begin(), multi_index.end());
  return(res);
}

template <typename T>
IntegerMatrix GridIndex<T>::multi_matrix(){
  const std::size_t n_rows = std::accumulate(dim_sizes.begin(), dim_sizes.end(), 1, std::multiplies<std::size_t>());
  IntegerMatrix res = no_init_matrix(n_rows, d);
  for (std::size_t i = 0; i < n_rows; ++i){
    IntegerVector tmp = multi_from_flat_rcpp(i);
    res(i, _) = tmp;
  }
  return(res);
}


// [[Rcpp::export]]
List index_test(IntegerVector sizes) {
  
  // Indexing struct
  GridIndex<uint_fast8_t> multi_index_struct = GridIndex<uint_fast8_t>(sizes);
  
  // Test the two main functions
  IntegerMatrix multi_indices = multi_index_struct.multi_matrix();
  const int n = multi_indices.nrow();
  IntegerVector flat_indices = IntegerVector(n);
  for (std::size_t i = 0; i < n; ++i){
    IntegerVector tmp = multi_indices(i, _);
    const int flat_idx = multi_index_struct.flat_from_multi_rcpp(tmp);
    flat_indices.at(i) = flat_idx;
  }
  
  List res = List::create(
    _["flat_idx"] = flat_indices, 
    _["multi_idx"] = multi_indices
  );
  return(res);
}

template struct GridIndex<uint_fast8_t>; 

/*** R
# x <- sample(1L:5L, size = 15, replace = TRUE)
index_test(c(3, 3))
index_test(c(3, 3, 3))

index_test(c(3, 4, 5))
*/
