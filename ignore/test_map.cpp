#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
#include <random>

typedef std::vector<uint_fast8_t> index_t;

index_t make_index(const IntegerVector idx){
  index_t res = index_t(idx.size()); 
  std::copy(idx.begin(), idx.end(), res.begin());
  return(res);
}

std::vector< index_t > create_index_set(const IntegerMatrix& index_set){
  std::vector< index_t > index_key_set; 
  const int n = index_set.nrow();
  for (std::size_t i = 0; i < n; ++i){
    index_t idx = make_index(index_set.row(i));
    index_key_set.push_back(idx);
  }
  return(index_key_set);
}

std::map<index_t, std::size_t> createMap(const IntegerMatrix& index_set){
  std::vector< index_t > indexes = create_index_set(index_set);
  std::map<index_t, std::size_t> res = std::map<index_t, std::size_t>();
  std::size_t i = 0; 
  std::transform(indexes.begin(), indexes.end(), std::inserter(res, res.end()), [=, &i](index_t idx){
    return(std::map<index_t, std::size_t>::value_type(idx, ++i));
  });
  return(res);
}

// std::unordered_map<index_t, std::size_t> createUMap(const IntegerMatrix& index_set){
//   std::vector< index_t > indexes = create_index_set(index_set);
//   std::unordered_map<index_t, std::size_t> res = std::unordered_map<index_t, std::size_t>();
//   std::size_t i = 0; 
//   for (int j = 0; j < indexes.size(); ++j){
//     res.insert(std::unordered_map<index_t, std::size_t>::value_type(indexes.at(j), ++i));
//   }
//   // std::transform(indexes.begin(), indexes.end(), std::inserter(res, res.begin()), 
//   //  [=, &i](index_t idx){
//   //   return(std::unordered_map<index_t, std::size_t>::value_type(idx, ++i));
//   // });
//   return(res);
// }

// [[Rcpp::export]]
void test_map(const IntegerMatrix& index_set) {
  std::map<index_t, std::size_t> t1 = createMap(index_set);
  // std::unordered_map<index_t, std::size_t> t2 = createUMap(index_set);
  std::vector< index_t > cart_prod = create_index_set(index_set);
  
  // std::default_random_engine generator;
  // std::uniform_int_distribution<int> distribution(0, index_set.nrow() - 1);
  // 
  // const int n = 10000;
  // std::vector<std::size_t> res = std::vector<std::size_t>(n);
  // 
  // auto start = std::chrono::system_clock::now();
  // for (std::size_t i = 0; i < n; ++i){
  //   int number = distribution(generator);
  //   std::size_t testing = t1.at(cart_prod.at(number));
  //   res[i] = testing + 1;
  // }
  // auto end = std::chrono::system_clock::now();
  // auto elapsed = end - start;
  // Rcout << elapsed.count() << '\n';
  // 
  // 
  // start = std::chrono::system_clock::now();
  // for (std::size_t i = 0; i < n; ++i){
  //   int number = distribution(generator);
  //   std::size_t testing = t1.at(cart_prod.at(number));
  //   res[i] = testing + 1;
  // }
  // end = std::chrono::system_clock::now();
  // elapsed = end - start;
  // Rcout << elapsed.count() << '\n';

  // for (std::size_t s = 0; s < num_intervals.size(); ++s){
  //   num_intervals.
  //   // num_intervals
  //   // uint8_t 
  // }

}

/*** R
index_set <- as.matrix(expand.grid(1:5, 1:5, 1:5))
# test_map(index_set)
# test_map()
*/
