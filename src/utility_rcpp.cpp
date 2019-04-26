// utility_rcpp.cpp
// Utility functions using Rcpp. 

#include "utility_rcpp.h"

// namespace util {
  std::vector< std::size_t > to_vec(IntegerVector v){
    return std::vector< std::size_t >(v.begin(), v.end());
  }
  
  // rbindlist_int: Takes a list of integer vectors and rbind's them together.
  IntegerMatrix rbindlist_int(std::list<IntegerVector>& lst){
    std::size_t n = lst.size();
    if(n == 0) { Rcpp::stop("Invalid sized list."); }
    std::size_t d = lst.front().size();
    Rcpp::IntegerMatrix res = Rcpp::no_init(n, d);
    std::size_t i = 0;
    for (std::list<IntegerVector>::iterator it = lst.begin(); it != lst.end(); ++it, ++i) {
      if (static_cast<uidx_t>((*it).size()) != d) { Rcpp::stop("Invalid integer vector size."); }
      res(i, _) = *it;
    }
    return res;
  }
  
  // Resizes a list 
  List resize_list(const List& x, int n){
    const std::size_t lst_sz = x.size();
    List y(n);
    for(std::size_t i = 0; i < lst_sz; i++) { y[i] = x[i]; }
    return(y);
  }
  
  
  // IntegerMatrix make_cartesian_product(const List& vecs){
  //   std::vector< std::vector<int> > vv(vecs.size());
  //   std::size_t n_rows = 1;
  //   uidx_t vsize = vecs.size();
  //   for (std::size_t i = 0; i < vsize; ++i){
  //     IntegerVector v = vecs.at(i);
  //     vv.at(i) = as< std::vector<int> >(v);
  //     n_rows = n_rows * v.size();
  //   }
  //   
  //   // Copy to result
  //   std::size_t i = 0;
  //   IntegerMatrix res = IntegerMatrix(n_rows, vecs.size());
  //   CartesianProduct(vv, [&i, &res](std::vector<int> idx){
  //     res(i++, _) = as<IntegerVector>(wrap(idx));
  //   });
  //   return(res);
  // }

  // template <typename T> to_ivec< enable_int::type >();
  
// } // util namespace 


template IntegerVector to_ivec<uint8_t>(std::vector<uint8_t> vec);
// template std::vector<uint8_t> seq_ij<uint8_t>(uint8_t i, uint8_t j);
// template std::vector<sidx_t> merge_vectors<sidx_t>(const std::vector< std::vector<sidx_t>* >& vec);
  

