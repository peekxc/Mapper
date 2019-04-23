// utility_rcpp.h
// Utility functions using Rcpp. 
#ifndef UTIL_RCPP_H
#define UTIL_RCPP_H

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]   

#include <memory> //shared_ptr

// Choice of base types to use for indexing/storage
using uint = unsigned int;                      // unsigned int
using uint8_t = uint_fast8_t;                   // 8+ bit signed integer type
using sint8 = int_fast8_t;                      // 8+ bit signed integer type
using uidx_t = std::size_t;                     // signed index type 
using sidx_t = std::ptrdiff_t;                  // unsigned index type 

template <typename T>
using enable_int = typename std::enable_if<std::is_integral<T>::value>;

// Nice simple interface from: https://stackoverflow.com/questions/17973442/index-element-from-list-in-rcpp
template <typename WHAT>
class ListOf : public List {
public:
  template <typename T>
  ListOf( const T& x) : List(x){}
  WHAT operator[](int i){ return as<WHAT>( ( (List*)this)->operator[]( i) ) ; }
  WHAT at(int i){ return as<WHAT>( ( (List*)this)->at(i) ) ; }
};

// namespace util {

  IntegerVector to_ivec(std::vector< std::size_t > v);
  std::vector< std::size_t > to_vec(IntegerVector v);

  sidx_t index_lower_triangular(sidx_t from, sidx_t to, const sidx_t N);

  template <typename T>
  std::vector<T> seq_ij(const T i, const T j);
  
  template<typename ForwardIterator>
  std::map<int, int> get_unique_indices(ForwardIterator first, ForwardIterator last); 
  
  template <typename T> 
  std::vector<T> merge_vectors(const std::vector< std::vector<T>* >& vec);
    
  // Rcpp-specific utilities 
  IntegerMatrix make_cartesian_product(const List& vecs);
  bool any_is_in(const IntegerVector& x, const IntegerVector& y);
  IntegerMatrix rbindlist_int(std::list<IntegerVector>& lst);
  List resize_list(const List& x, int n);
  template <typename T> IntegerVector to_ivec(std::vector<T> vec);
// }

// utility_rcpp.hpp
// Contains pure-template implementations suitable for being inlined.  



// template <typename T> struct intersect_flag {       
//   typedef int iterator;
//   typedef typename T::const_reference const_reference;
//   bool flag; 
//   intersect_flag() : flag( false ) {}
//   iterator insert( iterator, const_reference ) { flag = true; return 0; }
// };
// Sort the ranges, rely on optimized intersection method with a flag to detect and break on first hit
// using v_type = std::vector<T>;
// v_type s0, s1;
// intersect_flag<v_type> intf;
// std::set_intersection( s0.begin(), s0.end(), s1.begin(), s1.end(), std::inserter( intf, 0 ));
// if ( intf.flag ) {
//   return true; 
// }
#endif /* UTIL_RCPP_H */

