// utility_rcpp.h
// Utility functions using Rcpp. 
#ifndef UTIL_RCPP_H
#define UTIL_RCPP_H

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]   

// Choice of base types to use for indexing/storage
using uint = unsigned int;                      // unsigned int
using uint8_t = uint_fast8_t;                   // 8+ bit signed integer type
using sint8 = int_fast8_t;                      // 8+ bit signed integer type
using uidx_t = std::size_t;                     // signed index type 
using sidx_t = std::ptrdiff_t;                  // unsigned index type 
template <typename T> 
using s_ptr = std::shared_ptr<T>;               // Shared pointer
template <typename T> 
using u_ptr = std::unique_ptr<T>;               // Unique pointer


template <typename T>
using enable_int = typename std::enable_if<std::is_integral<T>::value>;

// To use alloca portably
#include <cstdlib> // alloca
#ifdef __GNUC__
/* Includes GCC, clang and Intel compilers */
# undef alloca
# define alloca(x) __builtin_alloca((x))
#elif defined(__sun) || defined(_AIX)
/* this is necessary (and sufficient) for Solaris 10 and AIX 6: */
# include <alloca.h>
#endif

// Nice simple interface from: https://stackoverflow.com/questions/17973442/index-element-from-list-in-rcpp
template <typename WHAT>
class ListOf : public List {
public:
  template <typename T>
  ListOf( const T& x) : List(x){}
  WHAT operator[](int i){ return as<WHAT>( ( (List*)this)->operator[]( i) ) ; }
  WHAT at(int i){ return as<WHAT>( ( (List*)this)->at(i) ) ; }
};

// Given a flat 0-based index 'k' and the size 'N', returns the 0-based column index 'k' would represent 
// if 'k' indexed into the lower triangular portion of an (N x N) column-major matrix.
#ifndef INDEX_TO
#define INDEX_TO(k, n) n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5) // expects 0-based, returns 0-based
#endif

// Given a flat 0-based index 'k', the size 'N', and the column index 'i', returns the 0-based row index 'k' 
// would represent if 'k' indexed into the lower triangular portion of an (N x N) column-major matrix.
#ifndef INDEX_FROM
#define INDEX_FROM(k, n, i) k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2 // expects 0-based, returns 0-based
#endif

// namespace util {

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

#include "utility_rcpp.hpp"

#endif /* UTIL_RCPP_H */

