#include <Rcpp.h>
using namespace Rcpp;

// Allows indexing lower triangular (dist objects)
#include <cmath>
// #define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1) // expects 0-based, returns lower-triangular indices, must obey order to < from
// #define INDEX_TO(k, n) n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5) // expects 0-based, returns 0-based
// #define INDEX_FROM(k, n, i) k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2 // expects 0-based, returns 0-based

using std::size_t;

inline size_t INDEX_TF(const size_t N, const size_t to, const size_t from){ 
  return size_t((N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)); 
}

inline size_t INDEX_TO(const size_t k, const size_t n){ 
  return size_t(n - 2 - std::floor(double(std::sqrt(size_t(-8*k + 4*n*(n-1)-7))/2.0 - 0.5)));
}

inline size_t INDEX_FROM(const size_t k, const size_t n, const size_t i){ 
  return size_t(k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2);
}
  
inline double euc_dist(const NumericVector& v1, const NumericVector& v2){
  return(sqrt(sum(pow(v2 - v1, 2))));
}

NumericVector dist_subset(const NumericVector& dist, IntegerVector idx);
NumericVector combine(const NumericVector& t1, const NumericVector& t2);
NumericVector concatDist(const NumericVector& dist, const int n_new, const NumericVector& new_dists);
NumericVector dist_from_to(const NumericMatrix& X_query, const NumericMatrix& X_ref);
