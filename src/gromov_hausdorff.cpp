#include <Rcpp.h>
using namespace Rcpp;

#include "utility_rcpp.h"
#include <cmath>

// std::pair<NumericMatrix, NumericMatrix> make_constraints(const size_t N){
//   NumericMatrix A = NumericMatrix(2*N, pow(N, 2), 0);
//   NumericMatrix B = NumericMatrix(2*N, 1, 1);
//   for (size_t i = 0; i < N; ++i){
//     std::fill()
//   }
// }

// function [A,b] = constraints(N)
//   A=zeros(2*N,N^2);
// b=ones(2*N,1);
// for i=1:N
//   A(i,(i-1)*N+1:i*N)=ones(1,N);
// for j=1:N
//   A(i+N, (j-1)*N+i)=1;
// end
//   end
//   end

// Traverses a choose(n, 2) distance vector as-if it was applying a row operation to a full (n x n) matrix
template<typename Func>
void doColumn(const NumericVector& dist_x, const size_t n, const size_t j, Func f) {
  size_t i = 0; 
  std::vector<size_t> idx = std::vector<size_t>();
  idx.reserve(n);
  std::generate_n(idx.begin(), n, [&i, &j, &n](){ return(index_lower_triangular(i++, j, n)); });
  std::for_each(idx.begin(), idx.end(), [&dist_x, &f](const size_t ii){
    f(dist_x[ii]);
  });
}

// Traverses a choose(n, 2) distance vector as-if it was applying a row operation to a full (n x n) matrix
template<typename Func>
void doRow(const NumericVector& dist_x, const size_t n, const size_t i, Func f) {
  size_t j = 0; 
  std::vector<size_t> idx = std::vector<size_t>();
  idx.reserve(n);
  std::generate_n(idx.begin(), n, [&i, &j, &n](){ return(index_lower_triangular(i, j++, n)); });
  std::for_each(idx.begin(), idx.end(), [&dist_x, &f](const size_t ii){
    f(dist_x[ii]);
  });
}

// Apply a lambda to a distance vector representing column-major lower-triangular distances
// template<typename Func>
// void apply_lt(const NumericVector& dist_x, Func f, const bool diagonal=false) {
//   const size_t n = dist_x.attr("Size");
//   NumericVector::const_iterator c_dist = dist_x.begin();
//   for (int i = 0; i < n; ++i){
//     for (int j = 0; j < i + int(diagonal); ++j){
//       f(i, j, )
//     }
//   }
//   std::for_each(dist_x.begin(), dist_x.end(), [](){
//     
//   });
// }

NumericVector wrap_dist(NumericVector v, const size_t Size, std::string method = ""){
  v.attr("class") = "dist";
  v.attr("Size") = Size;
  v.attr("Diag") = false;
  v.attr("Upper") = false;
  if (method != ""){ v.attr("method") = method; } 
  return(v);
}

// Assume idx is 0-based
// [[Rcpp::export]]
IntegerMatrix gh_make_A(const IntegerMatrix& idx){
  const std::size_t n = idx.nrow(), m = idx.ncol(); 
  const std::size_t n_vars = n * m;
  
  // Fill in constraints row-by-row
  IntegerVector tmp(n_vars, 0);
  IntegerMatrix constraints = Rcpp::no_init_matrix(n + m, n_vars);
  
  // Set lambda indicators 
  for (std::size_t i = 0; i < n; ++i){
    const IntegerVector ii = idx.row(i);
    tmp[ii] = 1;
    constraints.row(i) = tmp;
    std::fill(tmp.begin(), tmp.end(), 0);
  }

  // Set nu indicators 
  for (std::size_t j = 0; j < m; ++j){
    const IntegerVector jj = idx.column(j);
    tmp[jj] = 1;
    constraints.row(n + j) = tmp;
    std::fill(tmp.begin(), tmp.end(), 0);
  }
  return(constraints);
}

// [[Rcpp::export]]
NumericMatrix gh_make_Q(const NumericMatrix& x_dist, const NumericMatrix& y_dist){
  using idx_t = std::vector< std::size_t >;
  std::size_t n_x = x_dist.nrow(), n_y = y_dist.nrow();
  const std::size_t N = n_x*n_y; 
  
  // Generate X and Y indices
  idx_t x_idx(n_x), y_idx(n_y);
  std::iota(x_idx.begin(), x_idx.end(), 0);
  std::iota(y_idx.begin(), y_idx.end(), 0);
  
  // Declare order to generate distance permutations
  idx_t order = { 0, 1, 2, 3 };
  std::vector< idx_t > to_permute = { x_idx, y_idx, x_idx, y_idx };
  
  // Use ordered cartesian product to fill Q row-wise
  // std::size_t ii = 0, jj = 0; 
  NumericMatrix Q = Rcpp::no_init_matrix(N, N);
  NumericMatrix::iterator q_it = Q.begin();
  CartesianProductOrdered(to_permute, [&q_it, &x_dist, &y_dist](const idx_t el){
    std::size_t i = el.at(0), j = el.at(1), k = el.at(2), l = el.at(3);
    (*q_it) = abs(x_dist(i, k) - y_dist(j, l));
    ++q_it;
  }, order);
  
  return(Q);
}


// This is more-or-less capped at around max(n) ~= 150^2
// Assumes X and Y are dist objects with 0 along the diagonals
// [[Rcpp::export]]
NumericMatrix all_correspondences(const NumericVector& X, const NumericVector& Y) {
  if (as< size_t >(X.attr("Size")) != as< size_t >(Y.attr("Size"))){ stop("X and Y must be the same size."); }
  const size_t n = X.attr("Size"); // number of points in X, Y
  const size_t n_sq = pow(n, 2); // size of distance matrix of X, Y
  size_t i, j, k, l;
  NumericMatrix gamma = no_init_matrix(n_sq, n_sq);
  for (i = 0; i < n; ++i){
    for (j = 0; j < n; ++j){
      for (k = 0; k < n; ++k){
        for(l = 0; l < n; ++l){
          size_t a = (i*n) + j, b = (k*n) + l;
          // Rprintf("{a:%d,b:%d}(i:%ld,j:%ld,k:%ld,l:%ld)\n", a, b, i, j, k, l);
          const double dist1 = i == k ? 0 : X.at(index_lower_triangular(i, k, n));
          const double dist2 = j == l ? 0 : Y.at(index_lower_triangular(j, l, n));
          gamma.at(a, b) = std::abs(dist1 - dist2);
        }
      }
    }
  }
  return(gamma);
  // Fill in diagonal
  // for (i = 0; i < n_sq; ++i){
  //   size_t ii = index_lower_triangular(i, i, n_sq);
  //   
  //   gamma(ii, ii) = std::abs(X.at(i) - Y.at(jj));
  // }
  // return(gamma);
  
  // **Vector version**
  // NumericVector gamma = no_init_vector(n_res);
  // const size_t n = X.attr("Size");
  // const size_t n_dist = X.size(); 
  // const size_t n_res = n_dist*(n_dist-1)/2;
  // size_t i = 0; 
  // std::generate_n(gamma.begin(), n_dist, [&i, &n, &n_dist, &X, &Y](){
  //   const size_t idx = i;
  //   const size_t to = INDEX_TO(idx, n_dist);
  //   const size_t from = INDEX_FROM(idx, n_dist, to);
  //   const size_t idx2 = index_lower_triangular(from, to, n);
  //   Rprintf("idx: %d, idx2: %d\n", i, idx2);
  //   i++;
  //   return(std::abs(X.at(idx2) - Y.at(idx2)));
  // });
  // return(wrap_dist(gamma, n_dist));
}



/*** R
x <- replicate(2, rnorm(20))
y <- replicate(2, rnorm(20))
wut <- Mapper:::gromov_hausdorff(dist(x), dist(y))
*/
