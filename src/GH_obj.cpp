#include <Rcpp.h>
using namespace Rcpp;
#include "utility_rcpp.h"

// struct GH_obj {
//   const size_t n_x, n_y; 
//   NumericVector dist_diff; 
//   IntegerVector idx1, idx2; 
//   GH_obj(const NumericMatrix& dist_x, const NumericMatrix& dist_y, const IntegerMatrix& idx) 
//     : n_x(dist_x.nrow()), n_y(dist_y.nrow()){
//     const size_t len = n_x*n_x*n_y*n_y;
//     dist_diff = NumericVector(len);
//     idx1 = IntegerVector(len), idx2 = IntegerVector(len);
//     size_t cc = 0; 
//     for (size_t i = 0; i < n_x; ++i){
//       for (size_t j = 0; j < n_y; ++j){
//         for (size_t k = 0; k < n_x; ++k){
//           for (size_t l = 0; l < n_y; ++l){
//             idx1[cc] = idx(i,j), idx2[cc] = idx(k,l);
//             double dist = abs(dist_x(i,k) - dist_y(j, l));
//             dist_diff[cc] = dist; 
//             cc++;
//           }
//         }
//       }
//     }
//   }
//   double objective(const NumericVector& mu){
//     return(sum(as<NumericVector>(mu[idx1]) * as<NumericVector>(mu[idx2]) * dist_diff));
//   }
//   
//   NumericMatrix makeQ(const NumericMatrix& dist_x, const NumericMatrix& dist_y){
//     std::map< std::pair<int, int>, size_t > idx_map_x = std::map< std::pair<int, int>, size_t >(); 
//     std::map< std::pair<int, int>, size_t > idx_map_y = std::map< std::pair<int, int>, size_t >(); 
//     
//     // Convert pairwise combinations of (0:n_x-1) and (0:n_y-1) to pair index maps
//     std::vector< std::vector<int> > x_idx, y_idx; 
//     std::vector<int> x_i = std::vector<int>(n_x), y_i = std::vector<int>(n_y);
//     std::iota(x_i.begin(), x_i.end(), 0), std::iota(y_i.begin(), y_i.end(), 0);
//     x_idx.push_back(x_i), x_idx.push_back(x_i);
//     y_idx.push_back(y_i), y_idx.push_back(y_i);
//     
//     size_t cc = 0; 
//     CartesianProduct(x_idx, [&cc, &idx_map_x](std::vector<int>& el){
//       using T = std::map< std::pair<int, int>, size_t >::key_type; 
//       T p = std::make_pair(el.at(0), el.at(1));
//       idx_map_x.emplace(p, cc++);
//     });
//     cc = 0; 
//     CartesianProduct(y_idx, [&cc, &idx_map_y](std::vector<int>& el){
//       using T = std::map< std::pair<int, int>, size_t >::key_type; 
//       T p = std::make_pair(el.at(0), el.at(1));
//       idx_map_y.emplace(p, cc++);
//     });
//     
//     // Loop through 
//     size_t n_x_sq = (n_x*n_x), n_y_sq = (n_y*n_y); 
//     NumericMatrix res = no_init_matrix(n_x_sq, n_y_sq);
//     for (auto& x_kv: idx_map_x){
//       for (auto& y_kv: idx_map_y){
//         std::pair<int, int> xx = x_kv.first;
//         std::pair<int, int> yy = y_kv.first;
//         double x_dist = dist_x.at(xx.first, xx.second);
//         double y_dist = dist_y.at(yy.first, yy.second);
//         res.at(x_kv.second, y_kv.second) = std::abs(x_dist - y_dist);
//       }
//     }
//     
//     return(res);
//   }
// };
// 
// 
// RCPP_MODULE(gh_module) {
//   class_<GH_obj>("GH_obj")
//   .constructor<NumericMatrix, NumericMatrix, IntegerMatrix>()
//   .field_readonly( "dist_diff", &GH_obj::dist_diff )
//   .field_readonly( "idx1", &GH_obj::idx1 )
//   .field_readonly( "idx2", &GH_obj::idx2 )
//   .method("objective", &GH_obj::objective)
//   .method("makeQ", &GH_obj::makeQ)
//   ;
// }

/*** R
# x <- replicate(2, rnorm(10))
# y <- replicate(2, rnorm(10))
# idx <- matrix(seq(nrow(x)^2), nrow = nrow(x))
# wut <- new(GH_obj, as.matrix(dist(x)), as.matrix(dist(y)), idx)
# test_mu <- rep(1L, nrow(x)^2)
# wut$objective(test_mu)
*/
