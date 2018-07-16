// #include <Rcpp.h>
// using namespace Rcpp;
// 
// const int N = 1e5;  // limit for array size
// #include <math.h>
// 
// 
// 
// // int query(int l, int r) {  // sum on interval [l, r)
// //   int res = 0;
// //   for (l += n, r += n; l < r; l >>= 1, r >>= 1) {
// //     if (l&1) res += t[l++];
// //     if (r&1) res += t[--r];
// //   }
// //   return res;
// // }
// 
// // A segment tree stores a balanced binary tree of n endpoints and (n - 1) intervals.
// struct SegmentTree {
//   int* tr; // the segment tree
//   const int N, n; // number of endpoints and intervals, respectively
//   const int h; // tree height
//   SegmentTree(const NumericVector& endpoints, const IntegerVector& idx) : N(endpoints.size()), n(N / 2), h(ceil(log2(n+1))) {
//     tr = new int[N];
//     // std::copy(idx.begin(), idx.end(), tr[N]); // insert the endpoints as leaves
//     build(); // build the internal nodes
//   }
// 
//   // Given a query point 'q', populates a vector 'intervals' with the indexes of the intervals 'q' intersects with
//   // void queryPoint(double q, std::vector<int>& intervals) {
//   //   for (int i = 0; i < )
//   //   tr[]
//   //   endpoints[]
//   //   return res;
//   // }
// 
//   // To query an (index-converted) interval
//   int queryInterval(int l, int r){
//     int res = 0;
//     for (l += n, r += n; l < r; l >>= 1, r >>= 1) {
//       if (l&1) res += tr[l++];
//       if (r&1) res += tr[--r];
//     }
//     return res;
//   }
// 
// 
//   // Given the leaves are defined, builds the internal nodes of the tree bottom-up
//   void build() {
//     for (int i = N - 1; i > 0; --i) { tr[i] = tr[i<<1] + tr[i<<1|1]; } // equivalent to: tr[i] = tr[2*i] + tr[2*i+1]
//   }
// };
// 
// // [[Rcpp::export]]
// NumericVector something(NumericVector x) {
//   return x * 2;
// }
// 
// 
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically
// // run after the compilation.
// //
// 
// /*** R
// */
