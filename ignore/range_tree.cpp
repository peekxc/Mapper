#include <Rcpp.h>
using namespace Rcpp;

// struct node {
//   double key;
//   std::vector<double>* sa;
//   node* left, *right, *aux;
// };
// struct range_tree {
//   node root; 
//   range_tree(int k, NumericMatrix S){ 
//     const int n = S.nrow();
//     if (n == 0){ return; } // base case 
//     if (k == 1){ // base case 
//       node v = node();
//       v.sa = new ();
//       
//       NumericVector* SA = =
//       std::partial_sort_copy (S.begin(), S.end(), SA.begin(), SA.end());
//     }
//     node v = node();
//     NumericVector S_c = S.column(k);
//     
//     v.key = median(S_c);
//     NumericVector S_l = S_c[S_c <= v.key];
//     NumericVector S_r = S_c[S_c > v.key];
//     
// 
//     
//     
//   }
// };

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
