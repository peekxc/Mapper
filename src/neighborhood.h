#include <Rcpp.h>
using namespace Rcpp;

#include "discrete.h"
#include "combinations.h"
#include "delegate.hpp"

using std::vector;
using std::size_t;
using set_f = std::function< void(vector< size_t >) >;
using del_f = delegate<void (set_f)>;

// Exports a generic 'neighborhood' for any cover, which is an iterator through 
// k-combinations of a combination tree
// n := number of sets
// k := k-length combinations to check
// [[Rcpp::export]]
SEXP generic_neighborhood(const size_t n, const size_t k){
  using idx_t = size_t;
  using comb_t = vector< idx_t >; 
  const del_f f = [n,k](set_f nerve_f) -> void {
    for_each_combination_idx(idx_t(n), idx_t(k), nerve_f);
  };
  SEXP rng_sexp = XPtr< del_f >(new del_f(f));
  return(rng_sexp);
}

// Converts a flat index to a multi index
void as_multi_idx(size_t idx, const size_t d, const vector< size_t >& ni, size_t* result){
  size_t rem;
  for (size_t d_i = 0; d_i < d; ++d_i) {
    rem = idx % ni[d_i];
    result[d_i] = rem;
    idx -= rem;
    idx /= ni[d_i];
  }
}
  
// Evaluates whether a given set of linear indices 'idx' is a valid combination 
// to test for nonempty intersection, depending on the distance between their indices.
// Rprintf(" === (%d) diff: %d, threshold: %d\n", d_i, diff, diff_thresh[d_i]);
// Rprintf("idx0: %d, idx1: %d, ni_size: %d\n", idx.at(ii.at(0)), idx.at(ii.at(1)), ni.size());
bool within_neighborhood(const vector< size_t > idx, const vector< size_t > ni, const vector< size_t > diff_thresh) {
  using It = typename vector< size_t >::iterator;
  bool is_valid = true; 
  
  // Use stack to store the multi-indices
  const size_t dim = static_cast< size_t >(ni.size());
  size_t* m1 = static_cast< size_t*>(alloca(sizeof(size_t)*dim));
  size_t* m2 = static_cast< size_t*>(alloca(sizeof(size_t)*dim));
  
  // Iterate through all combinations, invalidating the current idx if any pass the difference threshold 
  const size_t k = static_cast< size_t >(idx.size()); 
  for_each_combination(k, size_t(2), [&is_valid, &m1, &m2, diff_thresh, idx, ni, dim](It a, It b) -> bool {
	  as_multi_idx(idx[*a], dim, ni, m1);
	  as_multi_idx(idx[*(a+1)], dim, ni, m2);
	  for (size_t d_i = 0; d_i < dim; ++d_i){
	    size_t diff = m1[d_i] > m2[d_i] ? m1[d_i] - m2[d_i] : m2[d_i] - m1[d_i];
	    is_valid = is_valid && (diff <= diff_thresh[d_i]); 
	  }
	  return !is_valid;
	});
  
	return(is_valid);
}

// Exports the 'neighborhood' of the fixed interval cover, which is an function pointer that, when invoked,
// iterates through through k-fold combinations in a pruned combination tree. 
// ni           := number of intervals per dimension
// overlap      := overlap percentage 
// k            := dimension of simplices to consider
// filter_len   := range of each dimension
// [[Rcpp::export]]
SEXP fixed_interval_neighborhood(std::vector< size_t > ni, std::vector< double > overlap, const size_t k, std::vector< double > filter_len) {
  using idx_t = size_t;
  
  // Get the number of open sets + dimension of the space
  const size_t n_sets = std::accumulate( ::begin(ni), ::end(ni), 1, std::multiplies< size_t >());
  const size_t dim = ni.size(); 
  
  // Based on the covers ply, determine the maximum pairwise difference 
  // valid combinations of indices can have
  vector< size_t > diff_thresh = vector< size_t >(dim, 0);
  for (idx_t d_i = 0; d_i < dim; ++d_i){
    const double base_len = filter_len[d_i]/ni[d_i];
    const double prop_overlap = overlap[d_i]/100.0;
    const double c_len = base_len + (base_len * prop_overlap)/(1.0 - prop_overlap);
    for (size_t i = 0; i < ni[d_i]; ++i){
      double crit = base_len + ((base_len/2.0)*i)*2;
      if (crit < c_len){
        diff_thresh.at(d_i) = i+1;  
      }
    }
  }

  // Export a delegate to contain the function pointer
	del_f f = [n_sets, k, ni, diff_thresh](set_f nerve_f) -> void {
    for_each_combination_idx(size_t(n_sets), size_t(k+1), [&ni, &diff_thresh, &nerve_f](vector< size_t > comb){
      if (within_neighborhood(comb, ni, diff_thresh)){
        nerve_f(comb);
      }
    });
  };
  SEXP rng_sexp = XPtr< del_f >(new del_f(f), true);  
    
  return(rng_sexp);
}

void apply_neighborhood(SEXP neighborhood_f, Function f){

  using del_f = delegate<void (set_f)>;
  XPtr< del_f > rng_fun(neighborhood_f); // Unwrap
  del_f rng_gen = *rng_fun;
  // size_t res = 0.0;
  const set_f g = [&f](vector< size_t > pids) -> void {
    //res += std::accumulate(begin(pids), end(pids), 0);
    f(pids);
    return;
  };
  rng_gen(g);
}


/*** R
#ni <- c(25,25)
ni <- c(5,5,5)
g <- c(25,25,25)
k <- 1
len <- c(10, 7, 9)

## Testing multiple solution to iterating multiple combinations 
# f_red <- fixed_interval_neighborhood(ni,g,k+1,len)
# f_gen <- generic_neighborhood(prod(ni), k+2)
# g1 <- function() { apply(combn(prod(ni),k+2), 2, sum) }
# g2 <- function() { colSums(combn(prod(ni),k+2)) }
# 
# microbenchmark::microbenchmark({ g1() }, times = 10L)
# microbenchmark::microbenchmark({ g2() }, times = 10L)
# microbenchmark::microbenchmark({ apply_neighborhood(f_gen, function(x){}) }, times = 10L)
# microbenchmark::microbenchmark({ apply_neighborhood(f_red, function(x){}) }, times = 10L)






# res <- vector("list", length = choose(25,2))
# i <- 1L
# apply_neighborhood(f, function(x){
#   print(as.integer(x))
#   res[[i]] <<- x
#   i <<- i + 1L
# })
*/

  // XPtr< setHandler > rng_fun(neighborhood_f); // Unwrap
  // setHandler rng_gen = *rng_fun;
  // const set_f g = [&f](vector< size_t > pids) -> void {
  //   // std::size_t res = std::accumulate(begin(pids), end(pids), 0);
  //   for (auto el: pids){ Rcout << el << ","; };
  //   Rcout << std::endl;
  //   return;
  // };
  // (*rng_gen)(g);
  // const auto ni = vector< size_t >{ 5, 5 };
  // const auto diff_thresh = vector< size_t >{ 1, 1 };
  // for_each_combination_idx(size_t(25), size_t(2), [&ni, &diff_thresh, &f](vector< size_t > comb){
  //   if (within_neighborhood(comb, ni, diff_thresh)){
  //     // Rcout << "is valid" << std::endl;
  //     f(comb); 
  //   }
  // });
  
  // for_each_combination_idx(static_cast< size_t >(9), static_cast< size_t >(2), [](vector< size_t > comb){
  //   
  //   Rcout << "combinations: ";
  //   for (auto el: comb){ Rcout << el << ","; }
  //   Rcout << std::endl;
  //   // Rcout << "end" << std::endl;
  //   // if (valid_comb(comb)){
  //   //   Rcout << "is valid" << std::endl;
  //   //   nerve_f(comb); 
  //   // }
  // });


// Function ptr that accepts a lambda function operating 
// on size_t vectors that return void
// typedef void (*setHandler)( set_f );

// Problem Motivation: R accepts types SEXP, which is only method of going from C++ -> R -> C++ anonymously. 
// As a result, passing *functions* first class back to R must be done via function ptr. 
// Technical Problem: std::function's are not, and cannot, be function ptrs.
// Solution: std::functions may not have function ptrs. But they may be stored in structs, 
// and structs can have member function ptrs. So store the std::function in a struct,
// make a member function to invoke it w/ perfect forwarding.  
// template <const size_t _UniqueId, typename _Res, typename... _ArgTypes>
// struct fun_ptr_helper {
// public:
//     typedef std::function<_Res(_ArgTypes...)> function_type;
//     static void bind(function_type&& f) { instance().fn_.swap(f); }
//     static void bind(const function_type& f) { instance().fn_=f; }
//     static _Res invoke(_ArgTypes... args) { return instance().fn_(args...); }
//     typedef decltype(&fun_ptr_helper::invoke) pointer_type;
//     static pointer_type ptr() { return &invoke; }
// private:
//     static fun_ptr_helper& instance() {
//       static fun_ptr_helper inst_;
//       return inst_;
//     }
//     fun_ptr_helper() {}
//     function_type fn_;
// };
// 
// template <const size_t _UniqueId, typename _Res, typename... _ArgTypes>
// typename fun_ptr_helper<_UniqueId, _Res, _ArgTypes...>::pointer_type
// get_fn_ptr(const std::function<_Res(_ArgTypes...)>& f) {
//   fun_ptr_helper<_UniqueId, _Res, _ArgTypes...>::bind(f);
//   return fun_ptr_helper<_UniqueId, _Res, _ArgTypes...>::ptr();
// }
  // Evaluates whether the two given multi-indices have distance less than the per-dimension threshold 
  // using multi_idx_t = vector< size_t >; 
  // const auto is_potential = [diff_thresh](const multi_idx_t m1, const multi_idx_t m2) -> bool {
  //   const size_t dim = m1.size(); 
  //   for (size_t d_i = 0; d_i < dim; ++d_i){
  //     size_t diff = m1[d_i] > m2[d_i] ? m1[d_i] - m2[d_i] : m2[d_i] - m1[d_i];
  //     Rprintf(" === (%d) diff: %d, threshold: %d\n", d_i, diff, diff_thresh[d_i]);
  //     if (diff > diff_thresh[d_i]){ 
  //       return false; 
  //     };
  //   }
  //   return(true);
  // };
  
  // Evaluates whether all pairs of a *partial* combination meets the above requirement
//   using comb_t = vector< size_t >;
//   const auto valid_partial = [is_potential, ni](const comb_t& comb){
//     bool is_valid = true; 
//     for_each_combination_idx(size_t(comb.size()), size_t(2), [&is_potential, &is_valid, comb, ni](vector< size_t > idx) -> void {
// 		  Rprintf("idx0: %d, idx1: %d, ni_size: %d\n", idx[0], idx[1], ni.size());
// 		  multi_idx_t m1 = as_multi_idx(comb[idx[0]], ni.size(), ni);
// 		  multi_idx_t m2 = as_multi_idx(comb[idx[1]], ni.size(), ni);
// 		  is_valid = (is_valid && is_potential(m1, m2));
// 		  return;
// 		});
// 		return(is_valid);
//   };
  
  // Predicate to prune invalid combinations
    // const auto valid_comb = [valid_partial](const comb_t& comb) -> bool {
    //   if (comb.size() <= 1){ return(true); }
    //   return(valid_partial(comb));
    // };
  
  // Copy the combination tree into a function
  // Rcout << "k=" << k << std::endl;
  // const std::function< void(set_f) > f = [n_sets, k, ni, diff_thresh](set_f nerve_f) -> void {
  //   // Rcout << "Checking: " << n_sets << " choose " << k+1 << " combinations" << std::endl;
  //   for_each_combination_idx(n_sets, k+1, [&nerve_f, &ni, &diff_thresh](vector< size_t > comb){
  //     if (within_neighborhood(comb, ni, diff_thresh)){
  //       // Rcout << "is valid" << std::endl;
  //       nerve_f(comb);
  //     }
  //   });
  // };
  
  // SEXP rng_sexp = XPtr< setHandler >(new setHandler(get_fn_ptr<0>(f)));
    // for_each_combination_idx(n_sets, k+1, [&nerve_f, &ni, &diff_thresh](vector< size_t > comb){
    //   if (within_neighborhood(comb, ni, diff_thresh)){
    //     nerve_f(comb);
    //   }
    // });
//     Rprintf("n_sets: %d, k: %d\n", n_sets, k);
// 	  for (auto el: ni){ Rcout << el << ", "; };
// 	  Rcout << std::endl;
// 	  for (auto el: diff_thresh){ Rcout << el << ", "; };
// 	  Rcout << std::endl;
// 	  
    // const auto ni = vector< size_t >{ 5, 5 };
    // const auto diff_thresh = vector< size_t >{ 1, 1 };
