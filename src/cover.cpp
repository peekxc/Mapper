#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <array>
#include <limits>
#include <map>
#include "CombinationTree.hpp"
#include "discrete.h"
// [[Rcpp::plugins(cpp11)]]

using std::size_t;
using std::begin; 
using std::end;
using std::vector;
using std::pair;
using std::array; 

// Alias to a function that accepts simplices and returns void
using set_f = std::function< void(vector< size_t >) >;

// Function ptr that accepts a lambda function and returns void

// Alias to a function that accepts simplices and returns void
using set_f = std::function< void(vector< size_t >) >;

// Function ptr that accepts a lambda function operating 
// on size_t vectors that returns void
typedef void (*setHandler)( set_f );

// Problem Motivation: R accepts types SEXP, which is only method of going from C++ -> R -> C++. 
// As a result, passing *functions* first class back to R must be done via function ptr. 
// Technical Problem: std::function's are not, and cannot, be casted explicitly to function ptrs.
// Solution: std::functions may not have function ptrs. But they may be stored in structs, 
// and structs can have member functions, which can be resolved to function ptrs. So store the 
// std::function in a struct, make a member function to invoke it w/ perfect forwarding.  
template <const size_t _UniqueId, typename _Res, typename... _ArgTypes>
struct fun_ptr_helper {
public:
    typedef std::function<_Res(_ArgTypes...)> function_type;
    static void bind(function_type&& f) { instance().fn_.swap(f); }
    static void bind(const function_type& f) { instance().fn_=f; }
    static _Res invoke(_ArgTypes... args) { return instance().fn_(args...); }
    typedef decltype(&fun_ptr_helper::invoke) pointer_type;
    static pointer_type ptr() { return &invoke; }
private:
    static fun_ptr_helper& instance() {
      static fun_ptr_helper inst_;
      return inst_;
    }
    fun_ptr_helper() {}
    function_type fn_;
};

template <const size_t _UniqueId, typename _Res, typename... _ArgTypes>
typename fun_ptr_helper<_UniqueId, _Res, _ArgTypes...>::pointer_type
get_fn_ptr(const std::function<_Res(_ArgTypes...)>& f) {
  fun_ptr_helper<_UniqueId, _Res, _ArgTypes...>::bind(f);
  return fun_ptr_helper<_UniqueId, _Res, _ArgTypes...>::ptr();
}
 

// Grid iterator. 
// ni := vector representing the number of intervals per dimension
// k := highest order simplex to consider
// order := ply of the cover
// (return): a function that iterates through the k-fold combinations, constrained to the order of the cover.
// void grid_iterator(const std::vector< size_t > ni, const size_t k, const size_t order){
//   
//   // Get the number of open sets + dimension of the space
//   const size_t n_sets = std::accumulate( ::begin(ni), ::end(ni), 1, std::multiplies< size_t >());
//   const size_t dim = ni.size(); 
//   
//   // Maximum simplex dimension is the minimum between the ply and the k
//   const size_t simplex_dim = std::min(k, order);
//   
//   // Maximum different allowed between the indices 
//   const size_t max_dim_diff = std::log(order)/std::log(2);
//   
//   // Evaluates whether a given pair of multi-indices meet the threshold requirement
//   using multi_idx_t = vector< size_t >; 
//   const auto is_potential = [max_dim_diff](const multi_idx_t m1, const multi_idx_t m2) -> bool {
//     size_t max_diff = 0; 
//     for (size_t i = 0; i < m1.size() && max_diff <= max_dim_diff; ++i){
//       size_t diff = m1[i] > m2[i] ? m1[i] - m2[i] : m2[i] - m1[i];
//       if (diff > max_diff){ max_diff = diff; };
//     }
//     return(max_diff <= max_dim_diff);
//   };
//   
//   // Evaluates whether all pairs of a *partial* combination meets the above requirement
//   using comb_t = vector< int32_t >;
//   const auto valid_partial = [&is_potential, ni](const comb_t& comb){
//     const size_t partial_k = comb.size();
//     bool is_valid = true; 
// 		apply_combinations(comb.size(), 2, [&is_potential, &is_valid, ni, partial_k](vector< size_t > idx) -> bool {
// 		  multi_idx_t m1 = as_multi_idx(idx[0], partial_k, ni);
// 		  multi_idx_t m2 = as_multi_idx(idx[1], partial_k, ni);
// 		  is_valid = (is_valid && is_potential(m1, m2));
// 		  return(is_valid);
// 		});
// 		Rcout << "is valid: " << is_valid << std::endl; 
// 		return(is_valid);
//   };
//   
//   // Predicate to prune invalid combinations
//   const auto valid_comb = [valid_partial](const comb_t& comb) -> bool {
//     if (comb.size() <= 1){ return(true); }
//     return(valid_partial(comb));
//   };
//   
//   // TODO: Return closure which compares k-combination,moves the order 
//   // computation elsewhere
//   const auto tr = pruned_ct(n_sets, k+1, valid_comb);
//   IntegerVector cc = IntegerVector(k+1);
//   Rcout << cc.size() << std::endl; 
// 	for (auto& el: tr){
// 		std::copy(std::begin(el), std::end(el), std::begin(cc));
// 		Rcout << cc << std::endl;
// 	}
//   // for (size_t k_i = 1; k_i < simplex_dim; ++k_i){
//   //   
//   // }
// 
//   
//   // // Variables to alter in the loop 
//   // vector< bool > v = vector< bool >(n_sets, false);   // permutation vector 
//   // array< size_t, 2 > idx;                             // current k-length combination
//   
//   // // Initialize
//   // std::fill(v.begin(), v.begin() + 2, true);
//   // for (size_t i = 0, cc = 0; i < n_sets; ++i) { 
//   //   if (v[i]) { idx[cc++] = i; }
//   // }
//   //     
//   
// 
//   // std::prev_permutation(v.begin(), v.end()); 
//   // array< size_t, 2 > operator*() override { 
//   //       for (size_t i = 0, cc = 0; i < n; ++i) { 
//   //         if (v[i]) { idx[cc++] = i; }
//   //       }
//   //       return idx; 
//   //     }
//   // 
//   // iterator& begin() override { return *new iterator(n_sets, 0); }
//   // iterator& end() override { 
//   //   const size_t n_combs = ((n_sets*(n_sets-1))/2)-1;
//   // 
//   // std::string bitmask(K, 1); // K leading 1's
//   // bitmask.resize(N, 0); // N-K trailing 0's
//   // do {
//   //   for (int i = 0; i < N; ++i) {
//   //     if (bitmask[i]) Rcout << " " << i;
//   //   }
//   //   Rcout << std::endl;
//   // } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
// }



// Given a logical vector, returns an integer vector (1-based) of the positions in the vector which are true
IntegerVector which_true( LogicalVector x) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; i++) { if (x[i]) y.push_back(i+1); }
  return wrap(y);
}

// [[Rcpp::export]]
IntegerMatrix valid_pairs(const IntegerMatrix& level_set_pairs){
  const int n = level_set_pairs.nrow();
  const int d = level_set_pairs.ncol();
  std::vector< int > from_index = std::vector< int >();
  std::vector< int > to_index = std::vector< int >();
  from_index.reserve(n*(d - 1)), to_index.reserve(n*(d - 1)); // reserve maximum in case
  for (int i = 0; i < n; ++i){
    IntegerMatrix::ConstRow pairs = level_set_pairs(i, _);
    const int from = pairs[0];
    for (int j = 1; j < d; ++j){
      if (!IntegerVector::is_na(pairs[j])){
        from_index.push_back(from);
        to_index.push_back(pairs[j]);
      }
    }
  }
  // Copy results
  const int n_valid = from_index.size();
  IntegerMatrix result = Rcpp::no_init_matrix(n_valid, 2);
  for (int i = 0; i < n_valid; ++i){ result(i, _) = IntegerVector::create(from_index[i], to_index[i]); }
  return(result);
}

// Must be disjoint cover 
// [[Rcpp::export]]
IntegerVector constructLevelSetIndex(const NumericMatrix& x, const NumericMatrix& bnds){
  if (x.ncol() != (bnds.ncol()/2)){ Rcpp::stop("dimension of points != dimension of bounds matrix / 2."); }
  const int n_level_sets = bnds.nrow(), d = bnds.ncol()/2;
  LogicalVector level_set_test = LogicalVector(x.nrow(), true); // which pts lie in the set; use logical vector to shorten code and use vectorized &
  IntegerVector res = IntegerVector(x.nrow(), -1);
  double eps = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < n_level_sets; ++i){
    NumericMatrix::ConstRow ls_bnds = bnds.row(i); // Update level set bounds
    std::fill(level_set_test.begin(), level_set_test.end(), true);// Reset to all true prior to doing logical range checks
    for (int d_i = 0; d_i < d; ++d_i){
      level_set_test = level_set_test & ((x.column(d_i) >= (ls_bnds[d_i] - eps)) & (x.column(d_i) <= (ls_bnds[d + d_i] + eps)));
    }
    res[level_set_test] = i+1; // Record which level set each point lies in
  }
  return(res);
}

// Given an (n x 2d) matrix of min/max bounds of any 'iso-oriented' rectangles in the plane, i.e. whose edges are parallel to the coordinate axes, 
// and a matrix of point cloud data 'x', return a list of indices of points in x which fall in the level set bounds given by 'bnds'.
// [[Rcpp::export]]
List constructIsoAlignedLevelSets(const NumericMatrix& x, const NumericMatrix& bnds, bool save_bounds=true){
  if (x.ncol() != (bnds.ncol()/2)){ Rcpp::stop("dimension of points != dimension of bounds matrix / 2."); }
  const int n_level_sets = bnds.nrow(), d = bnds.ncol()/2;
  List level_sets = List(n_level_sets);
  LogicalVector level_set_test = LogicalVector(x.nrow(), true); // which pts lie in the set; use logical vector to shorten code and use vectorized &
  for (int i = 0; i < n_level_sets; ++i){
    NumericMatrix::ConstRow ls_bnds = bnds.row(i); // Update level set bounds
    std::fill(level_set_test.begin(), level_set_test.end(), true);// Reset to all true prior to doing logical range checks
    for (int d_i = 0; d_i < d; ++d_i){
      level_set_test = level_set_test & ((x.column(d_i) >= ls_bnds[d_i]) & (x.column(d_i) <= ls_bnds[d + d_i]));
    }

    // Save the level set
    IntegerVector ls = which_true(level_set_test);
    if (save_bounds){ ls.attr("bounds") = ls_bnds; }
    level_sets[i] = ls;
  }
  return(level_sets);
}

// [[Rcpp::export]]
List constructFixedLevelSets(const NumericMatrix& filter_values,
                             const IntegerMatrix& index_set,
                             const NumericVector& overlap,
                             const IntegerVector& number_intervals,
                             const NumericMatrix& filter_range,
                             const NumericVector& filter_len) {
  // Base variables
  const int n = index_set.nrow();
  const int d = index_set.ncol();
  NumericMatrix ls_bnds = NumericMatrix(2, d);
  List level_sets = List(n);
  LogicalVector level_set_test = LogicalVector(filter_values.nrow(), true);
  const NumericVector base_interval_length = filter_len/as<NumericVector>(number_intervals);
  const NumericVector filter_min = filter_range(0, _);
  const NumericVector interval_length = base_interval_length + (base_interval_length * overlap)/(1.0 - overlap);
  // Rcout << interval_length << std::endl;
  const NumericVector eps = interval_length/2.0;
  for (int i = 0; i < n; ++i){
    IntegerVector ls_idx  = index_set(i, _) - 1;

    // Compute the level sets bounds
    NumericVector level_set_centroid = filter_min + (as<NumericVector>(ls_idx) * base_interval_length) + base_interval_length/2.0;
    ls_bnds(0, _) = level_set_centroid - eps;
    ls_bnds(1, _) = level_set_centroid + eps;

    // Reset to all true prior to doing logical range checks
    std::fill(level_set_test.begin(), level_set_test.end(), true);
    for (int d_i = 0; d_i < d; ++d_i){
      level_set_test = level_set_test & ((filter_values.column(d_i) >= ls_bnds(0, d_i)) & (filter_values.column(d_i) <= ls_bnds(1, d_i)));
    }

    // Don't explicitly need to save the bounds, but they may useful later
    level_sets[i] = List::create(_["points_in_level_set"] = which_true(level_set_test), _["bounds"] = clone(ls_bnds));
  }
  return(level_sets);
}


// [[Rcpp::export]]
List constructRestrainedLevelSets(const NumericMatrix& filter_values, const IntegerMatrix& index_set, const NumericVector& interval_length, const NumericVector& step_size, const NumericVector& filter_min) {
  const int n = index_set.nrow();
  const int d = index_set.ncol();

  NumericMatrix ls_bnds = NumericMatrix(2, d);
  List level_sets = List(n);
  LogicalVector level_set_test = LogicalVector(filter_values.nrow(), true);
  for (int i = 0; i < n; ++i){
    IntegerVector ls_idx  = index_set(i, _) - 1;
    NumericVector level_set_min = filter_min + (as<NumericVector>(ls_idx) * step_size);
    ls_bnds(0, _) = level_set_min;
    ls_bnds(1, _) = level_set_min + interval_length;

    // Reset to all true prior to doing logical range checks
    std::fill(level_set_test.begin(), level_set_test.end(), true);
    for (int d_i = 0; d_i < d; ++d_i){
      level_set_test = level_set_test & ((filter_values.column(d_i) >= ls_bnds(0, d_i)) & (filter_values.column(d_i) <= ls_bnds(1, d_i)));
    }

    // Don't explicitly need to save the bounds, but they may useful later
    level_sets[i] = List::create(_["points_in_level_set"] = which_true(level_set_test), _["bounds"] = clone(ls_bnds));
  }
  return(level_sets);
}

// Creates a mapping in the form of the tuple (i, j) representing whether level set i in Mapper #1
// is a subset of (or equal to) level set j of Mapper #2. Returns an integer matrix, where each
// row delimits this mapping.
// Expects as input:
// ls1 := list of 2xd numeric matrixes representing the bounds of each level set for the first cover.
// ls2 := list of 2xd numeric matrixes representing the bounds of each level set for the second cover.
// d := the dimensionality of the covers
// [[Rcpp::export]]
NumericMatrix createCoverMap(const List& ls1, const List& ls2, const int d){
  const int n_ls_1 = ls1.length();
  const int n_ls_2 = ls2.length();
  const int min_maps = std::max(n_ls_1, n_ls_2);

  // Create the vector to represent the lsfi maps
  std::vector<int> from_lsfi = std::vector<int>();
  std::vector<int> to_lsfi = std::vector<int>();
  from_lsfi.reserve(min_maps); // Every subset of the first cover will map to at least one subset of the second cover
  to_lsfi.reserve(min_maps); // Every subset of the second cover will get mapped to by at least one subset of the first cover
  for (int i = 0; i < n_ls_1; ++i){
    for (int j = 0; j < n_ls_2; ++j){
      const NumericMatrix& ls_bnds_i = ls1[i];
      const NumericMatrix& ls_bnds_j = ls2[j];
      bool all_intersect = true;
      for (int d_i = 0; d_i < d && all_intersect; ++d_i){
        all_intersect = all_intersect & bool(ls_bnds_i(0, d_i) <= ls_bnds_j(1, d_i) && ls_bnds_i(1, d_i) >= ls_bnds_j(0, d_i));
      }
      if (all_intersect){ from_lsfi.push_back(i), to_lsfi.push_back(j); }
    }
  }
  const int n = from_lsfi.size();
  NumericMatrix res = NumericMatrix(n, 2);
  for (int i = 0; i < n; ++i){
    res(i, _) = NumericVector::create(from_lsfi[i], to_lsfi[i]);
  }
  return(res);
}

// [[Rcpp::export]]
List edgelist_to_adjacencylist(const IntegerMatrix& el){
  const int n = el.nrow();
  std::map<int, std::vector<int> > vertex_map = std::map< int, std::vector<int> >();
  for (int i = 0; i < n; ++i){
    const int from = el(i, 0), to = el(i, 1);
    std::map<int, std::vector<int> >::iterator it = vertex_map.lower_bound(from);
    if (it != vertex_map.end() && it->first == from){
      it->second.push_back(to); // add to index to adjacency list
    } else {
      std::vector<int> to_ids = std::vector<int>();
      to_ids.push_back(to);
      vertex_map.emplace_hint(it, std::make_pair(from, to_ids));
    }
  }
  return(wrap(vertex_map));
}

// Equivalent to R's which, but in C++
IntegerVector which_cpp(const IntegerVector& x, int value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] == value) y.push_back(i); }
  return wrap(y);
}

// Finds the first index (1-based) of x which equals the integer value given.
// For some reason, I couldn't find an efficient version of this in R, and Position/Find is slow.
// [[Rcpp::export]]
int findFirstEqual(const IntegerVector& x, int value){
  IntegerVector::const_iterator start = x.begin();
  return(std::distance(start, std::find(start, x.end(), value)) + 1);
}



// Creates an adjacency list connecting two lists of nodes...
// IntegerMatrix nodeMap(const IntegerVector& node_lsfi1, const IntegerVector& node_lsfi2, const IntegerMatrix& cover_map){
//   const int n = cover_map.nrow();
//   std::vector<int> from_node_comb = std::vector<int>();
//   std::vector<int> to_node_comb = std::vector<int>();
//   for (int i = 0; i < n; ++i){
//     const int from_idx = cover_map(i, 0), to_idx = cover_map(i, 1);
//     IntegerVector from_nodes = which_cpp(node_lsfi1, from_idx);
//     IntegerVector to_nodes = which_cpp(node_lsfi2, to_idx);
//     for (IntegerVector::const_iterator f = from_nodes.begin(); f != from_nodes.end(); ++f){
//       for (IntegerVector::const_iterator t = to_nodes.begin(); t != to_nodes.end(); ++t){
//         from_node_comb.push_back(*f);
//         to_node_comb.push_back(*t);
//       }
//     }
//   }
//   const int n2 = from_node_comb.size();
//   IntegerMatrix res = IntegerMatrix(n2, 2);
//   for (int i = 0; i < n2; ++i){
//     res(i, _) = IntegerVector::create(from_node_comb[i], to_node_comb[i]);
//   }
//   return(res);
// }

// Computes the absolute distances from a given point to the closest endpoint of each level set. This distance should 
// represent 1/2 the smallest interval length the target level set would have to be (via expansion) to intersect the given point.
// [[Rcpp::export]]
List dist_to_boxes(const IntegerVector& positions, const double interval_length, const int num_intervals, const NumericVector& dist_to_lower, const NumericVector& dist_to_upper) {
  
  // Sequence from 1 - < number of intervals >  
  std::vector<int> all_positions = std::vector<int>(num_intervals);
  std::iota(std::begin(all_positions), std::end(all_positions), 1);
  
  // Iterators 
  IntegerVector::const_iterator pos_it = positions.begin();
  NumericVector::const_iterator dtl_it = dist_to_lower.begin(), dtu_it = dist_to_upper.begin();
  //IntegerVector::iterator k; 
  
  // Variables needed 
  const int n = positions.size();
  std::array<int, 1> current_position = { {0} };
  
  // To fill each iteration
  IntegerVector target_positions = no_init(num_intervals - 1);
  NumericVector target_distances = no_init(num_intervals - 1);
  
  // Outputs 
  IntegerMatrix res_pos = no_init_matrix(n, num_intervals - 1);
  NumericMatrix res_dist = no_init_matrix(n, num_intervals - 1);
  
  double dtl = 0.0, dtu = 0.0;
  for (int i = 0, pos = 0; i < n; ++i, ++pos_it, ++dtl_it, ++dtu_it){
    pos = *pos_it, dtl = *dtl_it, dtu = *dtu_it;
    current_position[0] = pos;
    
    // Only compute distances to level sets not intersecting the current point
    std::set_difference(all_positions.begin(), all_positions.end(), 
                        current_position.begin(), current_position.end(), 
                        target_positions.begin());
    
    // Distance calculation
    std::transform(target_positions.begin(), target_positions.end(), target_distances.begin(), 
       [interval_length, pos, dtl, dtu](int target_position){
         double offset = std::abs((target_position - pos - 1) * interval_length);
         return (target_position < pos ? dtl + offset : dtu + offset);
    });
    
    // Store results
    res_dist.row(i) = clone(target_distances);
    res_pos.row(i) = clone(target_positions);
  }
  return(List::create(_["target_pos"] = res_pos, _["target_dist"] = res_dist));
}



/*** R
# load("test.rdata")
# fv <- as.matrix(test$filter_values)
# is <- as.matrix(test$index_set)
# il <- as.numeric(test$interval_length)
# ss <- as.numeric(test$step_size)
# fmin <- as.numeric(test$filter_min)
# constructLevelSets(filter_values = fv, index_set = is, interval_length = il, step_size = ss, filter_min = fmin)
*/
