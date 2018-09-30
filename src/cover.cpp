#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <array>
#include <limits>

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

// Returns a list of
List getUpdateBlocks(const NumericMatrix& overlap, const IntegerMatrix& target_lsfi, const IntegerMatrix& block_idx,
                     const IntegerVector& block_ids){
  // const int d = g.length();
  const int n = block_idx.nrow();
  const int n_blocks = block_ids.size();
  // for (int di = 0; di < d; ++di){
    // const List params_di = g[di];
    // const NumericMatrix overlap = params_di["overlap"];
    // const IntegerMatrix target_lsfi = params_di["target_lsfi"];
  List output = List(n_blocks);
  // std::vector< std::list< IntegerVector > > level_sets_to_update =  std::vector< std::list< IntegerVector > >(n_blocks);
  // for (int i = 0; i < n; ++i){
  //   IntegerMatrix::ConstRow i_blocks = block_idx(i, _);
  //   IntegerMatrix::ConstRow i_lsfis = target_lsfi(i, _);
  //   IntegerMatrix::ConstRow::const_iterator i_block = i_blocks.begin();
  //   IntegerMatrix::ConstRow::const_iterator i_lsfi = i_lsfis.begin();
  //   for (; i_block != i_blocks.end(); ++i_block, ++i_lsfi){
  //     level_sets_to_update.at(*i_block)[*i_lsfi].push_back();
  //   }
  // }
  // }
  return(output);
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
      Rcout << "testing level set: min = " << ls_bnds[d_i] <<  ", max = " << ls_bnds[d + d_i] << std::endl;
      level_set_test = level_set_test & ((x.column(d_i) >= (ls_bnds[d_i] - eps)) & (x.column(d_i) <= (ls_bnds[d + d_i] + eps)));
    }
    Rcout << level_set_test << std::endl; 
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

#include <map>
// [[Rcpp::plugins(cpp11)]]
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

IntegerVector which_cpp(const IntegerVector& x, int value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] == value) y.push_back(i); }
  return wrap(y);
}

// Finds the first index (1-based) of x which equals the integer value given.
// For some reason, I couldn't find this in R, and Position/Find is slow.
// [[Rcpp::export]]
int findFirstEqual(const IntegerVector& x, int value){
  IntegerVector::const_iterator start = x.begin();
  return(std::distance(start, std::find(start, x.end(), value)) + 1);
}

// Creates an adjacency list representing the simplicial mapping
// [[Rcpp::export]]
IntegerMatrix nodeMap(const IntegerVector& node_lsfi1, const IntegerVector& node_lsfi2, const IntegerMatrix& cover_map){
  const int n = cover_map.nrow();
  std::vector<int> from_node_comb = std::vector<int>();
  std::vector<int> to_node_comb = std::vector<int>();
  for (int i = 0; i < n; ++i){
    const int from_idx = cover_map(i, 0), to_idx = cover_map(i, 1);
    IntegerVector from_nodes = which_cpp(node_lsfi1, from_idx);
    IntegerVector to_nodes = which_cpp(node_lsfi2, to_idx);
    for (IntegerVector::const_iterator f = from_nodes.begin(); f != from_nodes.end(); ++f){
      for (IntegerVector::const_iterator t = to_nodes.begin(); t != to_nodes.end(); ++t){
        from_node_comb.push_back(*f);
        to_node_comb.push_back(*t);
      }
    }
  }
  const int n2 = from_node_comb.size();
  IntegerMatrix res = IntegerMatrix(n2, 2);
  for (int i = 0; i < n2; ++i){
    res(i, _) = IntegerVector::create(from_node_comb[i], to_node_comb[i]);
  }
  return(res);
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
