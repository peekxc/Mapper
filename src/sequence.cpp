#include <Rcpp.h>
using namespace Rcpp;

using std::vector;

// Assuming lst contains numeric vectors not necessarily of equal length, generates a matrix 
// containing values from each element in 'lst' that are nondecreasing 
// [[Rcpp::export]]
List nondecreasing_seq(const List& lst){
  using std::vector; 
  using std::size_t;
  
  // Dimension 
  const size_t d = lst.size();
  vector< size_t > max_sizes = vector< size_t >(d);
  for (size_t d_i = 0; d_i < d; ++d_i){
    NumericVector eps = as< NumericVector >(lst.at(d_i));
    max_sizes.at(d_i) = eps.size();
  }
  
  // Start off at index -1 across all the dimensions
  vector< size_t > idx = vector< size_t >(d, 0);
  
  // While not at the very last index of all the lists
  const auto finished = [&d, &max_sizes, &idx]() -> bool {
    bool is_finished = true; 
    for (size_t d_i = 0; d_i < d; ++d_i){
      is_finished &= (max_sizes.at(d_i) == (idx.at(d_i)+1));
    }
    return is_finished;
  };
  vector< vector< double > > eps_res; 
  vector< vector< size_t > > idx_res; 
  while (!finished()){
    
    // Populate current eps + idx values
    vector< double > c_eps = vector< double >(d);
    size_t d_i = 0; 
    std::transform(begin(idx), end(idx), begin(c_eps), [&lst, &d_i](const size_t ii){
      const NumericVector& eps = as< NumericVector >(lst.at(d_i++));
      return eps.at(ii);
    });
    
    // Record them
    eps_res.push_back(c_eps);
    idx_res.push_back(vector< size_t >(begin(idx), end(idx)));
    
    // Get next eps value
    vector< double > new_eps = vector< double >(d);
    d_i = 0; 
    std::transform(begin(idx), end(idx), begin(new_eps), [&lst, &d_i](const size_t ii){
      const NumericVector& eps = as< NumericVector >(lst.at(d_i++));
      if (ii+1 >= eps.size()){ return std::numeric_limits< double >::infinity(); }
      return eps.at(ii+1);
    });
    auto min_el = std::min_element(begin(new_eps), end(new_eps));
    size_t el_idx = std::distance(begin(new_eps), min_el);
    idx.at(el_idx) += 1; 
  }
  
  const size_t n = eps_res.size();
  NumericMatrix eps_out(n, d);
  IntegerMatrix idx_out(n, d);
  for (size_t i = 0; i < n; ++i){
    NumericVector c_eps = wrap(eps_res.at(i));
    IntegerVector c_idx = wrap(idx_res.at(i));
    eps_out(i, _) = c_eps;
    idx_out(i, _) = c_idx;
  }
  return(List::create(_["eps"] = eps_out, _["idx"] = idx_out));
}
