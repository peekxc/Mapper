#include <Rcpp.h>
using namespace Rcpp;

#include <functional>
#include <algorithm>

using std::vector;
using std::begin; 
using std::end;
using std::size_t;


inline double max_dist(const NumericVector& x, const NumericVector& y){
  return(max(abs(x-y)));
}
  
inline double man_dist(const NumericVector& x, const NumericVector& y){
  NumericVector diff = abs(x-y);
  return(sum(diff));
}

inline double sq_dist(const NumericVector& x, const NumericVector& y){
  NumericVector diff = x-y;
  return(sqrt(sum(pow(diff, 2.0))));
}

// Chooses a distance metric form one of the common p-norms
auto choose_dist(size_t metric) -> std::function<double(const NumericVector&, const NumericVector&)> 
{
  if (metric == 0){ return(sq_dist); }
  else if (metric == 1){ return(man_dist); }
  else if (metric == 2){ return(max_dist); }
  return(sq_dist);
}

// Eccentricity
// Computes a distance measure from the barycenter of the data, or from a given center, if provided. 
// [[Rcpp::export]]
NumericVector eccentricity(const NumericMatrix& from, const NumericMatrix& x, const int type = 1) {
  if (from.ncol() != x.ncol()){ stop("Matrices must have identical number of columns."); }
  const size_t n_pts = x.nrow();
  const size_t n_src_pts = from.nrow();
  std::vector< double > pt_ecc(n_src_pts);
  std::vector< double > pt_dist(n_pts);
  for (size_t i = 0; i < n_src_pts; ++i){
    for (size_t j = 0; j < n_pts; ++j){
      pt_dist.at(j) = sq_dist(from.row(i), x.row(j));
    }
    pt_ecc.at(i) = accumulate(begin(pt_dist), end(pt_dist), 0.0, std::plus< double >())/n_pts;
  }
  return(wrap(pt_ecc));
}

// Uses the euclidean maxmin procedure to choose n landmarks. 
// x := point cloud (each row is one point)
// n := number of landmarks requested
// seed := initial point (default is point at index 0)
// [[Rcpp::export]]
IntegerVector maxmin_n(const NumericMatrix& x, const int n, const size_t metric = 0, const int seed = 0) {
  const size_t n_pts = x.nrow();
  if (seed < 0 || seed >= n_pts){ stop("Invalid seed point."); }
  
  // Get the distance function
  const auto dist_f = choose_dist(metric);
  
  // Preallocate distance vector for landmarks
  std::vector< double > lm_dist(n_pts, std::numeric_limits<double>::infinity());
  
  // Choose the initial landmark
  vector< size_t > lm = vector< size_t >();
  lm.reserve(n);
  lm.push_back(seed);
  
  // Generate the landmarks
  double new_min_dist = std::numeric_limits<double>::infinity();
  std::generate_n(back_inserter(lm), n-1, [&](){
    
    // Inductively, replace point-to-landmark distances if lower than previously computed
    size_t c_lm = lm.back(), i = 0;
    std::replace_if(begin(lm_dist), end(lm_dist), [&i, &dist_f, &x, &c_lm, &new_min_dist](const double c_dist){
      // new_min_dist = as< double >( dist(x.row(c_lm), x.row(i++)) );
      new_min_dist = dist_f(x.row(c_lm), x.row(i++));
      return(new_min_dist < c_dist);
    }, new_min_dist);
    
    // Find the point that maximizes said distances, move to next landmark
    auto max_landmark = std::max_element(begin(lm_dist), end(lm_dist));
    return(std::distance(begin(lm_dist), max_landmark));
  });
    
  // Wrap in R-memory and return 1-based indices
  IntegerVector res_lm = wrap(lm);
  return(res_lm+1);
}

// Uses the euclidean maxmin procedure to choose n landmarks. 
// x := point cloud (each row is one point)
// n := number of landmarks requested
// dist_f := distance function
// seed := initial point (default is point at index 0)
// [[Rcpp::export]]
IntegerVector maxmin_n_f(const NumericMatrix& x, const int n, const Function& dist_f, const int seed = 0) {
  const size_t n_pts = x.nrow();
  if (seed < 0 || seed >= n_pts){ stop("Invalid seed point."); }
  
  // Preallocate distance vector for landmarks
  std::vector< double > lm_dist(n_pts, std::numeric_limits<double>::infinity());
  
  // Choose the initial landmark
  vector< size_t > lm = vector< size_t >();
  lm.reserve(n);
  lm.push_back(seed);
  
  // Generate the landmarks
  double new_min_dist = std::numeric_limits<double>::infinity();
  std::generate_n(back_inserter(lm), n-1, [&](){
    
    // Inductively, replace point-to-landmark distances if lower than previously computed
    size_t c_lm = lm.back(), i = 0;
    std::replace_if(begin(lm_dist), end(lm_dist), [&i, &dist_f, &x, &c_lm, &new_min_dist](const double c_dist){
      // new_min_dist = as< double >( dist(x.row(c_lm), x.row(i++)) );
      new_min_dist = as< double >( dist_f(x.row(c_lm), x.row(i++)) );
      return(new_min_dist < c_dist);
    }, new_min_dist);
    
    // Find the point that maximizes said distances, move to next landmark
    auto max_landmark = std::max_element(begin(lm_dist), end(lm_dist));
    return(std::distance(begin(lm_dist), max_landmark));
  });
  
  // Wrap in R-memory and return 1-based indices
  IntegerVector res_lm = wrap(lm);
  return(res_lm+1);
}

// x := point cloud (each row is one point)
// eps := distance to stop finding landmarks at
// seed := initial point (default is point at index 0)
// [[Rcpp::export]]
IntegerVector maxmin_eps(const NumericMatrix& x, const double eps, const size_t metric = 0, const int seed = 0) {
  const size_t n_pts = x.nrow();
  if (seed < 0 || seed >= n_pts){ stop("Invalid seed point."); }
  
  // Get the distance function
  const auto dist_f = choose_dist(metric);
  
  // Preallocate distance vector for landmarks
  std::vector< double > lm_dist(n_pts, std::numeric_limits<double>::infinity());
  
  // Choose the initial landmark
  vector< size_t > lm = vector< size_t >();
  lm.push_back(seed);
  
  // Generate the landmarks
  double new_min_dist = std::numeric_limits<double>::infinity();
  while (true){
    // Inductively, replace point-to-landmark distances if lower than previously computed
    size_t c_lm = lm.back(), i = 0;
    std::replace_if(begin(lm_dist), end(lm_dist), [&i, &x, &dist_f, &c_lm, &new_min_dist](const double c_dist){
      new_min_dist = dist_f(x.row(c_lm), x.row(i++));
      return(new_min_dist < c_dist);
    }, new_min_dist);
    
    // Find the point that maximizes said distances, move to next landmark
    auto max_landmark = std::max_element(begin(lm_dist), end(lm_dist));
    if (*max_landmark > eps){
      lm.push_back(std::distance(begin(lm_dist), max_landmark));
    }
    else { break; }
  }
  
  // Wrap in R-memory and return 1-based indices
  IntegerVector res_lm = wrap(lm);
  return(res_lm+1);
}

// x := point cloud (each row is one point)
// eps := distance to stop finding landmarks at
// dist_f := distance function
// seed := initial point (default is point at index 0)
// [[Rcpp::export]]
IntegerVector maxmin_eps_f(const NumericMatrix& x, const double eps, const Function& dist_f, const int seed = 0) {
  const size_t n_pts = x.nrow();
  if (seed < 0 || seed >= n_pts){ stop("Invalid seed point."); }
  
  // Preallocate distance vector for landmarks
  std::vector< double > lm_dist(n_pts, std::numeric_limits<double>::infinity());
  
  // Choose the initial landmark
  vector< size_t > lm = vector< size_t >();
  lm.push_back(seed);
  
  // Generate the landmarks
  double new_min_dist = std::numeric_limits<double>::infinity();
  while (true){
    // Inductively, replace point-to-landmark distances if lower than previously computed
    size_t c_lm = lm.back(), i = 0;
    std::replace_if(begin(lm_dist), end(lm_dist), [&i, &x, &dist_f, &c_lm, &new_min_dist](const double c_dist){
      new_min_dist = as< double >( dist_f(x.row(c_lm), x.row(i++)) );
      return(new_min_dist < c_dist);
    }, new_min_dist);
    
    // Find the point that maximizes said distances, move to next landmark
    auto max_landmark = std::max_element(begin(lm_dist), end(lm_dist));
    if (*max_landmark > eps){
      lm.push_back(std::distance(begin(lm_dist), max_landmark));
    }
    else { break; }
  }
  
  // Wrap in R-memory and return 1-based indices
  IntegerVector res_lm = wrap(lm);
  return(res_lm+1);
}


