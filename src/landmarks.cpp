#include <Rcpp.h>
using namespace Rcpp;

#include <functional>
#include <algorithm>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using std::size_t;

inline double sq_dist(const NumericVector& x, const NumericVector& y){
  NumericVector diff = x-y;
  return(sum(pow(diff, 2.0)));
}


// Eccentricity
// [[Rcpp::export]]
NumericVector eccentricity(const NumericMatrix& from, const NumericMatrix& x, const int type = 1, bool display_progress=true) {

  // Compute eccentricity
  const size_t n_pts = x.nrow();
  const size_t n_src_pts = from.nrow();
  Progress p(n_src_pts, display_progress);
  std::vector< double > pt_ecc(n_src_pts);
  std::vector< double > pt_dist(n_pts);
  for (size_t i = 0; i < n_src_pts; ++i){
    if (Progress::check_abort()){ p.cleanup(); return(NumericVector::create(-1.0)); }
    for (size_t j = 0; j < n_pts; ++j){
      pt_dist.at(j) = sq_dist(from.row(i), x.row(j));
    }
    p.increment();
    pt_ecc.at(i) = accumulate(begin(pt_dist), end(pt_dist), 0.0, std::plus< double >())/n_pts;
  }
  return(wrap(pt_ecc));
}

// Uses the euclidean maxmin procedure to choose n landmarks. 
// [[Rcpp::export]]
IntegerVector landmark_maxmin(const NumericMatrix& x, const int n, const int seed = 0, bool display_progress=true) {
  const size_t n_pts = x.nrow();
  std::vector< double > lm_dist(n_pts, std::numeric_limits<double>::infinity());
  IntegerVector landmark_idx = no_init_vector(n);
  landmark_idx[seed] = 0; 
  IntegerVector::iterator c_landmark = landmark_idx.begin();
  double new_min_dist;
  Progress p(n, display_progress);
  std::generate(landmark_idx.begin()+1, landmark_idx.end(), [&](){
    size_t i = 0;
    new_min_dist = std::numeric_limits<double>::infinity();
    
    // Progress bar
    if (Progress::check_abort()) { p.cleanup(); Rcpp::stop("stopping"); }
    
    // Replace nearest landmark distances
    std::replace_if(lm_dist.begin(), lm_dist.end(), [&i, &x, &c_landmark, &new_min_dist](const double c_dist){
      new_min_dist = sq_dist(x.row(*c_landmark), x.row(i++));
      return(new_min_dist < c_dist);
    }, new_min_dist);
    
    // Find the point that maximizes said distances, move to next landmark
    c_landmark++;
    p.increment();
    auto max_landmark = std::max_element(lm_dist.begin(), lm_dist.end());
    return(std::distance(lm_dist.begin(), max_landmark));
  });
  return(landmark_idx+1); // 1-based return
}
