#include "MultiSegmentTree.h"

MultiSegmentTree::MultiSegmentTree(const List& endpoints) : d(endpoints.size()) {
  for (std::size_t i = 0; i < d; ++i){
    NumericVector epts_di = endpoints.at(i);
    n = epts_di.size();
    std::shared_ptr<SegmentTree> sptr(new SegmentTree(epts_di));
    tr.push_back(sptr);
  }
}

// empty destructor
MultiSegmentTree::~MultiSegmentTree(){ }
  
void MultiSegmentTree::insert_pts(const NumericMatrix& pts){
  for (std::size_t i = 0; i < d; ++i){
    NumericVector pts_di = pts.column(i);
    tr.at(i)->insert_pts(pts_di);
  }
}
  
// Query the trees
IntegerVector MultiSegmentTree::query(const IntegerVector& l, const IntegerVector& r){
  if (l.size() != d || r.size() != d){ stop("Invalid dimensionality of requested query bounds."); }
  std::vector< int > final_res;
  std::vector<int>::iterator it;
  for (std::size_t i = 0; i < d; ++i){
    std::vector<int> leaf_res = tr.at(i)->queryInterval(l[i], r[i]);
    std::sort(leaf_res.begin(), leaf_res.end());
    if (i == 0){
      std::copy(leaf_res.begin(), leaf_res.end(), std::back_inserter(final_res));
    } else {
      it = std::set_intersection(final_res.begin(), final_res.end(), leaf_res.begin(), leaf_res.end(), final_res.begin());
      final_res.resize(it-final_res.begin());
    }
  }
  // return(List::create(_["all_queries"] = wrap(prelim_res), _["intersection"] = wrap(final_res)));
  return(wrap(final_res));
}

IntegerVector MultiSegmentTree::query_dimension(const int l, const int r, const int d_i){
  if (d_i < 0 || d_i >= d){ stop("Invalid dimension."); }
  return(wrap(tr.at(d_i)->queryInterval(l, r)));
}


RCPP_MODULE(multi_segment_tree_module) {
  Rcpp::class_<MultiSegmentTree>("MultiSegmentTree")
  .constructor<Rcpp::List>()
  .method( "insert_pts", &MultiSegmentTree::insert_pts)
  .method( "query", &MultiSegmentTree::query)
  .method( "query_dimension", &MultiSegmentTree::query_dimension)
  ;
}
