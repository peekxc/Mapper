#include <Rcpp.h>
using namespace Rcpp;

#include "SegmentTree.h"

struct MultiSegmentTree {
  std::vector< std::shared_ptr<SegmentTree> > tr; // the trees composing the tree (multidimensional case only)
  const int d;
  int n;
  MultiSegmentTree(const List& endpoints);
  ~MultiSegmentTree();
  void insert_pts(const NumericMatrix& pts);
  IntegerVector query(const IntegerVector& l, const IntegerVector& r);
  IntegerVector query_dimension(const int l, const int r, const int d_i);
};
