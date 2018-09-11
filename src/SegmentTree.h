#include <Rcpp.h>
using namespace Rcpp;

#include <iterator> // std::back_inserter
#include <chrono>
#include <algorithm> // std::remove

// Simple helper function to load an R function into C++
Rcpp::Function getFunctionR(std::string fname, std::string ns = "package:base");

// A segment tree stores a balanced binary tree of n endpoints and (n - 1) intervals.
// This specific segment tree class is specialized to handle interval queries which return the integer indices of points within a given query interval.
// The following code used for the building and query functions is based on Ai.Cash's blog post: https://codeforces.com/blog/entry/18051 
//' @export SegmentTree
struct SegmentTree {
  NumericVector s_endpts; // sorted endpoints
  IntegerVector o_endpts; // original order of the unsorted endpoints
  std::vector< std::vector<int> > tr; // the nodes composing the tree
  std::vector< std::vector<int> > leaves; // the leaves
  const int n, h; // number of endpoints, height of tree
  
  SegmentTree(const NumericVector& endpoints);
  ~SegmentTree();
  SEXP as_XPtr();
  void build();
  void printTree();
  void insert_pts(const NumericVector& pts);
  List getInnerNodes();
  List getLeafNodes();
  List getAllNodes();
  std::vector<int> query(int l, int r);
  std::vector<int> queryInterval(int l, int r);
  int find_leaf(const int pt_idx, int l, int r);
  void swap(int from, int to, const int idx);
};
