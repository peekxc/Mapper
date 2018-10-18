#include <Rcpp.h>
using namespace Rcpp;

#include <iterator> // std::back_inserter
#include <chrono>
#include <algorithm> // std::remove
#include <queue> // std::priority_queue
#include "SegmentTree.h"

// Simple helper function to load an R function into C++
Rcpp::Function getFunctionR(std::string fname, std::string ns){
  Rcpp::Environment ns_env(ns); // Obtain environment containing function
  Rcpp::Function f = ns_env[fname]; // Make function callable from C++ 
  return(f);
}

// Merges a vector of sorted vectors into one sorted vector 
template <typename T> 
std::vector<T> sorted_merge(const std::vector< std::vector<T> >& vec){
  // Get the size of the output vector 
  std::vector< typename std::vector<T>::const_iterator > its; 
  std::size_t total_vec_size = 0;
  std::for_each(vec.begin(), vec.end(), [&](const std::vector<T>& v){
    total_vec_size += v.size();
    its.push_back(v.begin());
  });
  
  // Dereference the constant iterators, extract the minimum element (and position) of the iterators, and increment the 
  // iterator that had the smallest element. If an iterator reaches its end, return the maximum numerical limit of T.
  std::vector<T> elements = std::vector<T>(vec.size());
  auto get_min = [&its, &elements, &vec](){
    std::size_t cc = 0;
    std::transform(its.begin(), its.end(), elements.begin(), [&](typename std::vector<T>::const_iterator& v_it){
      return(v_it == vec[cc++].end() ? std::numeric_limits<T>::max() : *v_it );
    });
    std::size_t min_pos = std::distance(elements.begin(), std::min_element(elements.begin(),elements.end()));
    ++its[min_pos];
    return(elements[min_pos]);
  };
  
  // Apply the minimum operation 
  std::vector< T > final_res = std::vector< T >();
  final_res.reserve(total_vec_size);
  for (std::size_t i = 0; i < total_vec_size; ++i){
    final_res.push_back(get_min());
  }
  return(final_res);
}

template <typename T> 
std::vector<T> sorted_merge2(const std::vector< std::vector<T> >& vec){
  std::size_t total_vec_size = 0;
  std::for_each(vec.begin(), vec.end(), [&](const std::vector<T>& v){ total_vec_size += v.size(); });
  std::vector< T > final_res = std::vector< T >();
  final_res.reserve(total_vec_size);
  std::for_each(vec.begin(), vec.end(), [&](const std::vector<T>& v){ std::copy(v.begin(), v.end(), std::back_inserter(final_res)); });
  std::sort(final_res.begin(), final_res.end());
  return(final_res);
}




SegmentTree::SegmentTree(const NumericVector& endpoints) : n(endpoints.size()), h(ceil(log2(n + 1))) {
  if (endpoints.size() % 2 != 0){ stop("There must be an even number of endpoints to represent intervals."); }
  s_endpts = clone(endpoints).sort(); // store the (sorted) endpoints 
  o_endpts = match(s_endpts, endpoints); // save the original ordering for interval queries; TODO: pass ordering in
  tr = std::vector< std::vector<int> >(n); // only inner nodes
  leaves = std::vector< std::vector<int> >(n); // leaf nodes contain indices of points
  build(); // build the internal nodes
}
  
// Empty destructor; all containers should automatically deallocate
SegmentTree::~SegmentTree(){ }

SEXP SegmentTree::as_XPtr(){
  Rcpp::XPtr< SegmentTree> p(this, false); // does *not* register a finalizer
  return(p);
}
  
// Given the leaves, builds the internal nodes of the tree bottom-up
// Rcout << "left child " << lc << " is leaf, pushing leaf index " << lc << " to node " << i << std::endl;
void SegmentTree::build() {
  for (int i = n - 1; i > 0; --i) {
    const int lc = i << 1, rc = i << 1|1; // equivalent to: [2*i], [2*i+1]
    if (lc >= n){ tr.at(i).push_back(lc); } // left child is a leaf
    else { std::copy(tr.at(lc).begin(), tr.at(lc).end(), std::back_inserter(tr.at(i))); } // left child is an inner node 
    if (rc >= n){ tr.at(i).push_back(rc); } // right child is a leaf
    else { std::copy(tr.at(rc).begin(), tr.at(rc).end(), std::back_inserter(tr.at(i))); } // right child is an inner node 
  } 
} // end build
  
void SegmentTree::printTree(){

  Rcout << "Segment Tree statistics:" << std::endl;
  Rcout << n << " endpoints, " << n - 1 << " inner nodes" << std::endl; 
  Rcout << "total number of nodes: " << tr.size() << std::endl; 
  Rcout << "endpoints order: " << o_endpts << std::endl; 
}

// Maps a set of point coordinates to their corresponding intervals 
void SegmentTree::insert_pts(const NumericVector& pts){
 
  // Load R functions
  Function findInterval = getFunctionR("findInterval"), rle = getFunctionR("rle"), order = getFunctionR("order");
    
  // Sort the point / retrieve the original ordering
  IntegerVector o_pts = order(pts); // sorted ordering
  NumericVector s_pts = pts[o_pts - 1]; // sort coordinates
  
  // Partition the points into the given intervals
  IntegerVector pt_int = findInterval(_["x"] = s_pts, _["vec"] = s_endpts, _["rightmost.closed"] = false, _["all.inside"] = true);
  
  // Use run-length encoding to get the interval sizes
  const List rle_pt_int = rle(pt_int);
  const IntegerVector int_len = rle_pt_int["lengths"];
  const IntegerVector int_val = rle_pt_int["values"];
  
  // Assign the point indices to the leaves. The indices are sorted on copy.
  const int ne_int = int_len.size(); // number of non-empty intervals
  for (int i = 0, ci = 0; i < ne_int; ++i){
    const int int_index = int_val.at(i) - 1, interval_len = int_len.at(i);
    leaves.at(int_index).resize(interval_len);
    std::partial_sort_copy(o_pts.begin() + ci, o_pts.begin() + ci + interval_len, 
                           leaves.at(int_index).begin(), leaves.at(int_index).end());
    ci += interval_len;
  }
}
  
// Returns a list of all of the non-leaves in the tree
List SegmentTree::getInnerNodes(){
  List inner_nodes = List(n - 1);
  for (int i = 1; i < n; ++i){ inner_nodes.at(i - 1) = tr.at(i); }
  return(inner_nodes);
}
  
// Returns a list of all of the leaves in the tree
List SegmentTree::getLeafNodes(){ 
  return(wrap(leaves));
}
  
// Returns a list of all of the nodes in the tree
List SegmentTree::getAllNodes(){
  List nodes = List(tr.size());
  for (int i = 1; i < tr.size(); ++i){  nodes.at(i - 1) = tr.at(i); }
  return(nodes);
}

// Given a range of leaf indices [l, r], finds the leaf within that range containing 
// the given point index, or -1 otherwise. 
int SegmentTree::find_leaf(const int pt_idx, int l, int r){
  for (int i = l; i <= r; i++){
    std::vector<int>& leaf = leaves.at(i);
    if (std::find(leaf.begin(), leaf.end(), pt_idx) != leaf.end()){
      return(i);
    }
  }
  return(-1);
}

// Query an (index-converted) interval in [l, r)
// Returns the 0-based indices of the leaves which are contained in the interval.
std::vector<int> SegmentTree::query(int l, int r){ // query on 
  if (l < 0){ l = 0; }
  if (r >= n){ r = n - 1; }
  std::vector<int> nids = std::vector<int>(); // the node indices of the leaves containing the interval [l, r)
  nids.reserve(h);
  for (l += n, r += n; l < r; l >>= 1, r >>= 1) {
    if (l & 1) { // if l & 1, it is a right child, thus include l and swap to right of parent. Otherwise, move to l's parent.
      if (l >= n){  nids.push_back(l); }
      else { std::copy(tr[l].begin(), tr[l].end(), std::back_inserter(nids)); }
      ++l;
    }
    if (r & 1) { // if r & 1, it is a right child, thus include r and swap to left of parent. Otherwise, move to r's parent.
      if (r >= n){ nids.push_back(r-1); }
      else { std::copy(tr[r-1].begin(), tr[r-1].end(), std::back_inserter(nids)); }
      --r;
    }
  }
  return(nids);
}
  
// Queries on the interval [l, r).
// Unlike query(l, r) above, this function returns the original point indices of the points inserted into the tree. 
std::vector<int> SegmentTree::queryInterval(int l, int r){ 
  std::vector<int> query_res = query(l, r);
  std::size_t n_l = query_res.size(), total_length = 0;
  for (std::size_t i = 0; i < n_l; ++i) { total_length += leaves.at(query_res[i] - n).size(); }
  std::vector<int> output = std::vector<int>();
  output.reserve(total_length);
  for (std::size_t i = 0; i < n_l; ++i)
  {
    const std::vector<int>& el = leaves.at(query_res[i] - n);
    std::copy(el.begin(), el.end(), std::back_inserter(output));
  }
  return(output);
}

// std::vector<int> SegmentTree::queryInterval_sorted(int l, int r){ 
//   std::vector<int> query_res = query(l, r);
//   std::vector< std::vector<int> > all_queries; 
//   std::transform(query_res.begin(), query_res.end(), all_queries.begin(), [&](const int leaf_idx){
//     return(leaves.at(leaf_idx - n));
//   });
//   return(sorted_merge2(all_queries));
// }

// Find point <idx> in leaf <from>, swap it to leaf <to>
// TODO: make this much more efficient
void SegmentTree::swap(int from, int to, const int idx){
  std::vector<int>& from_segment = leaves.at(from);
  const int idx_exists = std::count(from_segment.begin(), from_segment.end(), idx);
  if (idx_exists == 0){
    Rprintf("point %d does not exist in from segment %d, trying to transer to segment %d", idx, from, to);
    stop("");
  }
  from_segment.erase(std::remove(from_segment.begin(), from_segment.end(), idx), from_segment.end());
  leaves.at(to).push_back(idx);
}

RCPP_MODULE(segment_tree_module) {
  Rcpp::class_<SegmentTree>("SegmentTree")
  .constructor<Rcpp::NumericVector>()
  .field( "s_endpts", &SegmentTree::s_endpts)
  .field( "o_endpts", &SegmentTree::o_endpts)
  .field( "tr", &SegmentTree::tr)
  .field_readonly( "n", &SegmentTree::n)
  .method( "as_XPtr", &SegmentTree::as_XPtr)
  .method( "printTree", &SegmentTree::printTree )
  .method( "getAllNodes", &SegmentTree::getAllNodes )
  .method( "getLeafNodes", &SegmentTree::getLeafNodes )
  .method( "getInnerNodes", &SegmentTree::getInnerNodes )
  .method( "insert_points", &SegmentTree::insert_pts )
  .method( "build", &SegmentTree::build )
  .method( "query", &SegmentTree::query )
  .method( "queryInterval", &SegmentTree::queryInterval )
  .method( "swap", &SegmentTree::swap)
  ;
}

// [[Rcpp::export]]
IntegerVector test_merge(const List& data){
  std::vector< std::vector<int> > my_vec;
  for (int i = 0; i < data.size(); ++i){
    IntegerVector v = data.at(i);
    my_vec.push_back(as< std::vector<int> >(v));
  }
  std::vector<int> res = sorted_merge(my_vec);
  return(wrap(res));
}

// [[Rcpp::export]]
IntegerVector test_merge2(const List& data){
  std::vector< std::vector<int> > my_vec;
  for (int i = 0; i < data.size(); ++i){
    IntegerVector v = data.at(i);
    my_vec.push_back(as< std::vector<int> >(v));
  }
  std::vector<int> res = sorted_merge2(my_vec);
  return(wrap(res));
}


/*** R
test_x <- replicate(15, { sort(as.integer(runif(10000)*10000)) }, simplify = F)
microbenchmark::microbenchmark({ invisible(Mapper:::test_merge(test_x)) })
microbenchmark::microbenchmark({ invisible(Mapper:::test_merge2(test_x)) })
microbenchmark::microbenchmark({ invisible(sort(unlist(test_x))) })

set.seed(1234)
# pts <- runif(26)
# intervals <- cbind(quantile(pts, c(0, 0.25, 0.50)), 
#                    quantile(pts, c(0.35, 0.65, 1)))
# s_endpts <- sort(as.vector(t(intervals)))
# pt_order <- order(pts)
# pt_cuts <- findInterval(pts[pt_order], vec = s_endpts, rightmost.closed = TRUE)
test_wut(3, 11, 16)

## Benchmarking 



## Testing multi segment tree 
set.seed(1234)
xy <- cbind(runif(10), runif(10))
xy_endpts <- cbind(c(0, 0.5, 0.5, 1), c(0, 0.5, 0.5, 1))

xy_endpts <- cbind(c(0, 0.5, 0.5, 1), c(0.0, 0.5, 0.5, 1))
y <- kdtools::kd_sort(xy_endpts)
kdtools::kd_range_query(x = y, l = c(0.0, 0.0), u = c(0.5, 0.5))

mstree <- Mapper:::MultiSegmentTree$new(xy_endpts)
mstree$insert_pts(xy)
mstree$query(c(0, 0), c(2, 2))

plot(xy, xlim = range(xy[, 1])+c(-0.5, 0.5), ylim=range(xy[, 2])+c(-0.5, 0.5))
text(xy, labels = 1:10, pos = 3)
abline(v = xy_endpts[, 1], col = "red")
abline(h = xy_endpts[, 2], col = "red")

stree <- Mapper::segment_tree(xy_endpts[, 1])
stree$insert_points(xy[,1])
stree$query(0, 6) # [l, r)


new_pts <- runif(n = 10000, min = 0, max = 16)
new_endpts <- sort(runif(16, min = 0, max = 16))
new_cuts <- findInterval(x = sort(new_pts), vec = new_endpts, all.inside = TRUE, rightmost.closed = TRUE)
  
stree <- Mapper::segment_tree(new_endpts)
stree$insert_points(new_pts)
stree$query(0, 6) # [l, r)

## Generate interval queries  
n_intervals <- length(new_endpts)
q <- matrix(sample(0:(n_intervals-1L), size = 1000*2, replace = TRUE), ncol = 2)
q <- cbind(pmin(q[, 1], q[, 2]), pmax(q[, 1], q[, 2])) # 0-based 

I <- q[3, ]
print(sprintf("Querying pts between: [%f, %f)", new_endpts[I[1]+1L], new_endpts[I[2]+1L]))
print(sort(new_pts[which((new_pts >= new_endpts[I[1]+1]) & (new_pts < new_endpts[I[2]+1]))]))
print(sort(new_pts[stree$queryInterval(I[1], I[2])]))

## Segment Tree allows the following query 
aug_endpts <- new_endpts
aug_endpts[c(1, length(aug_endpts))] <- c(-Inf, Inf)
microbenchmark::microbenchmark({ 
  apply(q, 1, function(I){ 
    invisible(which((new_pts >= aug_endpts[I[1]+1]) & (new_pts < aug_endpts[I[2]+1]))) 
  })
})

## Using sorted ranges w/ kdtools
library("kdtools")
x <- matrix(new_pts)
y <- matrix_to_tuples(x)
kdtools::kd_sort(x = y, inplace = TRUE)
microbenchmark::microbenchmark({ 
  apply(q, 1, function(I){ 
    invisible(kdtools::kd_range_query(y, l = aug_endpts[I[1]+1], u = aug_endpts[I[2]+1]))
  })
})

## Using ANN 
ANN::

## Segment Tree query equivalent
microbenchmark::microbenchmark({ 
  apply(q, 1, function(I){ invisible(stree$queryInterval(I[1], I[2])) })
})


## Testing multisegment tree 
new_pts <- replicate(2, runif(n = 10000, min = 0, max = 16), simplify = TRUE)
new_endpts <- matrix(rep(sort(runif(16, min = 0, max = 16)), 2), ncol = 2)
aug_endpts <- new_endpts
aug_endpts[1,] <- -Inf
aug_endpts[nrow(aug_endpts),] <- Inf

## The range queries
q <- matrix(sample(0:(nrow(new_endpts)-1L), size = 1000*2, replace = TRUE), ncol = 2)
q <- cbind(pmin(q[, 1], q[, 2]), pmax(q[, 1], q[, 2])) # 0-based 

## Testing kdtools 
x <- matrix(new_pts, ncol = 2)
y <- matrix_to_tuples(x)
kdtools::kd_sort(x = y, inplace = TRUE)
microbenchmark::microbenchmark({ 
  apply(q, 1, function(I){ 
    invisible(kdtools::kd_range_query(y, l = aug_endpts[rep(I[1]+1, 2)], u = aug_endpts[rep(I[2]+1, 2)]))
  })
})

## Testing multisegment tree
mstree <- Mapper:::MultiSegmentTree$new(new_endpts)
mstree$insert_pts(new_pts)
microbenchmark::microbenchmark({ 
  apply(q, 1, function(I){ invisible(mstree$query(rep(I[1], 2), rep(I[2], 2))) })
})

# 
# microbenchmark::microbenchmark({ 
#   apply(q, 1, function(I){ invisible(which(new_cuts %in% (I[1]+1):(I[2]+1))) })
# })
# 
# ## Testing out swapping 
# X <- c(x_1=0.2, x_2=4/3, x_3=1.6)
# endpts <- c(0, 1, 1, 2)
# stree <- Mapper::segment_tree(endpts)
# stree$insert_points(X)
# stree$queryInterval(4-4, 6-4)
# stree$queryInterval(4-4, 7-4) 
# 
# stree$swap(6-4, 5-4, 2)
# 
# stree$queryInterval(4-4, 6-4)
# stree$queryInterval(6-4, 7-4) 
# stree$queryInterval(5-4, 7-4) 
*/
