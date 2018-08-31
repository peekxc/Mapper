#include <Rcpp.h>
using namespace Rcpp;

#include <iterator> // std::back_inserter
#include <chrono>
// Simple helper function to load an R function into C++
Rcpp::Function getFunctionR(std::string fname, std::string ns = "package:base"){
  Rcpp::Environment ns_env(ns); // Obtain environment containing function
  Rcpp::Function f = ns_env[fname]; // Make function callable from C++ 
  return(f);
}


// A segment tree stores a balanced binary tree of n endpoints and (n - 1) intervals.
// This specific segment tree class is specialized to handle interval queries which return the integer indices of points within a given query interval.
// The following code used for the building and query functions is based on Ai.Cash's blog post: https://codeforces.com/blog/entry/18051 
//' @export SegmentTree
struct SegmentTree {
  NumericVector s_endpts; // sorted endpoints
  IntegerVector o_endpts; // original order of the unsorted endpoints
  std::vector< std::vector<int> > tr; // the nodes composing the tree
  std::vector< std::vector<int> > leaves; // the leaves
  // List tr; // the nodes composing the tree
  const int n; // number of endpoints
  
  // Construct the tree
  SegmentTree(const NumericVector& endpoints) : n(endpoints.size()) { // , n(N - 1), h(ceil(log2(n+1)))
    s_endpts = clone(endpoints).sort(); // store the (sorted) endpoints 
    o_endpts = match(s_endpts, endpoints); // save the original ordering for interval queries; TODO: pass ordering in
    // tr = List(2 * n); // allocate memory for the tree
    tr = std::vector< std::vector<int> >(n); // only inner nodes
    leaves = std::vector< std::vector<int> >(n);
    build(); // build the internal nodes
  }
  
  // empty destructor 
  ~SegmentTree(){ }
  
  // Given the leaves, builds the internal nodes of the tree bottom-up
  // Rcout << "left child " << lc << " is leaf, pushing leaf index " << lc << " to node " << i << std::endl;
  void build() {
    for (int i = n - 1; i > 0; --i) {
      const int lc = i << 1, rc = i << 1|1; // equivalent to: [2*i], [2*i+1]
      if (lc >= n){ tr.at(i).push_back(lc); } // left child is a leaf
      else { std::copy(tr.at(lc).begin(), tr.at(lc).end(), std::back_inserter(tr.at(i))); } // left child is an inner node 
      if (rc >= n){ tr.at(i).push_back(rc); } // right child is a leaf
      else { std::copy(tr.at(rc).begin(), tr.at(rc).end(), std::back_inserter(tr.at(i))); } // right child is an inner node 
    } 
  } // end build
  
  void printTree(){
    Rcout << "Segment Tree statistics:" << std::endl;
    Rcout << n << " endpoints, " << n - 1 << " inner nodes" << std::endl; 
    Rcout << "total number of nodes: " << tr.size() << std::endl; 
    Rcout << "endpoints order: " << o_endpts << std::endl; 
  }

  // Maps a set of point coordinates to their corresponding intervals 
  void insert_pts(const NumericVector& pts){
   
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
    
    // Assign the point indices to the leaves
    const int ne_int = int_len.size(); // number of non-empty intervals
    for (int i = 0, ci = 0; i < ne_int; ++i){
      // const int int_index = n + int_val.at(i), interval_len = int_len.at(i);
      // std::vector<int>& interval = tr.at(int_index);
      // interval.reserve(interval_len);
      // std::copy(o_pts.begin() + ci, o_pts.begin() + ci + interval_len, std::back_inserter(interval));
      // ci += interval_len;
      const int int_index = int_val.at(i) - 1, interval_len = int_len.at(i);
      std::copy(o_pts.begin() + ci, o_pts.begin() + ci + interval_len, std::back_inserter(leaves[int_index]));
      ci += interval_len;
    }

  }
  
  // Returns a list of all of the non-leaves in the tree
  List getInnerNodes(){
    List inner_nodes = List(n - 1);
    for (int i = 1; i < n; ++i){ inner_nodes.at(i - 1) = tr.at(i); }
    return(inner_nodes);
  }
  
  // Returns a list of all of the leaves in the tree
  List getLeafNodes(){
    // List leaves = List(n);
    // for (int i = 0; i < n; ++i){ leaves.at(i) = tr.at(n + i); }
    return(wrap(leaves));
  }
  
  // Returns a list of all of the nodes in the tree
  List getAllNodes(){
    List nodes = List(tr.size());
    for (int i = 1; i < tr.size(); ++i){  nodes.at(i - 1) = tr.at(i); }
    return(nodes);
  }

  // Query an (index-converted) interval in [l, r)
  // Returns the 0-based indices of the leaves which are contained in the interval.
  std::vector<int> query(int l, int r){ // query on 
    if (l < 0){ l = 0; }
    if (r >= n){ r = n - 1; }
    // auto start = std::chrono::steady_clock::now();
    std::vector<int> nids = std::vector<int>(); // the node indices of the leaves containing the interval [l, r)
    nids.reserve(int(log2(n)) + 1);
    for (l += n, r += n; l < r; l >>= 1, r >>= 1) {
      // Rcout << "l=" << l <<  std::endl; 
      if (l & 1) { // if l & 1, it is a right child, thus include l and swap to right of parent. Otherwise, move to l's parent.
        if (l >= n){ 
          // Rcout << "Pushing l=" << l << std::endl; 
          nids.push_back(l); 
        }
        else { 
          // Rcout << "Copying l=" << l << std::endl; 
          std::copy(tr[l].begin(), tr[l].end(), std::back_inserter(nids)); 
        }
        ++l;
      }
      if (r & 1) { // if r & 1, it is a right child, thus include r and swap to left of parent. Otherwise, move to r's parent.
        if (r >= n){ 
          // Rcout << "Pushing r=" << r-1 << std::endl; 
          nids.push_back(r-1); 
        }
        else { 
          // Rcout << "Copying r=" << r-1 << std::endl; 
          std::copy(tr[r-1].begin(), tr[r-1].end(), std::back_inserter(nids)); 
        }
        --r;
      }
    }
    return(nids);
    // auto end = std::chrono::steady_clock::now();
    // auto diff1 = end - start;
    // start = std::chrono::steady_clock::now();
    // std::sort(nids.begin(), nids.end()); // return points in order
    // std::vector<int> pt_ids = std::vector<int>();
    // for (int i = 0; i < nids.size(); ++i){
    //   std::copy(tr.at(nids.at(i)).begin(), tr.at(nids.at(i)).end(), std::back_inserter(pt_ids));
    // }
    // end = std::chrono::steady_clock::now();
    // auto diff2 = end - start;
    
    // std::cout << std::chrono::duration <double, std::nano> (diff1).count() << " ns" << std::endl;
    // std::cout << std::chrono::duration <double, std::nano> (diff2).count() << " ns" << std::endl;
    // double t1 = std::chrono::duration <double, std::milli> (diff1).count();
    // double t2 = std::chrono::duration <double, std::milli> (diff2).count();
    // std::vector<double> ret = { t1, t2 };
    // return ret;
  }
  
  // Queries on the interval [l, r).
  // Unlike query(l, r) above, this function returns the original point indices of the points inserted into the tree. 
  IntegerVector queryInterval(int l, int r){ 
    std::vector<int> query_res = query(l, r);
    std::size_t n_l = query_res.size(), total_length = 0, index = 0;;
    for (std::size_t i = 0; i < n_l; ++i) total_length += leaves.at(query_res[i] - n).size();
    IntegerVector output = no_init(total_length);
    for (std::size_t i = 0; i < n_l; ++i)
    {
      const std::vector<int>& el = leaves.at(query_res[i] - n);
      std::copy(el.begin(), el.end(), output.begin() + index);
      index += el.size();
    }
    return output;
  }
  
  // std::vector<int> queryInterval2(int l, int r){ 
  //   std::vector<int> res = std::vector<int>(); // the node indices of the leaves containing the interval [l, r)
  //   for (l += n, r += n; l < r; l >>= 1, r >>= 1) {
  //     //Rcout << "(l, r) := (" << l << "," << r << ")" << std::endl; 
  //     if (l & 1) { // if l & 1, it is a right child, thus include l and swap to right of parent. Otherwise, move to l's parent.
  //       if (l >= n){ std::copy(tr.at(l).begin(), tr.at(l).end(), std::back_inserter(res)); } 
  //       else { 
  //         for (int i = 0; i < tr.at(l).size(); ++i){
  //           std::copy(tr.at(tr.at(l).at(i)).begin(), tr.at(tr.at(l).at(i)).end(), std::back_inserter(res));
  //         }
  //       }
  //       ++l;
  //     }
  //     if (r & 1) { // if r & 1, it is a left child, thus include r and swap to left of parent. Otherwise, move to r's parent.
  //       if (r >= n){ //res.push_back(r); 
  //         std::copy(tr.at(r).begin(), tr.at(r).end(), std::back_inserter(res));
  //       } else { 
  //         for (int i = 0; i < tr.at(r-1).size(); ++i){
  //           std::copy(tr.at(tr.at(r-1).at(i)).begin(), tr.at(tr.at(r-1).at(i)).end(), std::back_inserter(res));
  //         }
  //       r--;
  //       }
  //     }
  //   }
  // }
  
  
  // std::vector<int> nodeToPoints(std::vector<int> nids){
  // 
  //   return(pt_ids);
  // }

};

// struct MultiSegmentTree {
//   std::vector< SegmentTree* > tr; // the trees composing the tree (multidimensional case only)
//  
//   std::vector< std::vector<int> > leaves; // the leaves
//   // List tr; // the nodes composing the tree
//   const int d, n;
//   
//   MultiSegmentTree(const NumericMatrix& endpoints) : d(endpoints.ncol()), n(endpoints.nrow()){
//     NumericVector ep_x = endpoints.column(0);
//       s_ep_x = clone(ep_x).sort(); // store the (sorted) endpoints 
//       o_ep_x = match(s_endpts, endpoints); // save the original ordering for interval queries; TODO: pass ordering in
//       // tr = List(2 * n); // allocate memory for the tree
//       tr = std::vector< SegmentTree* >(n); // only inner nodes
//       leaves = std::vector< std::vector<int> >(n);
//       build(); // build the internal nodes
//   }
//   
//   // empty destructor 
//   ~MultiSegmentTree(){ }
//   
//   void build() {
//     for (int i = n - 1; i > 0; --i) {
//       const int lc = i << 1, rc = i << 1|1; // equivalent to: [2*i], [2*i+1]
//       if (lc >= n){ tr.at(i) .push_back(lc); } // left child is a leaf
//       else { std::copy(tr.at(lc).begin(), tr.at(lc).end(), std::back_inserter(tr.at(i))); } // left child is an inner node 
//       if (rc >= n){ tr.at(i).push_back(rc); } // right child is a leaf
//       else { std::copy(tr.at(rc).begin(), tr.at(rc).end(), std::back_inserter(tr.at(i))); } // right child is an inner node 
//     } 
//   } // end build
// };

RCPP_MODULE(st_module) {
  Rcpp::class_<SegmentTree>("SegmentTree")
  .constructor<Rcpp::NumericVector>()
  .field( "s_endpts", &SegmentTree::s_endpts)
  .field( "o_endpts", &SegmentTree::o_endpts)
  .field( "tr", &SegmentTree::tr)
  .field_readonly( "n", &SegmentTree::n)
  .method( "printTree", &SegmentTree::printTree )
  .method( "getAllNodes", &SegmentTree::getAllNodes )
  .method( "getLeafNodes", &SegmentTree::getLeafNodes )
  .method( "getInnerNodes", &SegmentTree::getInnerNodes )
  .method( "insert_points", &SegmentTree::insert_pts )
  .method( "build", &SegmentTree::build )
  .method( "query", &SegmentTree::query )
  .method( "queryInterval", &SegmentTree::queryInterval )
  ;
}



/*** R
set.seed(1234)
# pts <- runif(26)
# intervals <- cbind(quantile(pts, c(0, 0.25, 0.50)), 
#                    quantile(pts, c(0.35, 0.65, 1)))
# s_endpts <- sort(as.vector(t(intervals)))
# pt_order <- order(pts)
# pt_cuts <- findInterval(pts[pt_order], vec = s_endpts, rightmost.closed = TRUE)
test_wut(3, 11, 16)




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

pt_order <- order(new_pts)



aug_endpts <- new_endpts
aug_endpts[c(1, length(aug_endpts))] <- c(-Inf, Inf)
wut1 <- apply(q, 1, function(I){ 
  invisible(which((new_pts >= aug_endpts[I[1]+1]) & (new_pts < aug_endpts[I[2]+1]))) 
})

wut2 <- apply(q, 1, function(I){ 
  stree$queryInterval(I[1], I[2])
}) 

which(sapply(wut1, length) == sapply(wut2, length))

microbenchmark::microbenchmark({ 
  apply(q, 1, function(I){ 
    invisible(which((new_pts >= aug_endpts[I[1]+1]) & (new_pts < aug_endpts[I[2]+1]))) 
  })
})
microbenchmark::microbenchmark({ 
  apply(q, 1, function(I){ invisible(stree$queryInterval(I[1], I[2])) })
})


microbenchmark::microbenchmark({ 
  apply(q, 1, function(I){ invisible(which(new_cuts %in% (I[1]+1):(I[2]+1))) })
})



sapply(1:nrow(q), function(i){ length(wut[[i]])}) == sapply(1:nrow(q), function(i){ length(wut2[[i]])})
*/
