#include <Rcpp.h>
using namespace Rcpp;

#include "SimplexTree.h"
#include "utility_rcpp.h"

// Fast partial-sort/binary-search check to see if the intersection between two given vectors has non-zero length
// Loosely based off of bugged version found here: https://stackoverflow.com/questions/21359432/a-c-version-of-the-in-operator-in-r
bool any_is_in(const IntegerVector& x, const IntegerVector& y){
  std::vector<int> y_sort(y.size());
  std::partial_sort_copy (y.begin(), y.end(), y_sort.begin(), y_sort.end()); // partial-sorted elements of y copied to y_sort
  const int nx = x.size();
  for (int i = 0; i < nx; ++i) {
    if (std::binary_search(y_sort.begin(), y_sort.end(), x[i])) {
      return(true); // end the search
    }
  }
  return(false);
}

// rbindlist_int: Takes a list of integer vectors and rbind's them together.
IntegerMatrix rbindlist_int(std::list<IntegerVector>& lst){
  unsigned int n = lst.size();
  if(n == 0) { Rcpp::stop("Invalid sized list."); }
  unsigned int d = lst.front().size();
  Rcpp::IntegerMatrix res = Rcpp::no_init(n, d);
  size_t i = 0;
  for (std::list<IntegerVector>::iterator it = lst.begin(); it != lst.end(); ++it, ++i) {
    if ((*it).size() != d) { Rcpp::stop("Invalid integer vector size."); }
    res(i, _) = *it;
  }
  return res;
}

// edgeList_int: Creates an edgelist forming edges between nodes w/ non-empty intersections
// Parameters:
//  ls_pairs := (n x 2) Integer Matrix of n level index pairs to consider
//  nodes := List of nodes (each element of which is a vector containing the point indices contained in the node)
//  node_map := List where each index corresponds to the ordered level set flat indices, and each element the indices of the nodes in that level set
// [[Rcpp::export]]
IntegerMatrix edgeList_int(const IntegerMatrix& ls_pairs, const List& nodes, const List& ls_node_map) {
  int n = nodes.length();
  std::list<IntegerVector> edgelist = std::list<IntegerVector>();
  for (int i = 0; i < ls_pairs.nrow(); ++i){

    // Get the current pair of level sets to compare; skip if either are empty
    const int ls_1 = ls_pairs(i, 0), ls_2 = ls_pairs(i, 1);
    if ( Rf_isNull(ls_node_map.at(ls_1 - 1)) || Rf_isNull(ls_node_map.at(ls_2 - 1))){
      continue;
    }
    const IntegerVector& nodes1 = ls_node_map.at(ls_1 - 1);
    const IntegerVector& nodes2 = ls_node_map.at(ls_2 - 1);

    // Compare the nodes within each level set
    for (IntegerVector::const_iterator n1 = nodes1.begin(); n1 != nodes1.end(); ++n1){
      for (IntegerVector::const_iterator n2 = nodes2.begin(); n2 != nodes2.end(); ++n2){
        // Retrieve point indices within each node
        const IntegerVector& n1_idx = nodes[*n1 - 1];
        const IntegerVector& n2_idx = nodes[*n2 - 1];
        // Add edge between the two if they share a data point
        bool intersect_check = any_is_in(n1_idx, n2_idx);
        if (intersect_check){
          edgelist.push_back(IntegerVector::create(*n1 - 1, *n2 - 1));
        }
      }
    }
  }
  return rbindlist_int(edgelist);
}


// adjacencyCpp: Creates an adjacency matrix forming edges between nodes w/ non-empty intersections
// Parameters:
//  ls_pairs := Integer Matrix of level index pairs to consider
//  nodes := List of nodes (each element of which is a vector containing the point indices contained in the node)
//  node_map := List where each index corresponds to the ordered level set flat indices, and each element the indices of the nodes in that level set
// [[Rcpp::export]]
IntegerMatrix adjacencyCpp(const IntegerMatrix& ls_pairs, const List& nodes, const List& ls_node_map) {
  int n = nodes.length();
  IntegerMatrix adj_mat(n, n);
  for (int i = 0; i < ls_pairs.nrow(); ++i){

    // Get the current pair of level sets to compare; skip if either are empty
    const int ls_1 = ls_pairs(i, 0), ls_2 = ls_pairs(i, 1);
    if ( Rf_isNull(ls_node_map.at(ls_1 - 1)) || Rf_isNull(ls_node_map.at(ls_2 - 1))){
      continue;
    }
    const IntegerVector& nodes1 = ls_node_map.at(ls_1 - 1);
    const IntegerVector& nodes2 = ls_node_map.at(ls_2 - 1);

    // Compare the nodes within each level set
    for (IntegerVector::const_iterator n1 = nodes1.begin(); n1 != nodes1.end(); ++n1){
      for (IntegerVector::const_iterator n2 = nodes2.begin(); n2 != nodes2.end(); ++n2){
        // Retrieve node indices
        const IntegerVector& n1_idx = nodes[*n1 - 1];
        const IntegerVector& n2_idx = nodes[*n2 - 1];
        // Add edge between the two if they share a data point
        bool intersect_check = any_is_in(n1_idx, n2_idx);
        if (intersect_check){
          adj_mat(*n1 - 1, *n2 - 1) = 1;
          adj_mat(*n2 - 1, *n1 - 1) = 1;
        }
      }
    }
  }
  return adj_mat;
}

// Builds the 0-skeleton by applying a given clustering function 'f' to a given set of level sets indexed by 'lsfis'.
// Only the level sets in 'which_levels' are updated. Returns a new list containing the updated vertices.
// Parameters:
//  which_levels := An integer vector of 0-based indices whose level sets are to be updated
//  f := The clustering function, which takes in a data matrix 'X' and an integer vector of point indices 'idx' to cluster on.
//  X := The data matrix. 
//  level_sets := List of integer vectors containing the level sets =
//  vertices := List of the current vertices (each element of which is a vector containing the point indices contained in the vertex)
//  ls_vertex_map := Ordered list of size (n_level_sets) where each element the indices of the vertices in the level set at that index
//  stree := SimplexTree object
// [[Rcpp::export]]
List build_0_skeleton(const IntegerVector which_levels, const NumericMatrix& X, Function f, const List& level_sets, List& vertices, List& ls_vertex_map, SEXP stree){
  Rcpp::XPtr<SimplexTree> stree_ptr(stree); // Collect the simplex tree
  int n = vertices.length();
  
  // Shallow(?) copy the internals vectors 
  std::vector< IntegerVector > res(vertices.begin(), vertices.end());
  //return(wrap(res));
  
  // Loop through the levels to update 
  for (const int index: which_levels){
    
    Rcout << "index: " << index << std::endl; 
    
    // Get which vertices need to be updated 
    IntegerVector c_vertices = as<IntegerVector>(ls_vertex_map.at(index));
    Rcout << "Current vertices to be removed: " << c_vertices << std::endl; 
    
    // Remove any vertices associated with the current level set
    stree_ptr->remove_vertices(c_vertices); // remove in the simplex tree
    std::for_each(c_vertices.begin(), c_vertices.end(), [&res](const int v_idx){
      Rcout << "Removing: " << v_idx << std::endl; 
      res.at(v_idx) = NULL;
    });
    
    // Get the point indices we're working with
    const IntegerVector ls_pt_idx = level_sets.at(index);
    Rcout << "point indices to cluster: " <<  ls_pt_idx << std::endl; 
    
    if (ls_pt_idx.size() > 0){
      
      // Apply the clustering
      const IntegerVector cl_results = f(_["X"] = X, _["idx"] = ls_pt_idx);
      const IntegerVector cl_idx = self_match(cl_results) - 1;
      const IntegerVector ids = unique(cl_idx);
      
      Rcout << "Clustering results: " << cl_results << std::endl; 
      Rcout << "Cluster ids: " << ids << std::endl; 
      Rcout << "id size: " << ids.size() << std::endl; 
      
      // Update the map with new vertex map with the new indices
      IntegerVector new_0_simplexes = stree_ptr->vertex_available(ids.size());
      Rcout << "Adding new simplexes: " << new_0_simplexes << " to map at " << index << std::endl; 
      Rcout << new_0_simplexes.size() << std::endl;  
      ls_vertex_map.at(index) = new_0_simplexes;
      
      // Update the simplex tree by inserting the 0-simplices
      std::for_each(new_0_simplexes.begin(), new_0_simplexes.end(), [&stree_ptr](const int v_i){
        std::vector<unsigned int> v = { static_cast<unsigned int>(v_i) };
        stree_ptr->insert_simplex(v);
      });
      // Rcout << "New vertex ids: " << as<IntegerVector>(ls_vertex_map.at(index)) << std::endl; 
      
      // Create the new vertices
      Rprintf("creating %d new vertices\n", ids.size()); 
      std::vector< IntegerVector > new_vertices(ids.size()); 
      IntegerVector::const_iterator c_i = cl_idx.begin();
      for (int i = 0; c_i != cl_idx.end(); ++c_i, ++i){
        Rprintf("c_i = %d, i = %d\n", *c_i, i);
        Rcout << ls_pt_idx << std::endl; 
        new_vertices.at(*c_i).push_back(ls_pt_idx.at(i));
      }
      
      // Replace the old vertices with the new ones
      std::size_t i = 0; 
      Rcout << "Replacing old vertices" << std::endl; 
      std::for_each(new_0_simplexes.begin(), new_0_simplexes.end(), [&i, &index, &res, &new_vertices](const int v_i){
        if (res.size() < v_i+1){ res.resize(v_i+1); } 
        new_vertices.at(i).attr("level_set") = index+1; // 1-based
        res.at(v_i) = new_vertices.at(i++);
      }); 
    }
  }
  
  return(wrap(res));
}

// Builds the 1-skeleton by inserting 1-simplexes for each pair of vertices whose points have non-empty 
// intersections. Will only compare vertices given by the 'ls_pairs' matrix.
//  ls_pairs := (n x 2) Integer Matrix of n level set index pairs (by flat index) to consider
//  vertices := List of nodes (each element of which is a vector containing the point indices contained in the node)
//  ls_vertex_map := List where each index corresponds to the ordered level set flat indices, and each element the indices of the nodes in that level set
//  stree := SimplexTree object
// [[Rcpp::export]]
void build_1_skeleton(const IntegerMatrix& ls_pairs, const List& vertices, const List& ls_vertex_map, SEXP stree){
  Rcpp::XPtr<SimplexTree> stree_ptr(stree); // Collect the simplex tree
  int n = vertices.length();
  for (int i = 0; i < ls_pairs.nrow(); ++i){
    
    // Get the current pair of level sets to compare; skip if either are empty
    const int ls_1 = ls_pairs(i, 0), ls_2 = ls_pairs(i, 1);
    if ( Rf_isNull(ls_vertex_map.at(ls_1)) || Rf_isNull(ls_vertex_map.at(ls_2))){
      continue;
    }
    const IntegerVector& nodes1 = ls_vertex_map.at(ls_1);
    const IntegerVector& nodes2 = ls_vertex_map.at(ls_2);
    
    // Compare the nodes within each level set
    for (IntegerVector::const_iterator n1 = nodes1.begin(); n1 != nodes1.end(); ++n1){
      for (IntegerVector::const_iterator n2 = nodes2.begin(); n2 != nodes2.end(); ++n2){
        // Retrieve point indices within each node
        const IntegerVector& n1_idx = vertices.at(*n1);
        const IntegerVector& n2_idx = vertices.at(*n2);
        
        // Add edge between the two if they share a data point. This also retrieves the size of the intersection.
        // int intersect_size = std::count_if(n1_idx.begin(), n1_idx.end(), [&](int k) { 
        //   return(std::find(n2_idx.begin(), n2_idx.end(), k) != n2_idx.end());
        // });
        bool intersect_check = std::any_of(n1_idx.begin(), n1_idx.end(), [&](int k){
          return(std::find(n2_idx.begin(), n2_idx.end(), k) != n2_idx.end());
        });
        if (intersect_check){
          std::vector<uint> simplex = { uint(*n1), uint(*n2) };
          stree_ptr->insert_simplex(simplex);
        }
      }
    }
  }
  // return rbindlist_int(edgelist);
}



// Returns a graph representing the intersection between two lists of nodes.
// Given two lists of nodes, checks each pairwise combination of nodes to see if they intersect.
// If they contain a non-empty intersection, an edge is created.
// [[Rcpp::export]]
List intersectNodes(const List& nodes1, const List& nodes2, const IntegerVector& node_ids1, const IntegerVector& node_ids2) {

  // Compres two groups of nodes for non-empty intersection
  std::list<IntegerVector> edgelist = std::list<IntegerVector>();
  int i = 0, j = 0;
  for (List::const_iterator n1 = nodes1.begin(); n1 != nodes1.end(); ++n1, ++i){
    j = 0;
    for (List::const_iterator n2 = nodes2.begin(); n2 != nodes2.end(); ++n2, ++j){
      // Retrieve point indices within each node
      const IntegerVector& n1_idx = *n1;
      const IntegerVector& n2_idx = *n2;

      // Add edge between the two if they share a data point
      bool intersect_check = any_is_in(n1_idx, n2_idx);
      if (intersect_check){
        edgelist.push_back(IntegerVector::create(node_ids1[i], node_ids2[j]));
      }
    }
  }
  return(wrap(edgelist));
}


/*** R
# load("test.rdata")
# adjacencyCpp(ls_pairs = test$ls_pairs, nodes = test$nodes, ls_node_map = test$ls_node_map)
*/
