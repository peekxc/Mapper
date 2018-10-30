// skeleton.cpp
// Primary functions for building the k-skeletons
// Includes exported functions for building the skeletons both with and without the simplex tree. 
#include "skeleton.h"

std::vector< IntegerVector > apply_clustering(const NumericMatrix& X, const IntegerVector& level_set, const Function f){
  // Apply the clustering
  const IntegerVector cl_results = f(_["X"] = X, _["idx"] = level_set);
  const IntegerVector cl_idx = self_match(cl_results) - 1;
  const IntegerVector ids = unique(cl_idx);
  
  // Collect the points into vertices  
  std::vector< IntegerVector > new_vertices(ids.size()); 
  IntegerVector::const_iterator c_i = cl_idx.begin();
  for (int i = 0; c_i != cl_idx.end(); ++c_i, ++i){
    // Rprintf("vertex 0-based index c_i = %d getting pt %d from level set idx = %d\n", *c_i, level_set.at(i), i);
    new_vertices.at(*c_i).push_back(level_set.at(i));
  }

  // Return the newly clustered vertices 
  return(new_vertices);
}


// Updates a single level set.
void update_level_set(const int ls_flat_index, const IntegerVector level_set, const NumericMatrix& X, const Function f, std::map< int, IntegerVector >& v_map, List& ls_vertex_map, SEXP stree){
  Rcpp::XPtr<SimplexTree> stree_ptr(stree); // Collect the simplex tree
    
  // Get which vertices need to be updated 
  IntegerVector c_vertices = as<IntegerVector>(ls_vertex_map.at(ls_flat_index));
  // Rcout << "Current vertices to be removed (0-based): " << c_vertices << std::endl; 
    
  // Remove any vertices associated with the current level set
  // NOTE: The ls_vertex_map maps vertex *ids*, not indices! 
  std::for_each(c_vertices.begin(), c_vertices.end(), [&v_map](const int v_id){
    v_map.erase(std::map<int, IntegerVector>::key_type(v_id));
    // Rprintf("Removing vertex id=%d\n", v_id);
  });
  stree_ptr->remove_vertices(c_vertices); // remove them from the simplex tree
  ls_vertex_map.at(ls_flat_index) = IntegerVector::create(); // Remove from the level set map
    
  if (level_set.size() > 0){
    // Apply the clustering to get the new vertices
    std::vector< IntegerVector > new_vertices = apply_clustering(X, level_set, f);
    const std::size_t n_new_v = new_vertices.size(); 
      
    // Update the map with new vertex map with the new indices
    IntegerVector new_0_simplexes = stree_ptr->vertex_available(n_new_v);
    // Rcout << "Adding new simplexes: " << new_0_simplexes << " to map for ls: " << ls_flat_index << std::endl; 
    ls_vertex_map.at(ls_flat_index) = new_0_simplexes;
    
    // Update the simplex tree by inserting the 0-simplices
    std::for_each(new_0_simplexes.begin(), new_0_simplexes.end(), [&stree_ptr](const int v_i){
      std::vector<unsigned int> v = { static_cast<unsigned int>(v_i) };
      stree_ptr->insert_simplex(v);
    });
    
    // Replace the old vertices with the new ones
    std::size_t i = 0; 
    std::for_each(new_0_simplexes.begin(), new_0_simplexes.end(), [&i, &ls_flat_index, &v_map, &new_vertices](const int v_i){
      new_vertices.at(i).attr("level_set") = ls_flat_index+1; // 1-based
      // Rprintf("Emplacing new vertex %d w/ pts: ", v_i);
      // Rcout << new_vertices.at(i) << std::endl; 
      v_map[v_i] = new_vertices.at(i++);
    }); 
  }
} 

// Given a vector of level sets to check, retrieves their corresponding vertices, and using those vertices 
// finds adjacency relations. Return the level sets of these adjacent vertices. 
// [[Rcpp::export]]
IntegerVector check_connected(const IntegerVector ls_to_check, const List& ls_vertex_map, const List& vertices, SEXP stree){
  Rcpp::XPtr< SimplexTree > stree_ptr(stree);
  IntegerVector res = IntegerVector(); 
    
  // For each level set to check, get the vertices, and for each of those vertices, get it's adjacent vertices, 
  // and for each of those vertices, get their parent level set...
  std::for_each(ls_to_check.begin(), ls_to_check.end(), [&](const int c_ls){
    const IntegerVector c_vertices = ls_vertex_map.at(c_ls); // vertex ids
    std::for_each(c_vertices.begin(), c_vertices.end(), [&](const int v_i){
      const IntegerVector adj_vertices = stree_ptr->adjacent_vertices(v_i);
      // Rcout << "Checking adjacent vertices: " << adj_vertices << std::endl; 
      std::for_each(adj_vertices.begin(), adj_vertices.end(), [&]( const int adj_v ){
        const IntegerVector& tmp = vertices[ std::to_string(adj_v) ];
        res.push_back(tmp.attr("level_set")); // pushes back 1-based level set  
      });
    });
  });
  
  // Return
  return(unique(res));
}

// Turns vertex list into a mapping between vertex id --> points 
std::map< int, IntegerVector > vertices_to_map(List& ls_vertex_map, List& vertices, SEXP stree){
  Rcpp::XPtr<SimplexTree> stree_ptr(stree); // Extract the simplex tree
  const std::size_t n_level_sets = ls_vertex_map.size();
  std::map< int, IntegerVector > res; 
  
  // Loop through each level set, creating the mapping along the way
  for (std::size_t i = 0; i < n_level_sets; ++i){
    const IntegerVector v_ids = ls_vertex_map.at(i);
    
    // Save the vertices into the map with their corresponding ids. 
    std::for_each(v_ids.begin(), v_ids.end(), [&stree_ptr, &vertices, &res](const int v_id){
      const std::size_t v_idx = stree_ptr->find_vertex(v_id);
      const IntegerVector& pts = as<IntegerVector>(vertices.at(v_idx)); // okay to access by index
      res.emplace(v_id, pts);
    });
  }
  return(res); 
}

// Builds the 0-skeleton by applying a given clustering function 'f' to a given set of level sets indexed by 'lsfis'.
// Only the level sets in 'which_levels' are updated. Returns a new list containing the updated vertices.
// Parameters:
//  which_levels := An integer vector of 0-based indices whose level sets are to be updated
//  f := The clustering function, which takes in a data matrix 'X' and an integer vector of point indices 'idx' to cluster on.
//  X := The data matrix. 
//  level_sets := List of integer vectors containing the level sets =
//  vertices := List of the current vertices (each element of which is a vector containing the point indices contained in the vertex)
//  ls_vertex_map := Ordered list of size (n_level_sets) where each element comprises the vertex ids contained in the level set at that index
//  stree := SimplexTree object
// [[Rcpp::export]]
List build_0_skeleton(const IntegerVector which_levels, const NumericMatrix& X, Function f, const List& level_sets, List& vertices, List& ls_vertex_map, SEXP stree){
  
  // std::vector< IntegerVector > res(vertices.begin(), vertices.end());
  std::map< int, IntegerVector > v_map = vertices_to_map(ls_vertex_map, vertices, stree);
  
  // Loop through the levels to update. This will update the vertices, the vertex map, and the simplex tree
  for (const int index: which_levels){
    // Rcout << "UPDATING LEVEL SET: " << index << std::endl;
    const IntegerVector ls_pt_idx = level_sets.at(index);
    update_level_set(index, ls_pt_idx, X, f, v_map, ls_vertex_map, stree);
  }
  
  // Return result wrapped as a List
  return(wrap(v_map));
  // return(List::create());
}



// Builds the 1-skeleton by inserting 1-simplexes for each pair of vertices whose points have non-empty 
// intersections. Will only compare vertices given by the 'ls_pairs' matrix.
//  ls_pairs := (n x 2) Integer Matrix of n level set index pairs (by flat index) to consider (0-based)
//  min_weight := integer representing the minimum intersection size to make an edge.
//  vertices := List of nodes (each element of which is a vector containing the point indices contained in the node)
//  ls_vertex_map := List where each index corresponds to the ordered level set flat indices, and each element the indices of the nodes in that level set
//  stree := SimplexTree object
// [[Rcpp::export]]
void build_1_skeleton(const IntegerMatrix& ls_pairs, const int min_sz, const List& vertices, const List& ls_vertex_map, SEXP stree){
  Rcpp::XPtr<SimplexTree> stree_ptr(stree); // Collect the simplex tree
  const bool check_sz = min_sz > 1;
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
        const IntegerVector& n1_idx = vertices[ std::to_string(*n1) ]; // access by vertex id
        const IntegerVector& n2_idx = vertices[ std::to_string(*n2) ]; // access by vertex id
        
        if (check_sz){
          // Add edge between the two if they share a data point. This also retrieves the size of the intersection.
          int intersect_size = std::count_if(n1_idx.begin(), n1_idx.end(), [&](int k) {
            return(std::find(n2_idx.begin(), n2_idx.end(), k) != n2_idx.end());
          });
          if (intersect_size >= min_sz){
            std::vector<uint> simplex = { uint(*n1), uint(*n2) };
            stree_ptr->insert_simplex(simplex);
          }
        } else {
          // Add edge between the two if they have a non-empty intersection.
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


// Multiscale-version of updating the vertices in the level sets. Unlike the regular version, 
// the level sets are not explicitly stored in memory, but are rather extracted from the segments dynamically. 
// [[Rcpp::export]]
List update_level_sets(const IntegerVector which_levels, SEXP ms, const NumericMatrix& X, const Function f, List& vertices, List& ls_vertex_map, SEXP stree){
  Rcpp::XPtr< MultiScale > ms_ptr(ms); // Cast the pointer to get the multiscale structure
  
  // // Shallow(?) copy the internals vectors 
  // std::vector< IntegerVector > res(vertices.begin(), vertices.end());
  std::map< int, IntegerVector > v_map = vertices_to_map(ls_vertex_map, vertices, stree);

  // Update the vertices in the level sets dynamically
  for (const int index: which_levels){
    const IntegerVector level_set = ms_ptr->extract_level_set(index) + 1; // convert to 1-based for R
    update_level_set(index, level_set, X, f, v_map, ls_vertex_map, stree);
  }
  
  // Return the new vertices
  return(wrap(v_map));
}

