// skeleton.cpp
// Primary functions for building the k-skeletons
// Includes exported functions for building the skeletons both with and without the simplex tree. 
#include "skeleton.h"

// Checks for empty intersection between sorted vectors
// bool empty_intersection(const IntegerVector& x, const IntegerVector& y) {
//   IntegerVector::const_iterator i = x.begin(), j = y.begin();
//   while (i != x.end() && j != y.end()){
//     if (*i<*j) ++i; else if(*j<*i) ++j; else return false;
//   }
//   return true;
// }


// Given a set of pullback indices and a simplex tree, for each pullback id, find the ids the
// of the other pullback sets whose corresponding vertices participate in a coface of the vertices 
// in the source pullback. 
// [[Rcpp::export]]
List connected_pullbacks(std::vector< std::string > pullback_ids, const List& pullback, SEXP stree){
  Rcpp::XPtr< SimplexTree > stree_ptr(stree);

  // Get all the pullback ids
  vector< std::string > all_pids = as< vector< std::string > >(pullback.names());
  
  // Reverse the pullback maps (key, value) pairs into a new map
  std::map< size_t, std::set< std::string > > v_to_pb;  
  for (std::string& pid: all_pids){
    const IntegerVector vids = as< IntegerVector >(pullback[pid]);
    for (size_t v: vids){ v_to_pb[v].insert( pid ); }
  }
  
  // Find pullback ids of adjacent vertices
  std::map< std::string, std::set< std::string > > pb_to_adj;
  for (std::string& pid: pullback_ids){
    IntegerVector vids = pullback[pid];
    for (int vid: vids){
      vector< idx_t > v_adj = stree_ptr->adjacent_vertices(vid);
      for (idx_t v: v_adj){
        pb_to_adj[pid].insert(v_to_pb[v].begin(), v_to_pb[v].end());
      }
    }
  }
  return wrap(pb_to_adj);
}

// Given the data set X, the current set of indices level_set, and a clustering function f, 
// apply the clustering function to the subset of X given by the level set, and returns a 
// vector of index vector representing which points fell into which partition. 
vector< IntegerVector > apply_clustering(std::string pid, const Function& level_set_f, const Function& cluster_f){
  const IntegerVector preimages = level_set_f(pid);
  if (preimages.size() == 0){ 
    vector<IntegerVector> res; 
    return(res); 
  } 
  const IntegerVector cl_results = cluster_f(pid, preimages);
  // Rcout << "cluster results: " << cl_results << std::endl;
  const IntegerVector cl_idx = self_match(cl_results) - 1; // guarenteed to be 0-based contiguous indices
  // Rcout << "self match results: " << cl_idx << std::endl;
  const IntegerVector ids = unique(cl_idx); 
  // Rcout << "cluster ids: " << ids << std::endl;
  
  // Collect the points into vertices  
  size_t i = 0; 
  vector< IntegerVector > new_vertices( ids.size() );
  for (const int ci: cl_idx){
    new_vertices.at(ci).push_back(preimages.at(i++));
  }
  // IntegerVector::const_iterator c_i = cl_idx.begin();
  // for (int i = 0; c_i != cl_idx.end(); ++c_i, ++i){
  //   new_vertices.at(*c_i).push_back(level_set.at(i));
  // }
  
  // Return the newly clustered vertices 
  return(new_vertices);
}

// Given a vector of integer ids, converts them to strings, and removes any vertices which 
// have an id as in the given id list. 
List remove_by_id(IntegerVector ids, const List& vertices){
  if (Rf_isNull(vertices.names()) || vertices.size() == 0){ return vertices; }
  if (ids.size() == 0){ return vertices; }
  const size_t n = ids.size(); 
  vector< std::string > vnames = as< vector< std::string > >(vertices.names());
  LogicalVector subset = LogicalVector(vertices.size(), true);
  for (size_t i = 0; i < n; ++i){
    std::string id_str = std::to_string(ids.at(i));
    auto vn_it = std::find(begin(vnames), end(vnames), id_str);
    if (vn_it != vnames.end()){
      size_t to_remove = std::distance(begin(vnames), vn_it);
      subset.at(to_remove) = false; 
    }
  }
  return vertices[subset];
  // erase_partition(vertices, begin(vnames), end(vnames), [](){
  //   
  // })
  // vector< std::string > vnames = as< vector< std::string > >(vertices.names());
  // vector< size_t > v_idx = seq_ij< size_t >(0, vertices.size()-1);
  // 
  // erase_partition(v_idx, begin(v_idx), end(v_idx), [&vnames, &ids](const size_t vid){
  //   auto id_it = std::find(begin(ids), end(ids), vid);
  //   if ()
  //   return bool(std::find(begin(vnames), end(vnames), std::to_string(vid)) != vnames.end());
  // });
  // IntegerVector tmp_idx = IntegerVector(begin(v_idx), end(v_idx));
  // vertices = vertices[ tmp_idx ]; 
}



// Updates the connected components for the given pullback ids. 
// This involves, for each pullback id 'pid': 
// 1) Removing the vertices (and cofaces) of the components mapped to by pid
// 2) Retrieving the level set for pid. 
// 3) Breaking the level set given by (2) into connected components. 
// 4) Inserting the new vertices into the simplex tree, as well as updating the current vertex list 
// Note that higher order simplices must be rebuilt after this method is called in conjunction with, e.g. connected_pullback.
// void update_pullback(
//     const std::string pid, 
//     const Function level_set_f, 
//     const NumericMatrix& X, 
//     const Function cluster_f, 
//     List& pullback, 
//     List& vertices, 
//     SEXP stree)
// {
//   Rcpp::XPtr<SimplexTree> stree_ptr(stree); // Collect the simplex tree
//   
//   // Get which vertices need to be updated 
//   Rcout << "here1" << std::endl;
//   IntegerVector c_vertices = as<IntegerVector>(pullback[pid]);
//   
//   // Remove existing vertices 
//   Rcout << "here2" << std::endl;
//   if (c_vertices.size() > 0 && !Rf_isNull(vertices.names())){
//     vector< std::string > vnames = as< vector< std::string > >(vertices.names());
//     vector< size_t > v_idx = seq_ij< size_t >(0, vertices.size()-1);
//     erase_partition(v_idx, begin(v_idx), end(v_idx), [&vnames](const size_t vid){
//       c_vertices
//       return bool(std::find(begin(vnames), end(vnames), std::to_string(vid)) != vnames.end());
//     });
//     IntegerVector tmp_idx = IntegerVector(begin(v_idx), end(v_idx));
//     vertices = vertices[ tmp_idx ]; 
//   }
//   
//   // Clear the pullback (NOTE: The pullback maps vertex *ids*, not indices!)
//   Rcout << "here3" << std::endl;
//   pullback[pid] = IntegerVector::create();
//   
//   // Remove the vertices from the simplex tree and the vertex list
//   for (int vid: c_vertices){ stree_ptr->remove_simplex({ size_t(vid) }); };
//     
//   // Construct new components
//   Rcout << "here4" << std::endl;
//   const IntegerVector level_set = level_set_f(pid);
//   if (level_set.size() > 0){
//     // Apply the clustering to get the new vertices
//     vector< IntegerVector > new_vertices = apply_clustering(X, level_set, cluster_f);
//       
//     // Update both the simplex tree and the pullback map with the new vertex ids
//     vector< idx_t > new_0_simplexes = stree_ptr->generate_ids(new_vertices.size());
//     pullback[pid] = wrap(new_0_simplexes);
//     for (idx_t v: new_0_simplexes){ stree_ptr->insert_simplex({ v }); }
//     
//     // Insert point indices into vertices
//     for (size_t i = 0; i < new_vertices.size(); ++i){
//       Rcout << "Inserting: " << std::to_string(new_0_simplexes.at(i)) << std::endl;
//       Rcout << new_vertices.at(i) << std::endl;
//       new_vertices.at(i).attr("level_set") = pid; 
//       vertices[std::to_string(new_0_simplexes.at(i))] = new_vertices.at(i);
//     }
//   } else { pullback[pid] = IntegerVector::create(); }
// } 

// Given a vector of level sets to check, retrieves their corresponding vertices, and using those vertices 
// finds adjacency relations. Return the level sets of these adjacent vertices. 
// IntegerVector check_connected(const IntegerVector ls_to_check, const List& ls_vertex_map, const List& vertices, SEXP stree){
//   Rcpp::XPtr< SimplexTree > stree_ptr(stree);
//   IntegerVector res = IntegerVector();
// 
//   // For each level set to check, get the vertices, and for each of those vertices, get it's adjacent vertices,
//   // and for each of those vertices, get their parent level set...
//   std::for_each(ls_to_check.begin(), ls_to_check.end(), [&](const int c_ls){
//     const IntegerVector c_vertices = ls_vertex_map.at(c_ls); // vertex ids
//     std::for_each(c_vertices.begin(), c_vertices.end(), [&](const int v_i){
//       const IntegerVector adj_vertices = wrap(stree_ptr->adjacent_vertices(v_i));
//       // Rcout << "Checking adjacent vertices: " << adj_vertices << std::endl;
//       std::for_each(adj_vertices.begin(), adj_vertices.end(), [&]( const int adj_v ){
//         const IntegerVector& tmp = vertices[ std::to_string(adj_v) ];
//         res.push_back(tmp.attr("level_set")); // pushes back 1-based level set
//       });
//     });
//   });
// 
//   // Return
//   return(unique(res));
// }

// Converts the List pullback to a C++ map
// std::map< int, IntegerVector > vertices_to_map(List& pullback, List& vertices){
//   const size_t n_pullbacks = pullback.size();
//   std::map< int, IntegerVector > pullback_map; 
//   
//   // Loop through each level set, creating the mapping along the way
//   for (size_t i = 0; i < n_pullbacks; ++i){
//     const IntegerVector v_ids = pullback.at(i);
//     
//     // Save the vertices into the map with their corresponding ids. 
//     std::for_each(v_ids.begin(), v_ids.end(), [&vertices, &pullback_map](const int v_id){
//       std::string v_key = std::to_string(v_id);
//       const IntegerVector& pts = as<IntegerVector>(vertices[v_key]); // okay to access by index
//       pullback_map.emplace(v_id, pts);
//     });
//   }
//   return(pullback_map); 
// }

vector< size_t > smallest_not_in(const size_t n){
  vector< size_t > new_ids = vector< size_t >();
  new_ids.resize(n);
  std::iota(begin(new_ids), end(new_ids), 0);
  return(new_ids);
}

vector< size_t > smallest_not_in(const size_t n, vector< size_t > old_ids){
  // if (any(old_ids < 0).is_true()){ stop("Supplied ids must be non-negative integers."); }
  const size_t max = old_ids.size() + n;
  std::sort(begin(old_ids), end(old_ids));
  vector< size_t > new_ids = vector< size_t >();
  new_ids.resize(n);
  for (size_t i = 0, cc = 0; i < max; ++i){
    auto lb = std::find(begin(old_ids), end(old_ids), i); 
    if (lb == end(old_ids)){ new_ids[cc++] = i; } // i not in old_vids
  }
  return(new_ids);
}

std::unordered_map< std::string, IntegerVector > lst2map(const List& lst){
  if (Rf_isNull(lst.names()) || lst.size() == 0){ return std::unordered_map< std::string, IntegerVector >(); }
  vector< std::string > keys = as< vector< std::string > >(lst.names());
  std::unordered_map< std::string, IntegerVector > res; 
  for (auto c_key: keys){ res.emplace(c_key, lst[c_key]); }
  return(res);
}

// Decomposes the preimages into connected components according to a given clustering function.
// This modifies both the vertex list and the pullback mapping. 
// This method assumes the names in the 'pullback' list are constant! 
// Parameters:
//  pullback_ids := The pullback ids to update by this method. Only performs the pullback operation on these keys.
//  X := The data matrix. 
//  f := The clustering function, which takes in a data matrix 'X' and an integer vector of point indices 'idx' to cluster on.
//  level_set_f := Function which takes as input a pullback id and returns as output the indices of the point within that pullback set. 
//  vertices := List of the current vertices (each element of which is a vector containing the point indices contained in the vertex)
//  pullback := Named-indexed list mapping pullback ids to vertex ids
// [[Rcpp::export]]
List decompose_preimages(
    const StringVector pullback_ids, 
    const Function cluster_f, 
    const Function level_set_f, 
    const List& vertices, 
    List& pullback)
{
  using str_map = std::unordered_map< std::string, IntegerVector >; 
  
  // Extract the current set of vertices in C++ 
  str_map mod_vertices = lst2map(vertices);
  
  // Create lambda to extract the vertex ids as integers
  const auto vertex_ids = [&mod_vertices](){ 
    using v_type = std::unordered_map< std::string, IntegerVector >::value_type;
    if (mod_vertices.empty()){ return(vector< size_t >()); }
    vector< size_t > vids(mod_vertices.size());
    std::transform(begin(mod_vertices), end(mod_vertices), begin(vids), [](v_type v){
      return(std::stoi(v.first));
    });
    return(vids);
  };
  
  // Loop through the pullback ids to update. This will update the vertices, the pullback map, and the simplex tree
  vector< std::string > pids = as< vector< std::string > >(pullback_ids);
  for (std::string pid: pids){
    IntegerVector vids = pullback[pid];
    
    // Remove vids in vertex list, pullback, and simplex tree
    // Rcout << "Removings vertices: " << vids << std::endl;
    // vertices = remove_by_id(vids, vertices);
    for (auto id: vids){ mod_vertices.erase(std::to_string(id)); }

    pullback[pid] = IntegerVector::create();
    
    // Apply the clustering to get the new vertices
    vector< IntegerVector > new_vertices = apply_clustering(pid, level_set_f, cluster_f);
    
    if (new_vertices.size() == 0){
      pullback[pid] = IntegerVector::create();
    } else {
      // Update both the simplex tree and the pullback map with the new vertex ids
      // vector< idx_t > new_0_simplexes = stree_ptr->generate_ids(new_vertices.size());
      // vector< idx_t > new_0_simplexes = as< vector< idx_t > >(id_generator(new_vertices.size()));
      vector< size_t > c_vids = vertex_ids();
      vector< size_t > new_0_simplexes = smallest_not_in(new_vertices.size(), c_vids);
      pullback[pid] = wrap(new_0_simplexes);
      
      // Insert point indices into vertices. New vertex ids should be guarenteed to not be in vertex ids
      for (size_t i = 0; i < new_vertices.size(); ++i){
        new_vertices.at(i).attr("level_set") = pid; 
        std::string new_vid = std::to_string(new_0_simplexes.at(i));
        mod_vertices.emplace(new_vid, new_vertices.at(i));
      }
    }
  }
  return(wrap(mod_vertices));
}


//  vids := New vertex ids
//  st := SimplexTree object containing the underlying simplicial complex
// [[Rcpp::export]]
void build_0_skeleton(const IntegerVector vids, SEXP st){
  Rcpp::XPtr<SimplexTree> st_ptr(st); // Collect the simplex tree
  st_ptr->clear();
  //for (int vid: vids){ stree_ptr->remove_simplex({ size_t(vid) }); };
  for (idx_t v: vids){ st_ptr->insert_simplex({ v }); }
}

// Builds the 1-skeleton by inserting 1-simplexes for each pair of vertices whose points have non-empty 
// intersections. Will only compare vertices given by the 'ls_pairs' matrix.
//  ls_pairs := (n x 2) Integer Matrix of n level set index pairs (by flat index) to consider (0-based)
//  min_weight := integer representing the minimum intersection size to make an edge.
//  vertices := List of nodes (each element of which is a vector containing the point indices contained in the node)
//  ls_vertex_map := List where each index corresponds to the ordered level set flat indices, and each element the indices of the nodes in that level set
//  stree := SimplexTree object
// [[Rcpp::export]]
List build_1_skeleton(const CharacterMatrix& pullback_ids, const int min_sz, const List& vertices, const List& pullback, SEXP stree, const bool modify){
  if (pullback_ids.ncol() != 2){ stop("Expecting pullback ids to be (n x 2) matrix."); }
  Rcpp::XPtr<SimplexTree> stree_ptr(stree); // Collect the simplex tree
  const bool check_sz = min_sz > 1;
  
  // Track the simplices added to the complex
  vector< vector< idx_t > > simplices_added; 
  for (int i = 0; i < pullback_ids.nrow(); ++i){
    
    // Get the current pair of level sets to compare; skip if either are empty
    const std::string pid1 = as< std::string >(pullback_ids(i, 0)); 
    const std::string pid2 = as< std::string >(pullback_ids(i, 1));
    if ( Rf_isNull(pullback[pid1]) || Rf_isNull(pullback[pid2])){
      continue;
    }
    const IntegerVector& nodes1 = pullback[pid1];
    const IntegerVector& nodes2 = pullback[pid2];
    
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
            vector< size_t > simplex = { size_t(*n1), size_t(*n2) };
            if (modify){ stree_ptr->insert_simplex(simplex); } 
            else { simplices_added.push_back(simplex); }
          } 
        } else {
          // Add edge between the two if they have a non-empty intersection.
          bool intersect_check = std::any_of(n1_idx.begin(), n1_idx.end(), [&](int k){
            return(std::find(n2_idx.begin(), n2_idx.end(), k) != n2_idx.end());
          });
          if (intersect_check){
            vector<size_t> simplex = { size_t(*n1), size_t(*n2) };
            if (modify){ stree_ptr->insert_simplex(simplex); }
            else { simplices_added.push_back(simplex); }
          } 
        }
      }
    }
  }
  return(wrap(simplices_added));
}

// Builds the k-skeleton by inserting k-simplexes which have a nonempty intersection 
// between the pairs of vertices given by combinations of each row. 
// Assumes the faces of the k-simplexes to be inserted 
// have already been inserted via, e.g. construction of the (k-1)-skeleton. 
//  k := the dimension to compute. Must be above 1. 
//  stree := SimplexTree object
// [[Rcpp::export]]
List build_k_skeleton(CharacterMatrix pullback_ids, const List& pullback, List& vertices, int k, SEXP stree, const bool modify){
  if (k < 2){ stop("'build_k_skeleton' is meant for building K-complexes for k > 1."); }
  if (pullback_ids.ncol() != k+1){ stop("Expecting (n x k+1) matrix fo pullback ids to compare."); }
  Rcpp::XPtr<SimplexTree> stree_ptr(stree); // Collect the simplex tree
  using uint = unsigned int;

  // If requested, track which simplices were added
  vector< vector< idx_t > > simplices_added; 
  
  // Iterate through the pullbacks 
  const size_t n = pullback_ids.nrow();
  for (size_t i = 0; i < n; ++i){
    if ((i % 100) == 0){ Rcpp::checkUserInterrupt(); }
    CharacterVector tmp = pullback_ids(i, _);
    vector< std::string > pids = as< vector< std::string > >(tmp);
    
    // If any of the pullback covers are empty, skip
    bool any_empty = std::accumulate(begin(pids), end(pids), false, [&pullback](bool any_empty, const std::string pid) -> bool{
      if (any_empty || bool(Rf_isNull(pullback[pid]))){ return true; }
      switch(TYPEOF(pullback[pid])){
      case INTSXP: return(as< IntegerVector >(pullback[pid]).size() == 0);
      case REALSXP: return(as< NumericVector >(pullback[pid]).size() == 0);
      default: stop("SEXP type of pullback value not handled."); 
      }
    });
    if (any_empty){ continue; }

    // Get a vector of vertices per pullback id
    // if (pids.size() != k+1){ stop("pullback ids size != k+1"); }
    vector< vector< size_t > > sigma = vector< vector< size_t > >(k+1);
    std::transform(begin(pids), end(pids), begin(sigma), [&pullback](const std::string pid){
      vector< size_t > p_vertices = as< vector< size_t > >(pullback[pid]);
      return p_vertices;
    });
    
    // Compare every (k+1)-length simplex
    CartesianProduct(sigma, [&vertices, &stree_ptr, &k, &simplices_added, &modify](const vector< size_t > k_simplex){

      // Check that every combination of (k-1) simplices already exist. If this 
      // doesn't pass, the current k-simplex cannot be in the complex. 
      bool all_exists = apply_combinations(k_simplex.size(), k, [&k_simplex, &stree_ptr](const vector< size_t > idx){
        vector< idx_t > simplex = vector< idx_t >(idx.size());
        std::transform(begin(idx), end(idx), begin(simplex), [&k_simplex](const size_t i){ 
          return k_simplex.at(i); 
        });
        bool exists = stree_ptr->find_simplex(simplex);
        return(exists);
      });

      // If all the (k-1) faces exist, check for nonempty intersection
      if (all_exists){
        using Iter = IntegerVector::const_iterator; 
        
        // Const range pair iterators through the points in each vertex
        vector< std::pair<Iter, Iter> > v_rngs(k_simplex.size()); 
        std::transform(begin(k_simplex), end(k_simplex), begin(v_rngs), [&vertices](const size_t vid){
          const IntegerVector vertex = vertices[ std::to_string(vid) ];
          return std::make_pair(vertex.begin(), vertex.end());
        });
        
        // Add a simplex if there's a nonempty intersection between all vertices
        bool valid = nonempty_intersection(v_rngs);
        if (valid){
          if (modify){ stree_ptr->insert_simplex(k_simplex); }
          else { simplices_added.push_back(k_simplex); }
        }
      }
    });
  }
  return(wrap(simplices_added));
}

    // const IntegerVector& nodes1 = pullback.at(ls_1);
    // const IntegerVector& nodes2 = pullback.at(ls_2);
    //
    // // Compare the nodes within each level set
    // for (IntegerVector::const_iterator n1 = nodes1.begin(); n1 != nodes1.end(); ++n1){
    //   for (IntegerVector::const_iterator n2 = nodes2.begin(); n2 != nodes2.end(); ++n2){
    //     // Retrieve point indices within each node
    //     const IntegerVector& n1_idx = vertices[ std::to_string(*n1) ]; // access by vertex id
    //     const IntegerVector& n2_idx = vertices[ std::to_string(*n2) ]; // access by vertex id
    //
    //
    // // Enumerate the combinations. If any of subsets simplices don't exist as faces in the tree,
    // // skip.
    // apply_combinations(pids.size(), pids.size()-1, [](vector< size_t > idx){
    //
    // });
  // }
//   
//   // vector< node_ptr > simplices; 
//   // stree_ptr->traverse_max_skeleton(stree_ptr->root, [&simplices](node_ptr sigma, size_t d){
//   //   simplices.push_back(sigma);
//   // }, k-1);
//   // const size_t n_k_simplices = stree_ptr->n_simplexes.at(k - 1);
//   // apply_combinations(n_k_simplices, k+1, [&simplices](vector< size_t > idx) -> bool {
//   //   simplices.at()
//   // });
// }

// Builds the k-skeleton as a flag complex by inserting k-simplexes which have a 
// nonempty intersection between all pairs of vertices. 
//  k := the dimension to compute. Must be above 1. 
//  stree := SimplexTree object
// [[Rcpp::export]]
void build_flag_complex(const size_t k, SEXP stree){
  if (k < 2){ stop("'build_k_skeleton' is meant for building K-complexes for k > 1."); }
  Rcpp::XPtr<SimplexTree> stree_ptr(stree); // Collect the simplex tree
  using uint = unsigned int; 
  // Push paired connected edges into vector
  vector< size_t > connected; 
  stree_ptr->traverse_max_skeleton(stree_ptr->root, [&stree_ptr, &connected](const node_ptr e, const size_t d){
    vector< idx_t > edge = stree_ptr->full_simplex(e);
    std::pair< uint, uint > e_id = std::minmax(uint(edge.at(0)), uint(edge.at(1)));
    connected.push_back(szudzik_pair< uint, size_t >(e_id.first, e_id.second));
  }, 1);
  
  // Sort the connection pairs for faster checking
  std::sort(begin(connected), end(connected));
  
  // Shortcut lambda to check connectedness 
  const vector< idx_t > v = stree_ptr->get_vertices();
  const auto is_connected = [&connected](const size_t i) -> bool{
    return std::binary_search(begin(connected), end(connected), i);
  };
  const auto paired = [&is_connected, &v](const size_t i, const size_t j) -> bool{
    return is_connected(szudzik_pair< idx_t, size_t >(v.at(i), v.at(j)));
  };
  
  // Check all (n choose k) combinations of vertices
  size_t i = 0; 
  apply_combinations(v.size(), k+1, [&stree_ptr, &paired, &v, &i](vector< size_t > idx) -> bool {
    if ((i++) % 100){ Rcpp::checkUserInterrupt(); }
    const auto pw_nonempty = [&idx, &paired](vector< size_t > e) -> bool { return paired(idx[e[0]], idx[e[1]]); };
    const bool all_intersecting = apply_combinations(idx.size(), 2, pw_nonempty);
    
    // If all edges intersect, form a k-simplex
    if (all_intersecting){
      vector< size_t > sigma(idx.size());
      std::transform(begin(idx), end(idx), begin(sigma), [&v](const size_t ii){ return v[ii]; });
      stree_ptr->insert_simplex(sigma);
    }
    return true; 
  });
};


// Returns a graph representing the intersection between two lists of nodes.
// Given two lists of nodes, checks each pairwise combination of nodes to see if they intersect.
// If they contain a non-empty intersection, an edge is created.
// [[Rcpp::export]]
List intersectNodes(const List& nodes1, const List& nodes2, const IntegerVector& node_ids1, const IntegerVector& node_ids2) {
  using v_idx = vector< idx_t >;
  // Compress two groups of nodes for non-empty intersection
  std::list< IntegerVector > edgelist = std::list<IntegerVector>();
  int i = 0, j = 0;
  for (List::const_iterator n1 = nodes1.begin(); n1 != nodes1.end(); ++n1, ++i){
    j = 0;
    for (List::const_iterator n2 = nodes2.begin(); n2 != nodes2.end(); ++n2, ++j){
      // Retrieve point indices within each node
      const v_idx n1_idx = as< v_idx >(*n1);
      const v_idx n2_idx = as< v_idx >(*n2);

      // Add edge between the two if they share a data point
      bool intersect_check = nonempty_intersection(begin(n1_idx), end(n1_idx), begin(n2_idx), end(n2_idx));
      if (intersect_check){
        edgelist.push_back(IntegerVector::create(node_ids1[i], node_ids2[j]));
      }
    }
  }
  return(wrap(edgelist));
}


// Multiscale-version of updating the vertices in the level sets. Unlike the regular version, 
// the level sets are not explicitly stored in memory, but are rather extracted from the segments dynamically. 
// List update_level_sets(const IntegerVector which_levels, SEXP ms, const NumericMatrix& X, const Function f, List& vertices, List& ls_vertex_map, SEXP stree){
//   Rcpp::XPtr< MultiScale > ms_ptr(ms); // Cast the pointer to get the multiscale structure
//   
//   // // Shallow(?) copy the internals vectors 
//   // vector< IntegerVector > res(vertices.begin(), vertices.end());
//   std::map< int, IntegerVector > v_map = vertices_to_map(ls_vertex_map, vertices, stree);
// 
//   // Update the vertices in the level sets dynamically
//   for (const int index: which_levels){
//     const IntegerVector level_set = ms_ptr->extract_level_set(index) + 1; // convert to 1-based for R
//     update_level_set(index, level_set, X, f, v_map, ls_vertex_map, stree);
//   }
//   
//   // Return the new vertices
//   return(wrap(v_map));
// }

