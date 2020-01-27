// skeleton.cpp
// Primary functions for building the k-skeletons
// Includes exported functions for building the skeletons both with and without the simplex tree. 
#include "skeleton.h"

#include "nerve_utility.h"
#include "neighborhood.h"

using str = std::string;

// Given a set of pullback indices and a simplex tree, for each pullback id, find the ids the
// of the other pullback sets whose corresponding vertices participate in a coface of the vertices 
// in the source pullback. 
// [[Rcpp::export]]
List connected_pullbacks(StringVector pullback_ids, const List& pullback, SEXP stree){
  Rcpp::XPtr< SimplexTree > stree_ptr(stree);

  // Get all the pullback ids
  StringVector pn = pullback.attr("names");
  vector< str > all_pids = as< vector< str > >(pn);
  
  // Reverse the pullback maps (key, value) pairs into a new map
  std::map< size_t, std::set< str > > v_to_pb;  
  for (str& pid: all_pids){
    const IntegerVector vids = as< IntegerVector >(pullback[pid]);
    for (size_t v: vids){ v_to_pb[v].insert( pid ); }
  }
  
  // Find pullback ids of adjacent vertices
  vector< str > pids = as< vector< str > >(pullback_ids);
  std::map< str, std::set< str > > pb_to_adj;
  for (str& pid: pids){
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
vector< IntegerVector > apply_clustering(str pid, const Function& level_set_f, const Function& cluster_f){
  // Rcout << "Running pullbakc preimage" << std::endl;
  const IntegerVector preimages = level_set_f(pid);
  if (preimages.size() == 0){ 
    vector<IntegerVector> res; 
    return(res); 
  } 
  // Rcout << "Running clustering" << std::endl;
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

// Given a vector of integer ids, converts them to strs, and removes any vertices which 
// have an id as in the given id list. 
List remove_by_id(IntegerVector ids, const List& vertices){
  if (Rf_isNull(vertices.names()) || vertices.size() == 0){ return vertices; }
  if (ids.size() == 0){ return vertices; }
  const size_t n = ids.size(); 
  vector< str > vnames = as< vector< str > >(vertices.names());
  LogicalVector subset = LogicalVector(vertices.size(), true);
  for (size_t i = 0; i < n; ++i){
    str id_str = std::to_string(ids.at(i));
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
  // vector< str > vnames = as< vector< str > >(vertices.names());
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
//     const str pid, 
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
//     vector< str > vnames = as< vector< str > >(vertices.names());
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
//       str v_key = std::to_string(v_id);
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

// std::map< str, IntegerVector > lst2map(const List& lst){
//   if (Rf_isNull(lst.names()) || lst.size() == 0){ return std::unordered_map< str, IntegerVector >(); }
//   vector< str > keys = as< vector< str > >(lst.names());
//   std::unordered_map< str, IntegerVector > res; 
//   for (auto c_key: keys){ 
//     IntegerVector v = as< IntegerVector >(lst[c_key]);
//     res.emplace(c_key, v); 
//   }
//   return(res);
// }

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
  auto str_cmp = [](const str& a, const str& b) -> bool{ 
    return(bool(std::stoi(a) < std::stoi(b)));
  };
  using str_map = std::map< str, IntegerVector, decltype(str_cmp) >; 
  
  // Extract the current set of vertices in C++ 
  str_map mod_vertices = str_map(str_cmp); 
  if (!Rf_isNull(vertices.names()) && as<List>(vertices.names()).size() > 0 ){
    vector< str > keys = as< vector< str > >(vertices.names());
    for (auto c_key: keys){ 
      IntegerVector v = as< IntegerVector >(vertices[c_key]);
      mod_vertices.emplace(c_key, v); 
    }
  }
  
  // Create lambda to extract the vertex ids as integers
  const auto vertex_ids = [&mod_vertices](){ 
    using v_type = str_map::value_type;
    if (mod_vertices.empty()){ return(vector< size_t >()); }
    vector< size_t > vids(mod_vertices.size());
    std::transform(begin(mod_vertices), end(mod_vertices), begin(vids), [](v_type v){
      return(std::stoi(v.first));
    });
    return(vids);
  };
  
  // Loop through the pullback ids to update. This will update the vertices, the pullback map, and the simplex tree
  vector< str > pids = as< vector< str > >(pullback_ids);
  for (str pid: pids){
    IntegerVector vids = pullback[pid];
    // if (vids.size() == 0){ continue; }
    
    // Remove vids in vertex list, pullback, and simplex tree
    // vertices = remove_by_id(vids, vertices);
    for (auto id: vids){ 
      mod_vertices.erase(std::to_string(int(id))); 
    }
    pullback[pid] = IntegerVector::create();
    
    // Apply the clustering to get the new vertices
    // Rprintf("clustering: %s\n", pid.c_str());
    vector< IntegerVector > new_vertices = apply_clustering(pid, level_set_f, cluster_f);
    //Rcout << "done clustering" << std::endl;
    
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
        str new_vid = std::to_string(new_0_simplexes.at(i));
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
    const str pid1 = as< str >(pullback_ids(i, 0)); 
    const str pid2 = as< str >(pullback_ids(i, 1));
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

  // If requested, track which simplices were added
  vector< vector< idx_t > > simplices_added; 
  
  // Iterate through the pullbacks 
  const size_t n = pullback_ids.nrow();
  for (size_t i = 0; i < n; ++i){
    if ((i % 100) == 0){ Rcpp::checkUserInterrupt(); }
    CharacterVector tmp = pullback_ids(i, _);
    vector< str > pids = as< vector< str > >(tmp);
    
    // If any of the pullback covers are empty, skip
    bool any_empty = std::accumulate(begin(pids), end(pids), false, [&pullback](bool any_empty, const str pid) -> bool{
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
    std::transform(begin(pids), end(pids), begin(sigma), [&pullback](const str pid){
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


template< typename T , typename I > 
T subset(const T& v, vector< I > idx){
  static_assert(std::is_integral<I>::value, "Must be integral type");
  T subset_(idx.size()); 
  std::transform(begin(idx), end(idx), begin(subset_), [&v](size_t pos) {
    return(v[pos]);
  });
  return(subset_);
}

// Based on: https://stackoverflow.com/questions/14945223/map-function-with-c11-constructs
template <typename T, typename Func>
auto map_f(const T& iterable, Func&& func) ->
  vector<decltype(func(std::declval<typename T::value_type>()))> {
  // Some convenience type definitions
  typedef decltype(func(std::declval<typename T::value_type>())) value_type;
  typedef vector<value_type> result_type;
  
  // Prepares an output vector of the appropriate size
  result_type res(iterable.size());
  
  // Let std::transform apply `func` to all elements
  // (use perfect forwarding for the function object)
  std::transform(
    begin(iterable), end(iterable), res.begin(),
    std::forward<Func>(func)
  );
  
  return res;
}

// template < typename T, typename U, typename Func> 
// auto apply_keyed(const T& map, const U& keys, Func&& f) ->
//   vector<decltype(func(std::declval<typename T::value_type>()))> {
//   // Typedefs 
//   typedef decltype(std::declval< typename U::value_type >()) elem_type; 
//   typedef decltype(func(std::declval<typename T::value_type>())) value_type;
//   typedef vector< value_type > result_type;
//   
//   // Prepares an output vector of the appropriate size
//   result_type res = map_f(keys, [&map, &f](const elem_type key){ 
//     return(f(map[key]));
//   });
//   return(res);
// }

// Convert a named list of integer vectors to an unordered_map
auto convert_pullback(const List& pullback) 
  -> std::unordered_map< str, IntegerVector >{

  // Retrieve names 
  CharacterVector pids = pullback.names();
  vector< str > index_set = as< vector< str > >(pids);
  
  // Convert named list to unordered_map
  std::unordered_map< str, IntegerVector > pb_map; 
  for (const str pid: index_set){
    if ((TYPEOF(pullback[pid]) == INTSXP) || (TYPEOF(pullback[pid]) == REALSXP)){ 
      pb_map.emplace(pid, as< IntegerVector >(pullback[pid]));
    } else { stop("SEXP type of pullback value not handled."); }
  }
  return(pb_map);
}
  
// nerve_functor 
// Creates a closure which accepts a vector of pullback indices. On evaluation, given a k-length 
// vector of pullback indices, k-combinations of nodes whose within each pullback are tested for 
// a non-empty intersection, provided their (k-1, k-2, ..., 0) simplices exist. 
std::function< void(vector< size_t >) > nerve_functor(
    const List& pullback, 
    List& vertices, 
    SEXP stree, 
    bool modify, 
    int threshold, 
    List& res
){
  using map_t = std::unordered_map< str, IntegerVector >;
  
  // Obtain (converted) pullback map
  map_t pb_map = convert_pullback(pullback);
  
  // Given a map-type container and a vector of ids, check if any of the 
  // containers mapped to by the id vectors are empty
  const auto any_empty = [](const map_t& C, const vector< str > ids) -> bool {
    const auto is_empty = [&C](const str key) -> bool { 
      auto el = C.find(key);
      bool has_v = (el != C.end() && el->second.size() > 0); 
      return(!has_v);
    };
    return std::any_of(begin(ids), end(ids), is_empty);
  };
  
  // Given a map-type container and a vector of ids, return the size 
  // of the containers mapped to by the id vectors. 
  const auto size_nested = [](const map_t& C, const vector< str > ids) -> const vector< size_t > {
    const auto size_of_vec = [&C](const str pid) -> size_t { 
      const auto v_it = C.find(pid);
      return size_t(v_it == C.end() ? 0 : v_it->second.size());
    };
    return map_f(ids, size_of_vec);
  };

  // Given a map-type container and a vector of ids and a vector of indices, return
  // a vector of the elements 
  const auto extract_inner = [](const map_t& C, const vector< str > ids, const vector< size_t > indices) -> vector< size_t >{
    if (ids.size() != indices.size()) { stop("Invalid extraction."); }
    const auto element = [&C](const str id, const size_t idx) -> size_t { 
      const auto v_it = C.find(id);
      assert(v_it != C.end());
      return size_t(v_it == C.end() ? 0 : v_it->second.at(idx));
    };
    
    const size_t k = ids.size();
    vector< size_t > simplex(k);
    for (size_t i = 0; i < k; ++i){
      simplex[i] = element(ids[i], indices[i]);
    }
    return(simplex);
  };
  
  // Convert the index set to a vector of strings
  const vector< str > index_set = as< vector< str > >(pullback.names());

  // Extract reference to simplex tree
  Rcpp::XPtr< SimplexTree > st_ptr(stree); // Collect the simplex tree
  SimplexTree& st = (*st_ptr);
    
  enum INTERSECTION_TYPE { NE_INSERT, INT_INSERT, NE_RETURN, INT_RETURN };
  INTERSECTION_TYPE int_check;
  if (threshold == 1 && modify){
    int_check = NE_INSERT;
  }
  else if (threshold == 1 && !modify){
    int_check = NE_RETURN;
  }
  else if (threshold > 1 && modify){
    int_check = INT_INSERT;
  }
  else {
    int_check = INT_RETURN;
  }
    
  // Create the closure to add a simplex to the complex, given (0-based) indices yielding 
  // the subset of the pullback to compare
  std::function< void(vector< size_t >) > add_simplex = [&st, &vertices, &res, threshold, pb_map, any_empty, index_set, size_nested, extract_inner, int_check](vector< size_t > pid_idx){
    
    // Get the subset to work with
    vector< str > c_pids = subset(index_set, pid_idx);
    // if (any_empty(pb_map, c_pids)){ Rcout << "empty pb" << std::endl; return; }

    // Otherwise, get the number of vertices per index
    const vector< size_t > pb_sizes = size_nested(pb_map, c_pids);
    
    // Check the product of the connected components between k-pairs of sets
    cart_prod(pb_sizes, [&st, &vertices, &c_pids, &pb_map, &extract_inner, &res, int_check, threshold](vector< size_t > idx){

      // Get the potential simplex to add
      vector< size_t > k_simplex = extract_inner(pb_map, c_pids, idx);

      // Prior to adding a k-simplex, check all (k+1 choose k) (k-1)-simplices exist
      const size_t k = k_simplex.size()-1;
      bool all_exists = apply_combinations(k+1, k, [&k_simplex, &st](const vector< size_t > idx){
        return(st.find_simplex(subset(k_simplex, idx)));
      });

      // If all the (k-1) faces exist, check the intersection criterion
      if (all_exists){
        
        // Get the ranges to work on
        auto v_pts = map_f(k_simplex, [&vertices](const size_t vid){
          const IntegerVector vertex = vertices[ std::to_string(vid) ];
          return std::make_pair(vertex.begin(), vertex.end());
        });

        // Add a simplex if there's a nonempty intersection between all vertices
        switch(int_check){
          case NE_INSERT: {
            if (nfold_nonempty(v_pts)){
              st.insert_simplex(k_simplex);
            }
            break;
          }
          case NE_RETURN: {
            if (nfold_nonempty(v_pts)){
              res.push_back(k_simplex);
            }
            break; 
          }
          case INT_INSERT:
            if (nfold_intersection(v_pts).size() >= threshold){
              st.insert_simplex(k_simplex);
            }
            break;
          case INT_RETURN:
            if (nfold_intersection(v_pts).size() >= threshold){
              res.push_back(k_simplex);
            }
            break; 
        }
      }
    });
  };
  return(add_simplex);
}

// Accepts a generic matrix of pullback ids, then computes the Cech nerve.
// [[Rcpp::export]]
List build_k_skeleton_ids(CharacterMatrix pullback_ids, const List& pullback, List& vertices, SEXP stree, const bool modify, const size_t threshold){
  using simplex_f = std::function< void(vector< size_t >) >;

  // Extract an enclosed simplex insertion function
  List resulting_simplices = List();
  const simplex_f nerve_f = nerve_functor(pullback, vertices, stree, modify, threshold, resulting_simplices); 

  // Returns the 0-based index of the str pullback id
  vector< str > index_set = as< vector< str > >(pullback.names());
  const auto index_of = [&index_set](str pid) -> size_t {
    auto it = std::find(begin(index_set), end(index_set), pid);
    return(static_cast< size_t >(std::distance(begin(index_set), it)));
  };

  // Casts given row i as a vector of std::string's
  const auto extract_pids = [&pullback_ids](const size_t i) -> vector< str > {
    CharacterVector ids = pullback_ids.row(i);
    return(as< vector< str > >(ids));
  };

  // Loop through the pullback ids
  const size_t n = pullback_ids.nrow();
  for (size_t i = 0; i < n; ++i){
    vector< size_t > indices = map_f(extract_pids(i), index_of);
    nerve_f(indices);
  }

  return(resulting_simplices);
}

// Also computes the nerve of a cover, however allows for generic functions-producing iterators 
// to be passed as well. 
// [[Rcpp::export]]
List build_k_skeleton_gen(SEXP subset_sexp, const List& pullback, List& vertices, SEXP stree, const bool modify, const size_t threshold){
  XPtr< del_f > subset_fun(subset_sexp); 
  del_f subset_f = *subset_fun;
  
  // Extract an enclosed simplex insertion function
  List resulting_simplices = List();
  const set_f nerve_f = nerve_functor(pullback, vertices, stree, modify, threshold, resulting_simplices); 
  
  // Send the nerve function to the index generator
  subset_f(nerve_f);
  
  // Return the added simplices
  return(resulting_simplices);
}

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

