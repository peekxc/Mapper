// skeleton.h
// Primary functions for building the k-skeletons
// Includes exported functions for building the skeletons both with and without the simplex tree. 
#include <Rcpp.h>
using namespace Rcpp;
#include <unordered_map>

#include "SimplexTree.h"
#include "MultiScale.h"

// Updates a single level set, whose flat index is given. This updates all of the corresponding vertices. 
void update_level_set(
    const int ls_flat_index,                 // The level set flat index. 
    const IntegerVector level_set,           // The level set. 
    const NumericMatrix& X,                  // The data 
    const Function f,                        // The clustering function 
    std::vector< IntegerVector >& vertices,  // Vertices of the Mapper. Updated by this method. 
    List& ls_vertex_map,                     // Level-set-index-to-vertex-id map. Updated by this method. 
    SEXP stree                               // Simplex tree external pointer. Updated by this method. 
); 

// Updates the vertices for the specified level sets. If all of the level sets are to be updated, this is
// equivalent to building the full 0-skeleton. 
List build_0_skeleton( 
    const IntegerVector which_levels,       // Vector of flat indexes indicating which level sets to use 
    const NumericMatrix& X,                 // The data 
    const Function f,                       // The clustering function 
    const List& level_sets,                 // All the level sets 
    List& vertices,                         // Vertices of the Mapper. Updated by this method. 
    List& ls_vertex_map,                    // Level-set-index-to-vertex-id map. Updated by this method. 
    SEXP stree                              // Simplex tree external pointer. Updated by this method.
); 

// Updates the edges for the specified level set pairs. If all of the level set pairs are to be updated, this is 
// equivalent to building the full 1-skeleton. 
void build_1_skeleton(
    const IntegerMatrix& ls_pairs,          // Matrix of level set index pairs to update. 
    const List& vertices,                   // Vertices of the Mapper. 
    const List& ls_vertex_map,              // Level-set-index-to-vertex-id map.
    SEXP stree                              // Simplex tree external pointer. Updated by this method.
);
  
// Multscale version
List update_level_sets(
    const IntegerVector which_levels, 
    SEXP ms, 
    const NumericMatrix& X, 
    const Function f, 
    List& vertices, 
    List& ls_vertex_map, 
    SEXP stree
);
