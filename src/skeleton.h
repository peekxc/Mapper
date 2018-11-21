// skeleton.h
// Primary functions for building the k-skeletons
// Includes exported functions for building the skeletons both with and without the simplex tree. 
#ifndef SKELETON_H
#define SKELETON_H

#include <Rcpp.h>
using namespace Rcpp;
#include <unordered_map>

#include "SimplexTree.h"

// Updates a single level set, whose flat index is given. This updates all of the corresponding vertices. 
void update_level_set(
    const int ls_flat_index,                    // The level set flat index. 
    const IntegerVector level_set,              // The level set. 
    const NumericMatrix& X,                     // The data 
    const Function f,                           // The clustering function 
    std::map< int, IntegerVector >& vertices,   // Vertices of the Mapper. Updated by this method. 
    List& ls_vertex_map,                        // Level-set-index-to-vertex-id map. Updated by this method. 
    SEXP stree                                  // Simplex tree external pointer. Updated by this method. 
); 

// Converts a given list of vertices (integer vectors) to a C++ mapping to integer vectors. 
std::map< int, IntegerVector > vertices_to_map(
    List& ls_vertex_map,                        // Mapping between level set index --> vertex ids 
    List& vertices                              // Mapping between vertex ids --> points 
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
  
#endif
