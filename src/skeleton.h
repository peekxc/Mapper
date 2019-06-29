// skeleton.h
// Primary functions for building the k-skeletons
// Includes exported functions for building the skeletons both with and without the simplex tree. 
#ifndef SKELETON_H
#define SKELETON_H

#include <Rcpp.h>
using namespace Rcpp;

// Shim to get the simplextree definitions exposed 
// [[Rcpp::depends(simplextree)]]
#include "simplextree.h"
#include "utilities.h"
#include <unordered_map>

// Aliases
using std::vector; 
using std::begin; 
using std::end; 

// Updates a single level set, whose flat index is given. This updates all of the corresponding vertices. 
void update_pullback(
    const std::string pid, 
    const Function level_set_f, 
    const NumericMatrix& X, 
    const Function cluster_f, 
    List& pullback, 
    List& vertices, 
    SEXP stree
);
// // Converts a given list of vertices (integer vectors) to a C++ mapping to integer vectors. 
// std::map< int, IntegerVector > vertices_to_map(
//     List& ls_vertex_map,                        // Mapping between level set index --> vertex ids 
//     List& vertices                              // Mapping between vertex ids --> points 
// );

// Updates the vertices for the specified level sets. If all of the level sets are to be updated, this is
// equivalent to building the full 0-skeleton. 
List build_0_skeleton( 
    const StringVector pullback_ids,       // Vector of flat indexes indicating which level sets to use 
    // const NumericMatrix& X,                 // The data 
    const Function cluster_f,               // The clustering function 
    const Function level_set_f,             // Function which returns the level sets. 
    List& vertices,                         // Vertices of the Mapper. Updated by this method. 
    List& pullback,                         // Pullback map. Updated by this method. 
    SEXP stree                              // Simplex tree external pointer. Updated by this method.
); 

// Updates the edges for the specified level set pairs. If all of the level set pairs are to be updated, this is 
// equivalent to building the full 1-skeleton. 
void build_1_skeleton(const CharacterMatrix& pullback_ids, const int min_sz, const List& vertices, const List& pullback, SEXP stree);
  
#endif
