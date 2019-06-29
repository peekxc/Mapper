// MultiScale.h 
// Header file for the multiscale indexing structure. 
#ifndef MULTISCALE_H
#define MULTISCALE_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(simplextree)]]
#include "skeleton.h" // skeleton update functions, simplex tree, etc.
#include "utilities.h" // many useful utilities
#include <memory> // smart pointers
#include <cstdint> // integer types

// Typenames aliases
using std::size_t;
using std::vector; 
using int_t = int_fast32_t;
using uint_t = uint_fast32_t;
using v_int_t = vector< int_t >;         // vector for signed integers 
using v_uint_t = vector< uint_t >;       // vector for unsigned integers
using rle = std::pair< v_int_t, v_uint_t >;   // run length encoding structure
// using int_t = std::ptrdiff_t;

// Each point requires the following 4 pieces of information (per dimension)
struct path_info {
  uint_t k_idx;       // (constant) index representing key that returns the level set path
  uint_t c_idx;       // current index *into* path
  uint_t p_idx;       // 'previous' index in path; resets to c_idx after each move 
  uint_t c_segment;   // current segment the point lies in
}; 

// Struct that simplifies computing level set intersections
struct pt_update {
  const int_t idx;
  vector< uint_t > min_ls, max_ls;
  vector< uint_t > target_segment;
  pt_update(const int_t id, const uint_t d) : idx(id){
    min_ls = vector< uint_t >(d, std::numeric_limits< uint_t >::max() );
    max_ls = vector< uint_t >(d, std::numeric_limits< uint_t >::min() );
    target_segment = vector< uint_t >(d);
  }
  void update_min_max(uint_t lb, uint_t ub, const uint_t d_i){
    if (lb > ub){ std::swap(lb, ub); }
    if (lb < min_ls.at(d_i)){ min_ls.at(d_i) = lb; }
    if (ub > max_ls.at(d_i)){ max_ls.at(d_i) = ub; }
  }
};

// Used to hash a vector as a unique key for a set
template <typename T>
struct VectorHash {
  static_assert(std::is_integral<T>::value, "Integral-type required as a range storage type.");
  size_t operator()(const std::vector<T>& v) const {
    std::hash<T> hasher;
    size_t seed = 0;
    for (T i : v) { seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2); }
    return seed;
  }
};

// MultiScale class
// The MultiScale class effectively creates a filtration of point indices per dimension. Given a Mapper is already computed, 
// this class provides tools for extracting which components of the 0-skeleton and 1-skeleton need to be updated in order 
// to efficiently compute a Mapper at a different parameterization.  
struct MultiScale {
  // Constants
  size_t d, n;                                 // Dimensionality of filter space + size of data      
  v_uint_t num_intervals;                          // The number of intervals to distribute per dimension
  
  // Needed to know what state the cover is in.
  // TODO: replace below w/ O(2nd) version when overlap < 50%
  v_int_t current_index;                            // The current index into the filtration [ O(d) ]
  vector< v_int_t > indices;                 // The filtration index set [ O(nkd) ]
  vector< NumericVector > eps;                 // The distances corresponding to each filtration index [ O(nkd) ]  

  // The storage of the points
  std::map< v_uint_t, v_int_t > segments;            // Container mapping (segment) index --> pt ids [ O(k^d + n) ]
  
  // Needed to map segments -> level sets 
  vector< v_int_t > cc_offsets;                 // Canonical cover offset: partitions the points into canonical covers (via cumulative run-length encoding) [ O(kd) ]
  v_uint_t cc_index;                            // Canonical cover index: Tracks the index of the current canonical cover (used by the segment map) [ O(d) ] 
  vector< v_uint_t > canonical_cover;   // The canonical cover associated with 'current_index' mapping segment indices -> level sets [ O(3kd) ]

  // Needed to track paths followed by points
  vector< vector< path_info > > pt_info;        // Point position information [ O(n) ] 
  vector< vector< v_uint_t > > ls_paths;       // Unique paths taken by any given point [ O(2kd) ]

  // Global information not indexed by dimension
  // GridIndex< uint_t > ls_grid;                           // Provides mappings from LSMI --> LSFI and vice versa. [ O(k^d) ]
  // vector< uint_t > index_AoS; 
  
  v_uint_t get_previous_lsmi(uint_t);
  v_uint_t get_current_lsmi(uint_t);
  
  
  v_uint_t d_range;                                   // Stores the sequence of dimension indices [ O(d) ]
  
  // Vertices to update. Only kept in the struct for efficiency. 
  std::map< int, IntegerVector > vertices;                // Stores the points associated with each vertex. [ O(nk^d) ]
  
  // Precomputed indices to simplify index arithmetic 
  vector< size_t > cum_prod; 

  // Member functions
  MultiScale(const List);
  SEXP as_XPtr();
  vector< int_t > extract_level_set(const size_t lsfi);
  IntegerMatrix point_info(const uint_t d_i);
  IntegerMatrix uniq_paths(const uint_t d_i);
  List get_segments();
  List get_segment_table();
  
  template <typename T>
  std::string vec_to_string(vector<T>, std::string);
  
  // void create_filtration(const IntegerVector& f_idx, const NumericVector& intervals, const uint_t d_i);
  // void set_filtration_rle(const IntegerVector& ls_changes, const uint_t d_i);
  // void create_ls_paths(const IntegerMatrix& _ls_paths, const uint_t d_i);
  IntegerVector get_nearest_index(NumericVector intervals);

  
  v_uint_t get_segment_midx(const int_t);
  void insert_pts(const IntegerMatrix&, const IntegerMatrix&);

  // Utility functions to convert between flat and multi indexes
  size_t multi_to_flat(vector< uint_t >);
  vector< uint_t > flat_to_multi(size_t);
  
  // Computes row 'i' of the segment table for dimension d_i. Requires k and d. 
  v_uint_t segment_cover_idx(int_t i, uint_t d_i);
  
  
  void update_canonical_cover(const int_t, const uint_t);
  void update_segments(std::unordered_map< size_t, vector< uint_t > >&);
  vector< v_uint_t > resolve_paths(const size_t);
  std::unordered_map< size_t, vector< uint_t > > update_paths(const IntegerVector);
  void modified_indices(const std::unordered_map< size_t, vector< uint_t > >&, 
                        const size_t, 
                        std::unordered_set< size_t >&, 
                        std::unordered_set< vector< uint_t >, VectorHash< uint_t > >&, 
                        const bool);
  List update_index(const IntegerVector, const size_t);
  
  // void initialize_vertices(List& ls_vertex_map, List& vertices);
  // void update_vertices(const IntegerVector, const NumericMatrix&, const Function, List&, SEXP);
};

#endif
