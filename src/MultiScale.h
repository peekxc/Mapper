// MultiScale.h 
// Header file for the multiscale indexing structure. 
#ifndef MULTISCALE_H
#define MULTISCALE_H

#include <Rcpp.h>
using namespace Rcpp;

#include "utility_rcpp.h"
#include "GridIndex.h" // faster grid indexing structure
#include "skeleton.h" // skeleton update functions, simplex tree, etc.
#include <memory> // smart pointers

// Shortcut typenames
using v_uint8_t = std::vector< uint8_t >;       // vector for compact unsigned int storage
using v_sidx_t = std::vector< sidx_t >;         // vector for signed indexes 
using rle = std::pair< v_sidx_t, v_uint8_t >;   // run length encoding structure

// Each point requires the following 4 pieces of information
struct path_info {
  uint8_t k_idx;       // (constant) index representing key that returns the level set path
  uint8_t c_idx;       // current index *into* path
  uint8_t p_idx;       // 'previous' index in path; resets to c_idx after each move 
  uint8_t c_segment;   // current segment the point lies in
}; 

// Struct that simplifies computing level set intersections
struct pt_update {
  const sidx_t idx;
  std::vector< uint8_t > min_ls, max_ls;
  std::vector< uint8_t > target_segment;
  pt_update(const sidx_t _id, const uidx_t d) : idx(_id){
    min_ls = std::vector< uint8_t >(d, std::numeric_limits<uint8_t>::max() );
    max_ls = std::vector< uint8_t >(d, std::numeric_limits<uint8_t>::min() );
    target_segment = std::vector< uint8_t >(d);
  }
  void update_min_max(uint8_t lb, uint8_t ub, const uidx_t d_i){
    if (lb > ub){ std::swap(lb, ub); }
    if (lb < min_ls.at(d_i)){ min_ls.at(d_i) = lb; }
    if (ub > max_ls.at(d_i)){ max_ls.at(d_i) = ub; }
  }
};

// MultiScale class
// The MultiScale class effectively creates a filtration of point indices per dimension. Given a Mapper is already computed, 
// this class provides tools for extracting which components of the 0-skeleton and 1-skeleton need to be updated in order 
// to efficiently compute a Mapper at a different parameterization.  
struct MultiScale {
  // Constants
  const std::size_t d, n;                                 // Dimensionality of filter space + size of data      
  const v_uint8_t num_intervals;                          // The number of intervals to distribute per dimension
  
  // Needed to know what state the Mapper is in.
  // TODO: replace below w/ O(2nd) version when overlap < 50%
  v_sidx_t filt_idx;                                      // The current index into the filtration [ O(d) ]
  std::vector< v_sidx_t > filt_index_set;                 // The filtration index set [ O(nkd) ]
  std::vector< NumericVector > filt_dist;                 // The distances corresponding to each filtration index [ O(nkd) ]  

  // Needed to map level sets --> segments 
  v_uint8_t c_ls_segment_idx;                             // Tracks the index into the level set segment map [ O(d) ] 
  std::vector< v_uint8_t > ls_segment_idx;                // Vector giving the segment indices spanned by each level set [ O(3kd) ]
  std::vector< v_sidx_t >  ls_change_idx;                 // Tracks when the level set segment indices change (cumulative run-length encoding) [ O(kd) ]
  
  // Needed to track paths followed by points
  std::vector< std::vector< path_info > > pt_info;        // Point position information [ O(n) ] 
  std::vector< std::vector< v_uint8_t > > ls_paths;       // Unique paths taken by any given point [ O(2kd) ]

  // Global information not indexed by dimension
  GridIndex< uint8_t > ls_grid;                           // Provides mappings from LSMI --> LSFI and vice versa. [ O(k^d) ]
  std::map< v_uint8_t, v_sidx_t > segment_map;            // Mapping from segment index --> pt ids [ O(k^d + n) ]
  v_uint8_t d_range;                                      // Stores the sequence of dimension indices [ O(d) ]
  
  // Vertices to update. Only kept in the struct for efficiency. 
  std::map< int, IntegerVector > vertices;                // Stores the points associated with each vertex. [ O(nk^d) ]
  
  // Member functions
  MultiScale(const uidx_t n_pts, IntegerVector resolution);
  SEXP as_XPtr();
  IntegerVector extract_level_set(const uidx_t lsfi);
  IntegerMatrix point_info(const uidx_t d_i);
  IntegerMatrix uniq_paths(const uidx_t d_i);
  List get_segment_map();
  void create_filtration(const IntegerVector& f_idx, const NumericVector& intervals, const uidx_t d_i);
  IntegerVector get_nearest_filtration_index(NumericVector intervals);
  void set_filtration_rle(const IntegerVector& ls_changes, const uidx_t d_i);
  v_uint8_t extract_segment(const sidx_t i);
  void insert_pts(const IntegerMatrix& A, const IntegerMatrix& pt_ls_path);
  void create_ls_paths(const IntegerMatrix& _ls_paths, const uidx_t d_i);
  void update_ls_segment_idx(const sidx_t i, const uidx_t d_i);
  v_uint8_t compute_ls_segment_idx(sidx_t i, uidx_t d_i);
  List update_segments(const IntegerVector target_idx);
  void initialize_vertices(List& ls_vertex_map, List& vertices);
  void update_vertices(const IntegerVector, const NumericMatrix&, const Function, List&, SEXP);
};

#endif
