// MultiScale.h 
// Header file for the multiscale indexing structure. 

#include <Rcpp.h>
using namespace Rcpp;

#include "GridIndex.h"
#include "utility.h"
#include "utility_rcpp.h"
#include <memory> // smart pointers

// Using typenames 
using u8 = uint_fast8_t;
using index_t = std::vector<u8>;
using pdiff = std::ptrdiff_t;
using l_index_t = std::vector< std::ptrdiff_t >;
using rle = std::pair< l_index_t, index_t >;

template <typename T> 
using u_ptr = std::unique_ptr<T>;

// Each point requires the following 4 pieces of information
struct path_info {
  u8 k_idx;       // (constant) index into key map giving ls path
  u8 p_idx;       // 'previous' index in path; reset to c_idx after each move 
  u8 c_idx;       // current index in path
  u8 c_segment;   // current segment the point lies in
}; 

// Struct that simplifies computing level set intersections
struct pt_update {
  const pdiff idx;
  std::vector< u8 > min_ls, max_ls;
  std::vector< u8 > target_segment;
  pt_update(const pdiff _id, const int d) : idx(_id){
    min_ls = std::vector< u8 >(d, std::numeric_limits<u8>::max() );
    max_ls = std::vector< u8 >(d, std::numeric_limits<u8>::min() );
    target_segment = std::vector< u8 >(d);
  }
  void update_min_max(u8 lb, u8 ub, const int d_i){
    if (lb > ub){ std::swap(lb, ub); }
    if (lb < min_ls.at(d_i)){ min_ls.at(d_i) = lb; }
    if (ub > max_ls.at(d_i)){ max_ls.at(d_i) = ub; }
  }
};

// MultiScale class
// The MultiScale class effectively creates a filtration of point indices per dimension. Given a Mapper is already computed, 
// this class provides tools for extracting which components of the 0-skeleton and 1-skeleton need to be updated in order 
// to compute a Mapper at a different parameterization.  
struct MultiScale {
  // Constants
  const std::size_t d, n;                              // Dimensionality of filter space + size of data      
  const index_t num_intervals;                         // The number of intervals to distribute per dimension
  
  // Information needed to know what state the Mapper is in.
  l_index_t filt_idx;                                  // The current index into the filtration [ O(d) ]
  std::vector< l_index_t > filt_index_set;             // The filtration index set [ O(nkd) ]
  std::vector< NumericVector > filt_dist;              // The distances corresponding to each filtration index [ O(nkd) ]  TODO: replace w/ O(2nd) version

  // Information needed to map level sets --> segments 
  index_t c_ls_segment_idx;                            // Tracks the index into the level set segment map [ O(d) ] 
  std::vector< index_t > ls_segment_idx;               // Vector giving the segment indices spanned by each level set [ O(2kd) ]
  std::vector< l_index_t >  ls_change_idx;             // Tracks when the level set segment indices change (cumulative run-length encoding) [ O(kd) ]
  
  // Information needed to track where the points go as the filtration changes
  std::vector< std::vector< path_info > > pt_info;     // Point position information [ O(n) ] 
  std::vector< std::vector< index_t > > ls_paths;      // Unique paths taken by any given point [ O(2kd) ]

  // Global information not indexed by dimension
  GridIndex< u8 > ls_grid;                             // Provides mappings from LSMI --> LSFI and vice versa. [ O(k^d) ]
  std::map< index_t, l_index_t > segment_map;          // Mapping from segment index --> pt ids [ O(k^d + n) ]
  index_t d_range;                                     // Stores the sequence of dimension indices [ O(d) ]
  
  // // Edgelist to track to 
  // std::vector< std::pair< pdiff, pdiff > > 
  
  // Member functions
  MultiScale(const int n_pts, IntegerVector resolution);
  SEXP as_XPtr();
  IntegerVector extract_level_set(const int lsfi);
  List get_segment_map();
  void create_filtration(const IntegerVector& f_idx, const NumericVector& intervals, const int d_i);
  IntegerVector get_nearest_filtration_index(NumericVector intervals);
  void set_filtration_rle(const IntegerVector& ls_changes, const int d_i);
  index_t extract_segment(const std::size_t i);
  void insert_pts(const IntegerMatrix& A, const IntegerMatrix& pt_ls_path);
  void create_ls_paths(const IntegerMatrix& _ls_paths, const int d_i);
  void update_ls_segment_idx(const std::size_t i, const int d_i);
  index_t compute_ls_segment_idx(std::size_t i, std::size_t d_i);
  List update_segments(const IntegerVector target_idx);
};


