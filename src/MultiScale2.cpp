#include <Rcpp.h>
using namespace Rcpp;

#include "GridIndex.h"
#include "utility.h"
#include <memory> // smart pointers

// Using typenames 
using u8 = uint_fast8_t;
using index_t = std::vector<u8>;
using pdiff = std::ptrdiff_t;
using l_index_t = std::vector< std::ptrdiff_t >;
using rle = std::pair< l_index_t, index_t >;

template <typename T> 
using u_ptr = std::unique_ptr<T>;

// For each point, store these 3 
struct path_info {
  u8 k_idx; // (constant) index into key map giving ls path
  u8 p_idx; // 'previous' index in path; reset to c_idx after each move 
  u8 c_idx; // current index in path
  u8 c_segment; // current segment the point lies in
  
  // path_info(const int k) : k_idx(static_cast<u8>(k)){
  //   
  // }
  // extract
}; // O(n) of these needed

// Simple struct to organize what information needs to be updated
// template <std::size_t d>
// struct pt_update {
//   pdiff id;
//   u8 min_source_ls[d], max_target_ls[d];
//   u8 source_segment[d];
// };
// 
// 
// pt_update<>* make_pt_update(const int d){
//   switch(d){
//   case 15: return(pt_update<15>());
//   case 14: return(pt_update<14>());
//   case 13: return(pt_update<13>());
//   case 12: return(pt_update<12>());
//   case 11: return(pt_update<11>());
//   case 10: return(pt_update<10>());
//   case 9: return(pt_update<9>());
//   case 8: return(pt_update<8>());
//   case 7: return(pt_update<7>());
//   case 6: return(pt_update<6>());
//   case 5: return(pt_update<5>());
//   case 4: return(pt_update<4>());
//   case 3: return(pt_update<3>());
//   case 2: return(pt_update<2>());
//   case 1: return(pt_update<1>());
//   default: stop("Dimension must be <= 15.");
//   }
// }

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

struct MultiScale2{

    // TODO: investigate this 
  // Constants: dimensionality + size of data 
  const std::size_t d, n;
  const index_t num_intervals; 
    
  // All the information needed per dimension
  std::vector< index_t > pt_idx;                  // The point ids changing level sets at each filtration step [ O(nkd) ]
  std::vector< index_t > origin_ls;               // The original level sets each point started in [ O(nd) ]
  std::vector< index_t > source_ls;               // Cached vector storing the current level sets each point is associated with [ O(nd) ]
  std::vector< index_t > source_segment;          // Cached vector storing the current segment each point is associated with [ O(nd) ]
  std::vector< index_t > target_ls;               // the target level sets in the filtration [ O(nkd) ]
  std::vector< index_t > ls_segment_idx;          // vector giving the segment indices spanned by each level set [ O(2kd) ]
  index_t c_ls_segment_idx;                       // Tracks the index into the level set segment matching [ O(d) ] 
  std::vector< l_index_t >  ls_change_idx;        // tracks when the level set segment indices change (run-length encoding) [ O(kd) ]
  l_index_t filtration_idx;                       // the current index of the filtration [ O(d) ]
  
  // New idea
  std::vector< std::vector< path_info > > pt_info; // point position information [ O(n) of these needed ] 
  std::vector< std::vector< index_t > > ls_paths; // unique paths taken by any given point
  std::vector< l_index_t > filtration; // the indices corresponding to the filtration
  
  // Global information not indexed by dimension
  GridIndex< u8 > ls_grid; // Structure that provides quick mappings from LSMI --> LSFI and vice versa.
  std::map< index_t, l_index_t > segment_map;     // mapping from segment index --> pt ids [ O(k^d + n) ]
  index_t d_range; 
  
  MultiScale2(const int n_pts, IntegerVector resolution) 
    : n(n_pts), d(resolution.size()), num_intervals(resolution.begin(), resolution.end()), ls_grid(GridIndex< u8 >(resolution)) {
    d_range = index_t(d);
    std::iota(d_range.begin(), d_range.end(), 0);
    
    pt_info = std::vector< std::vector< path_info > >(d);
    ls_paths = std::vector< std::vector< index_t > >(d);
    filtration = std::vector< l_index_t >(d);
    ls_change_idx = std::vector< l_index_t >(d);
    ls_segment_idx = std::vector< index_t >(d);
    c_ls_segment_idx = index_t(d);
    
    filtration_idx = l_index_t(d, -1);
      
    for (auto d_i: d_range){ update_ls_segment_idx(0, d_i); }
    // Any point not in the map is in its default position
    // std::for_each(pt_ids.begin(), pt_ids.end(), [&](const int pt){
    //   std::map<int, int>::iterator it = pt_to_res.find(pt);
    //   if (it == pt_to_res.end()){ // pt wasn't found in the map
    //     // Extract first occurrence of point in point ids. It's from level set is where it came from originally.
    //     std::size_t first_pt_idx = std::distance(c_pt_idx.begin(), std::find(c_pt_idx.begin(), c_pt_idx.begin(), pt));
    //     int first_ls = from_ls_idx.at(d_i).at(first_pt_idx);
    //     pt_to_res.emplace_hint(it, pt, first_ls); // every point initially exists in the lower segment of its initial level set
    //   }
    // });
    
    // create_ls_paths()
    // insert_pts(A, pt_ls_path); // Create the points and their associated paths
    
  }
  
  // Expects 1-based RLE built from the point swap indices. Lengths should be cumulative.
  void set_filtration_rle(const IntegerVector& ls_changes, const int d_i){
    ls_change_idx.at(d_i) = l_index_t(ls_changes.begin(), ls_changes.end());
  }
  
  // Expects a 0-based vector of integers representing the total order that point indices 
  // along dimension 'd_i' change level sets.
  void create_filtration(const IntegerVector& f_idx, const int d_i){
    filtration.at(d_i) = l_index_t(f_idx.begin(), f_idx.end());
  }
  
  index_t extract_segment(const std::size_t i){
    if (i < 0 || i >= n){ stop("Invalid i"); }
    index_t res(d);
    for (auto d_i: d_range){ res.at(d_i) = pt_info.at(d_i).at(i).c_segment; }
    return(res);
  }
  
  // A contains the 1-based level set indexes each point lies in when the cover is disjoint
  // pt_ls_path contains the 1-based index of the unique path taken by 
  void insert_pts(const IntegerMatrix& A, const IntegerMatrix& pt_ls_path){
    if (A.nrow() != n || A.ncol() != d){ stop("A must be a matrix of size equal to the number of points in the filtration."); }
    if (pt_ls_path.nrow() != n || pt_ls_path.ncol() != d){ stop("Dimensionality of A must equal dimensionality of filtration."); }
    
    // Clears any points 
    for (auto d_i: d_range){
      pt_info.at(d_i).clear(); 
    }
    
    // Add the points to the structure 
    for (std::size_t i = 0; i < n; ++i){
      // IntegerVector tmp = A.row(i) - 1, tmp2 = pt_ls_path.row(i) - 1;
      // index_t c_ls = index_t(tmp.begin(), tmp.end()); // current level set 
      // l_index_t c_path = l_index_t(tmp2.begin(), tmp2.end()); // current ls path
      
      // Assemble the initial path information to each dimension
      path_info new_pt = path_info(); 
      for (std::size_t d_i = 0; d_i < d; ++d_i){
        new_pt.k_idx = static_cast<int>(pt_ls_path(i, d_i)) - 1; // the ls path the point follows along the current dimension
        new_pt.c_idx = new_pt.p_idx = 0; // every point starts at its originating level set
        new_pt.c_segment = (static_cast<int>(A(i, d_i)) - 1) * 2; // every point starts in its originating level sets lower segment
        pt_info.at(d_i).push_back(new_pt);
      }
      
      // Add point to the segment map 
      index_t pt_segment = extract_segment(i);
      segment_map[pt_segment].push_back(i);
    }
  }
  
  // Stores a vector of the distinct ls paths taken by an individual point along dimension d_i. 
  void create_ls_paths(const IntegerMatrix& _ls_paths, const int d_i){
    const std::size_t _n = _ls_paths.nrow(), _k = _ls_paths.ncol(); 
    if (_k != num_intervals.at(d_i)){ stop("Nope"); }
    ls_paths.at(d_i).clear();
    for (std::size_t i = 0; i < _n; ++i){
      IntegerVector tmp = _ls_paths.row(i);
      index_t c_path = index_t(tmp.begin(), tmp.end());
      ls_paths.at(d_i).push_back(c_path);
    }
  }
  
  // Given a global filtration index 'i' \in { -1, 0, 1, ..., n*d_i - 1 } and dimension 'd_i', 
  // updates the 'ls_segment_idx' and 'c_ls_segment_idx' variables.
  void update_ls_segment_idx(const std::size_t i, const int d_i){
    const l_index_t& cum_lengths = ls_change_idx.at(d_i);
    std::size_t ii;
    if ( i  == -1){
      ii = 0; 
    } else {
      auto ub = std::upper_bound(cum_lengths.begin(), cum_lengths.end(), i);
      ii = std::distance(cum_lengths.begin(), ub);
    }
    // if (ii > 0){ ii -= 1; }
    // Rprintf("swap idx: %d, c_segment_idx: %d \n", ii, c_ls_segment_idx.at(d_i));
    if (ii != c_ls_segment_idx.at(d_i)){
      ls_segment_idx.at(d_i) = compute_ls_segment_idx(ii, d_i);
      c_ls_segment_idx.at(d_i) = ii;
    }
    // Rcout << "here_extra" << std::endl;
  }
  
  // Creates a vector with the range [i, j]
  template <typename T> 
  std::vector<T> seq_ij(const int i, const int j){
    std::size_t sz = std::abs(j - i)+1;
    std::vector<T> rng = std::vector<T>(sz);
    std::iota(rng.begin(), rng.end(), i);
  }

  index_t compute_ls_segment_idx(std::size_t i, std::size_t d_i){
    if (d_i >= num_intervals.size()){ stop("Invalid index: Must be 0-based index less than filter dimensionality"); }
    const int n_level_sets = num_intervals.at(d_i); 
    const int n_endpts = n_level_sets*2;
    index_t ls_idx_res = index_t(n_endpts);
    int hc = n_level_sets - 1, lc = 0;
    for (std::size_t j = 0; j < n_level_sets; ++j){
      int j_tmp = n_endpts - j - 1;
      // Rprintf("j: %d, j_tmp: %d, i: %d, lc: %d, hc: %d\n");
      if (j % 2 == 0){
        ls_idx_res[j] = u8(i >= lc ? j - lc : j - i);
        ls_idx_res[j_tmp] = u8(i >= lc ? j_tmp + lc : j_tmp + i);
        ++lc;
      }
      else {
        ls_idx_res[j] = u8(i >= hc ? j + hc : j + i);
        ls_idx_res[j_tmp] = u8(i >= hc ? j_tmp - hc : j_tmp - i);
        --hc;
      }
    }
    return(ls_idx_res);
  }
  
  List update_segments2(const IntegerVector target_idx){
    
    // The information needed to update the filtration
    std::vector< index_t > updated_source_segments(d), updated_target_segments(d); 
    std::vector< index_t > updated_source_ls(d), updated_target_ls(d);
    
    // pt_update<d> test; 
    
    // The information needed to update the filtration
    std::map< pdiff, pt_update > pts_to_update = std::map< pdiff, pt_update >(); 
    // pt_update pt_to_update; 
    
    // Aggregate per-dimension information
    for (auto d_i: d_range){
      const pdiff c_idx = filtration_idx.at(d_i); // current filtration index
      const pdiff t_idx = target_idx.at(d_i); // target filration index
      if (c_idx == t_idx){ continue; } // Don't update if we're at the current index
      const bool expanding = c_idx < t_idx; // are we expanding or contracting the current dimension?
      const l_index_t c_filt_idx = filtration.at(d_i); // the filtration itself
    
      
      // Rprintf("d=%d: expanding? %d. Current index: %d, target index: %d\n", d_i, expanding, c_idx, t_idx);
      if (expanding){
        
        // Compute the filtration updates
        for(pdiff i = c_idx+1; i <= t_idx; ++i){
          pdiff f_i = c_filt_idx.at(i); // filtration index
          // pdiff k_i = std::ceil(f_i/n) - 1; // relative index into the path 
          pdiff pt_idx = (f_i % n); // the point changing at this step TODO: use a positive-only modulus
          
          // Debugging
          // Rprintf("(%d=%d), pt idx: %d\n", i, f_i, pt_idx);
          
          // Extract the path information
          path_info& c_path = pt_info.at(d_i).at(pt_idx);
          index_t ls_path = ls_paths.at(d_i).at(c_path.k_idx);
          
          // Debugging 
          IntegerVector tmp1 = IntegerVector(ls_path.begin(), ls_path.end()); 
          // Rcout << "LS Path: " << tmp1 << std::endl; 
          
          // Get the source and target level sets
          u8 source_ls = ls_path.at(c_path.c_idx);
          u8 target_ls = ls_path.at(c_path.c_idx+1);
        
          // Extract the (cached) source segment and compute the target segment 
          // u8 source_segment = c_path.c_segment;
          update_ls_segment_idx(i, d_i); // updates which segments the level sets encompass, if need be
          
          IntegerVector tmp = IntegerVector(ls_segment_idx.at(d_i).begin(), ls_segment_idx.at(d_i).end()); 
          
          // std::transform(ls_segment_idx.at(d_i).begin(), ls_segment_idx.at(d_i).end(), tmp.begin(), [](const u8 seg_idx){
          //   return(static_cast<int>(seg_idx));
          // }); 
          // Rcout << "Current segment indices: " << tmp << std::endl; 
          
          const bool intersecting_right = source_ls < target_ls;
          u8 target_segment = intersecting_right ? ls_segment_idx.at(d_i).at(target_ls*2) : ls_segment_idx.at(d_i).at((target_ls*2)+1)-1; 
          
          
          // Update the path information 
          // Rprintf("Updating path index (pt id %d):  %d --> %d\n", pt_idx+1, c_path.c_idx, c_path.c_idx+1);
          c_path.c_idx++;
          // c_path.c_segment = target_segment;
            
          // Debug
          int source_segment = -1;
          Rprintf("pt id %d is going from ls %d to ls %d (going right? %d via segments f=%d, t=%d)\n", 
                  int(pt_idx)+1, int(source_ls), int(target_ls), int(intersecting_right), int(source_segment), int(target_segment));
        
          // If point exists in the update map, update the min/max source/target ls bounds
          // otherwise create a new pt_update and set its current target segment
          std::map< pdiff, pt_update >::iterator pt_it = pts_to_update.lower_bound(pt_idx);
          if (pt_it != pts_to_update.end() && pt_it->first == pt_idx){
            pt_it->second.update_min_max(source_ls, target_ls, d_i);
            pt_it->second.target_segment.at(d_i) = target_segment;
          } else {
            pt_update c_pt = pt_update(pt_idx, d);
            c_pt.update_min_max(source_ls, target_ls, d_i);
            c_pt.target_segment.at(d_i) = target_segment;
            pts_to_update.emplace_hint(pt_it, pt_idx, c_pt);
          }

        }
      } else { // if (!expanding)
        for(pdiff i = c_idx; i >= (t_idx+1); --i){
          pdiff f_i = c_filt_idx.at(i); // filtration index
          pdiff pt_idx = (f_i % n); // the point changing at this step TODO: use a positive-only modulus
          
          // Debugging
          // Rprintf("(%d=%d), pt idx: %d\n", i, f_i, pt_idx);
          
          // Extract the path information
          path_info& c_path = pt_info.at(d_i).at(pt_idx);
          index_t ls_path = ls_paths.at(d_i).at(c_path.k_idx);
          
          // Debugging 
          IntegerVector tmp1 = IntegerVector(ls_path.begin(), ls_path.end()); 
          // Rcout << "LS Path: " << tmp1 << std::endl; 
          
          // Get the source and target level sets
          u8 source_ls = ls_path.at(c_path.c_idx); // where its coming from
          u8 target_ls = ls_path.at(c_path.c_idx-1); // where its going to
          
          // Extract the (cached) source segment and compute the target segment 
          // u8 source_segment = c_path.c_segment;
          update_ls_segment_idx(i, d_i); // updates which segments the level sets encompass, if need be
          
          IntegerVector tmp = IntegerVector(ls_segment_idx.at(d_i).begin(), ls_segment_idx.at(d_i).end()); 
          // Rcout << "Current segment indices: " << tmp << std::endl; 
          
          // TODO: chaneg this?
          const bool intersecting_right = source_ls < target_ls;
          // u8 target_segment = intersecting_right ? ls_segment_idx.at(d_i).at(target_ls*2) : ls_segment_idx.at(d_i).at((target_ls*2)+1)-1; 
          u8 target_segment = intersecting_right ? ls_segment_idx.at(d_i).at((source_ls*2)+1) : ls_segment_idx.at(d_i).at(source_ls*2) -1;  
            
          // Update the path information 
          //Rprintf("Updating path index (pt id %d):  %d --> %d\n", pt_idx+1, c_path.c_idx, c_path.c_idx-1);
          c_path.c_idx--;
          // c_path.c_segment = target_segment;
          
          // Debug
          int source_segment = -1;
          Rprintf("pt id %d is going from ls %d to ls %d (going right? %d via segments f=%d, t=%d)\n", 
                  int(pt_idx)+1, int(source_ls), int(target_ls), int(intersecting_right), int(source_segment), int(target_segment));
          
          std::map< pdiff, pt_update >::iterator pt_it = pts_to_update.lower_bound(pt_idx);
          if (pt_it != pts_to_update.end() && pt_it->first == pt_idx){
            pt_it->second.update_min_max(source_ls, target_ls, d_i);
            pt_it->second.target_segment.at(d_i) = target_segment;
          } else {
            pt_update c_pt = pt_update(pt_idx, d);
            c_pt.update_min_max(source_ls, target_ls, d_i);
            c_pt.target_segment.at(d_i) = target_segment;
            pts_to_update.emplace_hint(pt_it, pt_idx, c_pt);
          }
          
        }
      }
      
      // Update current filtration index
      filtration_idx.at(d_i) = t_idx; 
    } // (auto d_i: d_range)
    

    std::unordered_set< std::size_t > ls_pairs_to_update;
    
    // For all the points that need updating, ensure their default min/max ranges are set, 
    // then iterate through the combinations to create the updated list of LS pairs to generate
    for (auto& kv: pts_to_update){
      const pdiff pt_idx = kv.first;
      pt_update& c_update = kv.second; 
      
      // By default, initialize the min/max ls bounds to the points current LS indices
      for (auto& d_i: d_range){
        path_info& c_path = pt_info.at(d_i).at(pt_idx);
        u8 base_ls = ls_paths.at(d_i).at(c_path.k_idx).at(c_path.p_idx); // TODO: update p_idx
        c_update.update_min_max(base_ls, base_ls, d_i);
      }
      
      // Debugging 
      Rprintf("TO UPDATE: pt_idx=%d: ", pt_idx); 
      for (auto& d_i: d_range){
        Rprintf("d_i=%d, min_ls=%d, max_ls=%d, target_segment=%d\n",
                d_i, 
                c_update.min_ls.at(d_i),
                c_update.max_ls.at(d_i), 
                c_update.target_segment.at(d_i));
      }
      
      // Move the point to its corresponding new multi segment 
      index_t pt_source_segment = extract_segment(pt_idx);
      index_t pt_target_segment = c_update.target_segment;
      l_index_t& from_pts = segment_map.at(pt_source_segment);
      l_index_t::iterator from_end = std::remove_if(from_pts.begin(), from_pts.end(), [pt_idx](const int x_i){ return(x_i == pt_idx); });
      from_pts.resize(std::distance(from_pts.begin(), from_end));
      segment_map.at(pt_target_segment).push_back(pt_idx);
      
     
      // Step 1. Expand the level set indices in each direction
      // std::vector< index_t > ls_expansions(d);
      // for (auto& d_i: d_range){
      //   index_t ls_expansion = seq_ij< u8 >(c_update.min_ls.at(d_i), c_update.max_ls.at(d_i));
      //   ls_expansions.at(d_i) = ls_expansion;
      // }
      
      // Step 2. The cartesian product of the expanded indices comprise the level sets that need to be updated. 
      // Save their corresponding flat indices. 
      // std::vector< std::size_t > flat_ls = std::vector< std::size_t >();
      // CartesianProduct(ls_expansions, [&](const index_t lsmi){
      //   std::size_t lsfi = ls_grid.flat_from_multi(lsmi);
      //   flat_ls.push_back(lsfi);
      // });
      
      // Step 3. The pairwise combinations of the new level sets comprise the LS pairs that need to be recomputed. 
      // Save their corresponding (lower-triangular) flat indices. 
      // const std::size_t n_ls_pairs = ls_grid.n_multi_indices;
      // combine_pairwise(flat_ls.begin(), flat_ls.end(), [&n_ls_pairs, &ls_pairs_to_update](const std::size_t ls_i, const std::size_t ls_j){
      //   std::size_t ij_flat = index_lower_triangular(ls_i, ls_j, n_ls_pairs);
      //   ls_pairs_to_update.insert(ij_flat);
      // });
      
    }
    
    
    return(List::create()); 
  } // update_segments2 
  
  // List update_segments(const IntegerVector target_idx){
  //   
  //   // The information needed to update the filtration
  //   std::vector< index_t > updated_source_segments(d), updated_target_segments(d); 
  //   std::vector< index_t > updated_source_ls(d), updated_target_ls(d);
  //   
  //   // Aggregate per-dimension information
  //   for (std::size_t d_i = 0; d_i < d; ++d_i){
  //     const std::size_t c_idx = filtration_idx.at(d_i); // current filtration index
  //     const std::size_t t_idx = target_idx.at(d_i); // target filration index
  //     const bool expanding = c_idx < t_idx; // are we expanding or contracting the current dimension?
  //     const index_t c_pts = pt_idx.at(d_i); // the current points changing this dimension
  //     
  //     if (expanding){
  //       // Retrieve the unique points in the range ( filtration_idx(d_i), target_idx(d_i) ] and their associated (relative)
  //       // indices into the range. Since the goal is to detect the last unique point positions of the range, use a reverse iterator.  
  //       std::map<int, int> uniq_pts = get_unique_indices(c_pts.rbegin() + (n - t_idx), c_pts.rend()-c_idx-1);
  // 
  //       // Create outer block to capture temporaries
  //       {
  //         int c_pt_id, c_pt_idx; 
  //         u8 c_src_ls, c_target_ls, c_src_segment, c_target_segment;
  //         bool intersecting_right; 
  //         
  //         // Iterate through the points found in the update range
  //         for (auto& pt: uniq_pts){
  //           c_pt_id = pt.first; // point id
  //           c_pt_idx = t_idx - pt.second; // global filtration index
  //           
  //           // Extract the source and target level sets 
  //           c_src_ls = source_ls.at(d_i).at(c_pt_id - 1);
  //           c_target_ls = target_ls.at(d_i).at(c_pt_idx);
  //           
  //           // Use the points filtration index to infer the source/target segments
  //           update_ls_segment_idx(c_pt_idx, d_i); // updates which segments the level sets encompass, if need be
  //           intersecting_right = c_src_ls < c_target_ls;
  //           c_src_segment = source_segment.at(d_i).at(c_pt_id - 1);
  //           c_target_segment = ls_segment_idx.at(d_i).at(intersecting_right ? c_target_ls*2 : c_target_ls*2 + 1);
  //         
  //           // Save the information
  //           updated_source_ls.at(d_i).push_back(c_src_segment);
  //           updated_target_ls.at(d_i).push_back(c_target_segment);
  //           updated_source_segments.at(d_i).push_back(c_src_segment);
  //           updated_target_segments.at(d_i).push_back(c_target_segment);
  //         }
  //       }
  //       
  //       
  //     }
  //     
  //   }
  // 
  //   // get_unique_indices()
  //   return(List::create());
  // }
  // 
};


RCPP_MODULE(multiscale2_module) {
  Rcpp::class_<MultiScale2>("MultiScale2")
  .constructor<const int, IntegerVector>()
  .method( "create_filtration", &MultiScale2::create_filtration )
  .method( "insert_pts", &MultiScale2::insert_pts )
  .method( "create_ls_paths", &MultiScale2::create_ls_paths )
  .method( "set_filtration_rle", &MultiScale2::set_filtration_rle )
  .method( "compute_ls_segment_idx", &MultiScale2::compute_ls_segment_idx )
  .method( "update_segments2", &MultiScale2::update_segments2 )
  ;
}


/*** R
  Rcpp::loadModule("multiscale2_module", TRUE)
*/