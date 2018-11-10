// MultiScale.cpp
// Implementation file for the multiscale indexing structure. 
#include "MultiScale.h"

// using namespace util;

MultiScale::MultiScale(const uidx_t n_pts, IntegerVector resolution) : 
    d(resolution.size()), n(n_pts), 
    num_intervals(resolution.begin(), resolution.end()), 
    ls_grid(GridIndex< uint8_t >(resolution)) 
{
  // Rcout << "Creating new multiscale\n";
  d_range = v_uint8_t(d);
  std::iota(d_range.begin(), d_range.end(), 0);
  
  // Instantiate everything
  pt_info = std::vector< std::vector< path_info > >(d);
  ls_paths = std::vector< std::vector< v_uint8_t > >(d);
  filt_index_set = std::vector< v_sidx_t >(d);
  ls_change_idx = std::vector< v_sidx_t >(d);
  ls_segment_idx = std::vector< v_uint8_t >(d);
  filt_dist = std::vector< NumericVector >(d);
  filt_idx = v_sidx_t(d, -1);
    
  // Assume cover is composed of disjoint sets at initialization 
  c_ls_segment_idx = v_uint8_t(d, -1); 
  for (auto d_i: d_range){ update_ls_segment_idx(0, d_i); }
}
  
SEXP MultiScale::as_XPtr(){
  Rcpp::XPtr< MultiScale > p(this, false); // do not register finalizer
  return(p);
}
  
IntegerMatrix MultiScale::point_info(const uidx_t d_i){
  IntegerMatrix res = IntegerMatrix(n, 4);
  for (uidx_t i = 0; i < n; ++i){
    path_info c_path = pt_info.at(d_i).at(i);
    res(i, _) = IntegerVector::create(int(c_path.k_idx), int(c_path.p_idx), int(c_path.c_idx), int(c_path.c_segment));
  }
  colnames(res) = CharacterVector::create("path_index", "previous_index", "current_index", "current_segment");
  return(res);
}  

IntegerMatrix MultiScale::uniq_paths(const uidx_t d_i){
  const std::vector< v_uint8_t >& paths = ls_paths.at(d_i);
  IntegerMatrix res = IntegerMatrix(paths.size(), num_intervals.at(d_i));
  for (uidx_t i = 0; i < paths.size(); ++i){
    res(i, _) = to_ivec< uint8_t >(paths.at(i));
  }
  return(res);
}
  
IntegerVector MultiScale::extract_level_set(const uidx_t lsfi){
  
  // Step 1. Convert the flat index into a multi index 
  v_uint8_t lsmi = ls_grid.multi_from_flat(lsfi);
  
  // Rcout << "LSFI: " << lsfi << ", LSMI: " << to_ivec<uint8_t>(lsmi) << std::endl; 
  
  // Step 2. Expand the segment indices in the range given by the LS mapping in each direction
  std::vector< v_uint8_t > seg_expansions(d);
  std::size_t cc = 0;
  for (auto& d_i: d_range){
    v_uint8_t c_ls_idx = ls_segment_idx.at(d_i);
    // Rcout << "d_i: " << int(d_i) << ", LS Segment idx: " << to_ivec<uint8_t>(c_ls_idx) << std::endl;
    uint8_t begin_segment = c_ls_idx.at(lsmi.at(d_i)*2);
    uint8_t end_segment = c_ls_idx.at(lsmi.at(d_i)*2 + 1) - 1; // exclusive outer segment 
    v_uint8_t segment_expansion = seq_ij< uint8_t >(begin_segment, end_segment);
    // Rcout << "segments: " << to_ivec< uint8_t >(segment_expansion) << std::endl; 
    seg_expansions.at(d_i) = segment_expansion;
    cc += seg_expansions.size();
  }
  
  // Step 3. The cartesian product of these expansions gives the segments that intersect the given level set 
  std::vector< v_sidx_t* > segments = std::vector< v_sidx_t* >();  
  CartesianProduct(seg_expansions, [&](const v_uint8_t segment){
    if (segment_map.find(segment) != segment_map.end()){
      v_sidx_t& c_segment = segment_map.at(segment);
      if (c_segment.size() > 0){
        // Rcout << "adding segment: " << to_ivec< sidx_t >(c_segment) << std::endl; 
        segments.push_back(&c_segment);
      }
    }
  });
  
  // Step 4. Merge the points at each segment into the level set, and return
  std::vector< sidx_t > res = merge_vectors< sidx_t >(segments);
  return(wrap(res));
  //return(IntegerVector::create());
}
  
List MultiScale::get_segment_map(){
  List res = List();
  for (auto& key_pair: segment_map){
    // std::string key = index_to_str(key_pair.first); 
    std::string key = std::string();
    std::for_each(std::begin(key_pair.first), std::end(key_pair.first), [&key](const uint8_t i){
      int ii = static_cast< int >(i);
      key.append(std::to_string(ii));
    });
    res[key] = wrap(key_pair.second);
  }
  return(res);
}
  
// Expects 1-based RLE built from the point swap indices. Lengths should be cumulative.
void MultiScale::set_filtration_rle(const IntegerVector& ls_changes, const uidx_t d_i){
  ls_change_idx.at(d_i) = v_sidx_t(ls_changes.begin(), ls_changes.end());
}

  
// Expects a 0-based vector of integers representing the total order that point indices 
// along dimension 'd_i' change level sets.
void MultiScale::create_filtration(const IntegerVector& f_idx, const NumericVector& intervals, const uidx_t d_i){
  filt_index_set.at(d_i) = v_sidx_t(f_idx.begin(), f_idx.end());
  filt_dist.at(d_i) = intervals;
}
 
// Given a set of interval sizes (per dimension), retrieves the filtration index corresponding to 
// the largest interval length that is less than the given interval size
IntegerVector MultiScale::get_nearest_filtration_index(NumericVector intervals){
  IntegerVector res(d);
  for (auto& d_i: d_range){
    const NumericVector& R = filt_dist.at(d_i);
    auto it = std::upper_bound(R.begin(), R.end(), intervals.at(d_i));
    res.at(d_i) = std::distance(R.begin(), it)-1; // res allowed to be -1
  }
  return(res);
}
  
v_uint8_t MultiScale::extract_segment(const sidx_t i){
  if (i < 0 || i >= n){ stop("Invalid i"); }
  v_uint8_t res(d);
  for (auto d_i: d_range){ res.at(d_i) = pt_info.at(d_i).at(i).c_segment; }
  return(res);
}

// A contains the 1-based level set indexes each point lies in when the cover is disjoint
// pt_ls_path contains the 1-based index of the unique path taken by 
void MultiScale::insert_pts(const IntegerMatrix& A, const IntegerMatrix& pt_ls_path){
  if (A.nrow() != n || A.ncol() != d){ stop("A must be a matrix of size equal to the number of points in the filt_index_set."); }
  if (pt_ls_path.nrow() != n || pt_ls_path.ncol() != d){ stop("Dimensionality of A must equal dimensionality of filt_index_set."); }
  
  // Clears any points 
  for (auto d_i: d_range){ pt_info.at(d_i).clear();  }
  
  // Add the points to the structure 
  for (std::size_t i = 0; i < n; ++i){
    
    // Assemble the initial path information to each dimension
    path_info new_pt = path_info(); 
    for (std::size_t d_i = 0; d_i < d; ++d_i){
      new_pt.k_idx = static_cast<int>(pt_ls_path(i, d_i)) - 1; // the ls path the point follows along the current dimension
      new_pt.c_idx = new_pt.p_idx = 0; // every point starts at its originating level set
      new_pt.c_segment = (static_cast<int>(A(i, d_i)) - 1) * 2; // every point starts in its originating level sets lower segment
      pt_info.at(d_i).push_back(new_pt);
    }
    
    // Add point to the segment map 
    v_uint8_t c_segment = extract_segment(i);
    segment_map[c_segment].push_back(i);
  }
}
  
// Stores a vector of the distinct ls paths taken by an individual point along dimension d_i. 
void MultiScale::create_ls_paths(const IntegerMatrix& _ls_paths, const uidx_t d_i){
  const std::size_t _n = _ls_paths.nrow(), _k = _ls_paths.ncol(); 
  if (_k != num_intervals.at(d_i)){ stop("Nope"); }
  ls_paths.at(d_i).clear();
  for (std::size_t i = 0; i < _n; ++i){
    IntegerVector tmp = _ls_paths.row(i);
    v_uint8_t c_path = v_uint8_t(tmp.begin(), tmp.end());
    ls_paths.at(d_i).push_back(c_path);
  }
}
  
// Given a global filt_index_set index 'i' \in { -1, 0, 1, ..., n*d_i - 1 } and dimension 'd_i', 
// updates the 'ls_segment_idx' and 'c_ls_segment_idx' variables.
void MultiScale::update_ls_segment_idx(const sidx_t i, const uidx_t d_i){
  const v_sidx_t& cum_lengths = ls_change_idx.at(d_i);
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

v_uint8_t MultiScale::compute_ls_segment_idx(sidx_t i, uidx_t d_i){
  if (d_i >= num_intervals.size()){ stop("Invalid index: Must be 0-based index less than filter dimensionality"); }
  const uidx_t n_level_sets = num_intervals.at(d_i); 
  const uidx_t n_endpts = n_level_sets*2;
  v_uint8_t ls_idx_res = v_uint8_t(n_endpts);
  uidx_t hc = n_level_sets - 1, lc = 0;
  for (uidx_t j = 0; j < n_level_sets; ++j){
    uidx_t j_tmp = n_endpts - j - 1;
    // Rprintf("j: %d, j_tmp: %d, i: %d, lc: %d, hc: %d\n");
    if (j % 2 == 0){
      ls_idx_res[j] = uint8_t(i >= lc ? j - lc : j - i);
      ls_idx_res[j_tmp] = uint8_t(i >= lc ? j_tmp + lc : j_tmp + i);
      ++lc;
    }
    else {
      ls_idx_res[j] = uint8_t(i >= hc ? j + hc : j + i);
      ls_idx_res[j_tmp] = uint8_t(i >= hc ? j_tmp - hc : j_tmp - i);
      --hc;
    }
  }
  return(ls_idx_res);
}
  
List MultiScale::update_segments(const IntegerVector target_idx){
  
  // The information needed to update the filt_index_set
  // std::map< sidx_t, pt_update > pts_to_update = std::map< sidx_t, pt_update >(); 
  // 
  // std::unordered_set< std::size_t > pts_to_update2 = std::unordered_set< std::size_t >(); 
  std::unordered_map< std::size_t, std::vector< uint8_t > > update_list;
  
  std::vector< bool > expansion_status(d, true); // by default, points are expanding 
    
  // Aggregate per-dimension information
  for (auto d_i: d_range){
    const sidx_t c_idx = filt_idx.at(d_i); // current filtration index
    const sidx_t t_idx = target_idx.at(d_i); // target filtration index
    if (c_idx == t_idx){ continue; } // Don't update if we're at the current index
    const bool expanding = c_idx < t_idx; // are we expanding or contracting the current dimension?
    const v_sidx_t c_filt_idx = filt_index_set.at(d_i); // the filt_index_set itself
    expansion_status.at(d_i) = expanding;
      
    // Rprintf("d=%d: expanding? %d. Current index: %d, target index: %d\n", d_i, expanding, c_idx, t_idx);
    sidx_t s_i = expanding ? c_idx+1 : c_idx;
    sidx_t e_i = expanding ? t_idx : t_idx+1;
    auto updating = [&s_i, &e_i, &expanding](){ return(expanding ? s_i <= e_i : s_i >= e_i); };
    auto update = [&s_i, &expanding](){ return(expanding ? ++s_i : --s_i); };
    
    // Update to the target state of the filtration
    for (; updating(); update()){
      sidx_t f_i = c_filt_idx.at(s_i); // filt_index_set index
      std::size_t pt_idx = (f_i % n); // the point changing at this step TODO: use a positive-only modulus
      
      // Extract the path information
      path_info& c_path = pt_info.at(d_i).at(pt_idx);
      v_uint8_t ls_path = ls_paths.at(d_i).at(c_path.k_idx);
      
      // Debugging
      // IntegerVector tmp1 = IntegerVector(ls_path.begin(), ls_path.end());
      // Rprintf("pt_idx %d  path: k=%d, p=%d, c=%d, c_seg=%d\n", pt_idx, c_path.k_idx, c_path.p_idx, c_path.c_idx, c_path.c_segment);
      // Rcout << "LS Path: " << tmp1 << std::endl;
      
      // Get the source and target level sets
      uint8_t source_ls = ls_path.at(c_path.c_idx);
      uint8_t target_ls = ls_path.at(expanding ? c_path.c_idx+1 : c_path.c_idx-1);
      
      // Update the points position 
      c_path.c_idx = expanding ? c_path.c_idx+1 : c_path.c_idx-1;
     
      // Check if the segments comprising the level sets need to be updated
      update_ls_segment_idx(s_i, d_i); 
      
      // IntegerVector tmp = IntegerVector(ls_segment_idx.at(d_i).begin(), ls_segment_idx.at(d_i).end());
      // Rcout << "Current segment indices: " << tmp << std::endl;
      
      // Compute the target segment 
      const bool intersecting_right = source_ls < target_ls;
      uint8_t target_segment = 
        expanding ? 
          (intersecting_right ? ls_segment_idx.at(d_i).at(target_ls*2) : ls_segment_idx.at(d_i).at((target_ls*2)+1)-1) :
          (intersecting_right ? ls_segment_idx.at(d_i).at((source_ls*2)+1) : ls_segment_idx.at(d_i).at(source_ls*2)-1)
      ;
      
      // Insert the point into the list of points to update
      auto element = update_list.find(pt_idx);
      if (element != update_list.end()){
        (*element).second.at(d_i) = target_segment;
      } else {
        update_list[pt_idx] = std::vector< uint8_t >(d, 255); // initialize to -1 to signal the segment hasn't changed
        update_list[pt_idx].at(d_i) = target_segment;
      }
      
      // Debug
      // Rprintf("pt id %d is going from ls %d to ls %d (going right? %d via segments f=%d, t=%d)\n",
      //         int(pt_idx)+1, int(source_ls), int(target_ls), int(intersecting_right), int(c_path.c_segment), int(target_segment));
    } // for(; updating(); update())
    
    filt_idx.at(d_i) = t_idx; // Update current filtration index
  } // for (auto d_i: d_range)
  
  // Update information to collect from the range 
  std::unordered_set< std::size_t > ls_to_update;
  std::unordered_set< std::size_t > ls_pairs_to_update;
  
  // Iterate through the update list, accumulating the level sets and level set pairs that need to be 
  // updated for the current range 
  for(auto& kv: update_list){
    std::size_t pt_idx = kv.first; 
    v_uint8_t c_target_segment = update_list[pt_idx];
    
    // Update target segment dimensions unchanged by the current range 
    for (auto& d_i: d_range){
      path_info& c_path = pt_info.at(d_i).at(pt_idx);
      if (c_target_segment.at(d_i) == 255){ 
        c_target_segment.at(d_i) = c_path.c_segment; 
      }
    }
    
    // Move the point to its new segment
    v_uint8_t pt_source_segment = extract_segment(pt_idx);
    v_sidx_t& from_pts = segment_map.at(pt_source_segment);
    from_pts.erase(
      std::remove_if(from_pts.begin(), from_pts.end(), 
                     [pt_idx](const int x_i){ return(x_i == pt_idx); }), 
                     from_pts.end()
    );
    segment_map[c_target_segment].push_back(pt_idx);
    for (auto& d_i: d_range){ pt_info.at(d_i).at(pt_idx).c_segment = c_target_segment.at(d_i); }
    
    // Rprintf("Pt idx: %d, source segment: %s, target_segment: %s\n", pt_idx,
    //         ls_grid.multi_to_string(pt_source_segment).c_str(),
    //         ls_grid.multi_to_string(c_target_segment).c_str());
    
    
    // Step 1. Expand the level set indices in each direction
    std::vector< v_uint8_t > ls_expansions(d);
    std::size_t cc = 0;
    for (auto& d_i: d_range){
      path_info& c_path = pt_info.at(d_i).at(pt_idx);
      v_uint8_t c_ls_path = ls_paths.at(d_i).at(c_path.k_idx);
      
      // Rprintf("(d_i=%d) pt_idx %d path_info: k=%d, p=%d, c=%d, c_seg=%d, path: %s\n", d_i, pt_idx, c_path.k_idx, c_path.p_idx, c_path.c_idx, c_path.c_segment, ls_grid.multi_to_string(c_ls_path).c_str());
      
      // TODO: if contracting, maybe the begin should be the c_idx, and the end 
      // should be the p_idx, since p_idx > c_idx
      v_uint8_t::iterator begin, end;
      if (expansion_status.at(d_i)){
        begin = c_ls_path.begin(); //  + int(c_path.p_idx);
        end = c_ls_path.begin() + int(c_path.c_idx)+1;
      } else {
        begin = c_ls_path.begin(); // + int(c_path.c_idx);
        end = c_ls_path.begin() + int(c_path.p_idx)+1;
      }

      auto ls_range = std::minmax_element(begin, end);
      // Rprintf("(d_i=%d) = min ls: %d, max ls: %d\n", d_i, *ls_range.first, *ls_range.second);
      sidx_t s = (sidx_t) *ls_range.first, e = (sidx_t) *ls_range.second;
      v_uint8_t ls_expansion = seq_ij< uint8_t >(s, e);
      ls_expansions.at(d_i) = ls_expansion;
      cc += ls_expansion.size();
      
      // Reset the previous idx to the current index
      c_path.p_idx = c_path.c_idx; 
    }
    
    // Only add the updated level sets if the point moved level sets
    if (cc > d){ 
      
      // Step 2. The cartesian product of the expanded indices comprise the level sets that need to be recomputed. 
      // Save their corresponding flat indices. 
      std::vector< std::size_t > flat_ls = std::vector< std::size_t >();
      CartesianProduct(ls_expansions, [&](const v_uint8_t lsmi){
        std::size_t lsfi = ls_grid.flat_from_multi(lsmi);
        flat_ls.push_back(lsfi);
        ls_to_update.insert(lsfi);
      });
      
      // Step 3. The pairwise combinations of the new level sets to update comprise the LS pairs that need to be recomputed. 
      // Save their corresponding (lower-triangular) flat indices. 
      const std::size_t n_ls_pairs = ls_grid.n_multi_indices;
      combine_pairwise(flat_ls.begin(), flat_ls.end(), [&n_ls_pairs, &ls_pairs_to_update](const std::size_t ls_i, const std::size_t ls_j){
        std::size_t ij_flat = index_lower_triangular(ls_i, ls_j, n_ls_pairs);
        ls_pairs_to_update.insert(ij_flat);
      });
    } // if (cc > d){ 
   } // for(auto& pt_idx: update_list)
  
  // Convert the LS to update to an integer vector 
  IntegerVector ls_res = IntegerVector(ls_to_update.begin(), ls_to_update.end());
  
  // Unexpand the the LS pairs to update to an integer matrix 
  IntegerMatrix ls_pairs = no_init_matrix(ls_pairs_to_update.size(), 2);
  const std::size_t n_level_sets = ls_grid.n_multi_indices; 
  std::size_t i = 0; 
  std::for_each(ls_pairs_to_update.begin(), ls_pairs_to_update.end(), [&i, &ls_pairs, &n_level_sets](const std::size_t idx){
    std::size_t to = INDEX_TO(idx, n_level_sets);
    std::size_t from = INDEX_FROM(idx, n_level_sets, to);
    ls_pairs(i++, _) = IntegerVector::create(to, from);
  });
  
  // Return the results 
  return(List::create(_["ls_to_update"] = ls_res, _["ls_pairs_to_update"] = ls_pairs)); 
} // update_segments


RCPP_MODULE(multiscale_module) {
  Rcpp::class_<MultiScale>("MultiScale")
  .constructor<const int, IntegerVector>()
  .field_readonly( "filt_dist", &MultiScale::filt_dist )
  .field_readonly( "filt_idx", &MultiScale::filt_idx )
  .method( "as_XPtr", &MultiScale::as_XPtr ) 
  .method( "point_info", &MultiScale::point_info ) 
  .method( "uniq_paths", &MultiScale::uniq_paths ) 
  .method( "create_filtration", &MultiScale::create_filtration )
  .method( "insert_pts", &MultiScale::insert_pts )
  .method( "create_ls_paths", &MultiScale::create_ls_paths )
  .method( "set_filtration_rle", &MultiScale::set_filtration_rle )
  .method( "compute_ls_segment_idx", &MultiScale::compute_ls_segment_idx )
  .method( "update_segments", &MultiScale::update_segments )
  .method( "get_segment_map", &MultiScale::get_segment_map )
  .method( "get_nearest_filtration_index", &MultiScale::get_nearest_filtration_index )
  .method( "extract_level_set", &MultiScale::extract_level_set )
  ;
}