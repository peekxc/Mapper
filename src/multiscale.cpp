// MultiScale.cpp
// Implementation file for the multiscale indexing structure. 
// The goal of this structure is to:
// 1) Decompose a symmetric, interval-like cover in disjoint segments 
// 2) Speed up the construction of a cover from these segments, for varying interval sizes
// 3) Given two covers (U_e, U_e'), where e != e' and U_e is already constructed, 
//    return which open sets in the U_e would change going from e -> e'. 
#include "MultiScale.h"

template< typename T > 
IntegerVector to_ivec(vector<T> v){
  IntegerVector res(v.size());
  std::transform(begin(v), end(v), res.begin(), [](const T elem) -> int {
    return static_cast< int >(elem);
  });
  return(res);
}

MultiScale::MultiScale(const List params) {
  StringVector arg_names = params.names();
  vector< std::string > args = as< vector< std::string > >(arg_names);
  // for (std::string& st: args){
  //   Rcout << "arg: " << st << std::endl; 
  // }
  n = params["n"];
  // Rcout << n <<  std::endl; 
  IntegerVector ni = params["number_intervals"];
  // Rcout << ni <<  std::endl; 
  num_intervals = v_uint_t(ni.size());
  std::transform(ni.begin(), ni.end(), num_intervals.begin(), [](const int i){
    return( static_cast< uint_t >(i) );
  });
  // std::copy(begin(ni), end(ni), begin(num_intervals));
  
  // for (uint_t& ii: num_intervals){
  //   Rcout << int(ii) << ", ";
  // }
  // Rcout << std::endl;
  d = num_intervals.size();
  // Rcout << d <<  std::endl; 
  
  const IntegerMatrix pt_multi_idx = params["pt_multi_idx"];
  const IntegerMatrix pt_path_idx = params["pt_path_idx"];
  const List target_order = params["target_order"]; 
  const List eps_values = params["eps_values"];
  const List unique_paths = params["unique_paths"]; 
  const List canonical_cuts = params["canonical_cuts"];
  d_range = seq_ij< uint_t >(0, d-1); // useful for iteration
  
  // for (auto nn: num_intervals){   Rcout << nn; }
  // for (uint_t d_i: d_range){   Rcout << "d: " << d_i; }

  // Stores the cumulative product (backwards) to precompute offsets 
  // for faster multi-to-flat index generation
  cum_prod = vector< size_t >(d, 0);
  for (size_t d_i = d; d_i > 0; --d_i){
    cum_prod.at(d_i-1) = std::accumulate(num_intervals.begin(), num_intervals.begin()+d_i-1, 1, std::multiplies<std::size_t>());
  }

  // Instantiate everything
  pt_info = vector< vector< path_info > >(d);
  canonical_cover = vector< v_uint_t >(d);
  current_index = v_int_t(d, -1);
  
  // Inserts points into segments
  segments = std::map< v_uint_t, v_int_t >();
  insert_pts(pt_multi_idx, pt_path_idx);
  
  // Rcout << "d = " << d << std::endl;
  // Store interval sizes and the order of the target sets 
  {
    indices = vector< v_int_t >(d);
    eps = vector< NumericVector >(d);
    for (auto d_i: d_range){
      const IntegerVector to = target_order.at(d_i);
      indices.at(d_i) = v_int_t(begin(to), end(to));
      eps.at(d_i) = eps_values.at(d_i);
    }
  }
  
  // Create unique paths
  {  
    ls_paths = vector< vector< v_uint_t > >(d);
    for (auto d_i: d_range){
      IntegerMatrix u_paths = unique_paths.at(d_i);
      const size_t n_paths = u_paths.nrow(); 
      ls_paths.at(d_i).clear();
      for (size_t i = 0; i < n_paths; ++i){
        IntegerVector tmp = u_paths.row(i);
        v_uint_t c_path = v_uint_t(tmp.begin(), tmp.end());
        ls_paths.at(d_i).push_back(c_path);
      }
    }
  }
  
  // Setup indices where the canonical cover changes
  {
    cc_offsets = vector< v_int_t >(d);
    for (auto d_i: d_range){
      IntegerVector can_offsets = canonical_cuts.at(d_i);
      cc_offsets.at(d_i) = v_int_t(can_offsets.begin(), can_offsets.end());
    }
  }
    
  // Assume cover is composed of disjoint sets at initialization 
  cc_index = v_uint_t(d, -1); 
  for (auto d_i: d_range){ update_canonical_cover(-1, d_i); }
}
  
SEXP MultiScale::as_XPtr(){
  Rcpp::XPtr< MultiScale > p(this, false); // do not register finalizer
  return(p);
}
  
IntegerMatrix MultiScale::point_info(const uint_t d_i){
  IntegerMatrix res = IntegerMatrix(n, 4);
  for (uint_t i = 0; i < n; ++i){
    path_info c_path = pt_info.at(d_i).at(i);
    res(i, _) = IntegerVector::create(int(c_path.k_idx), int(c_path.p_idx), int(c_path.c_idx), int(c_path.c_segment));
  }
  colnames(res) = CharacterVector::create("path_index", "previous_index", "current_index", "current_segment");
  return(res);
}  

IntegerMatrix MultiScale::uniq_paths(const uint_t d_i){
  const vector< v_uint_t >& paths = ls_paths.at(d_i);
  IntegerMatrix res = IntegerMatrix(paths.size(), num_intervals.at(d_i));
  for (uint_t i = 0; i < paths.size(); ++i){
    IntegerVector tmp(begin(paths.at(i)), end(paths.at(i))); 
    res(i, _) = tmp;
  }
  return(res);
}
  
// Given a multi index, returns the equivalent flat index
size_t MultiScale::multi_to_flat(vector< uint_t > index){
  if (index.size() != d){ stop("Given index vector does not match dimensionality."); }
  size_t flat_index = 0, d_i = d;
  // cum_prod = vector< size_t >(d, 0);
  // for (size_t d_i = d; d_i > 0; --d_i){
  //   cum_prod.at(d_i-1) = std::accumulate(num_intervals.begin(), num_intervals.begin()+d_i-1, 1, std::multiplies<std::size_t>());
  // }
  std::for_each(index.rbegin(), index.rend() - 1, [&](uint_t dim_idx){
    flat_index += cum_prod.at(--d_i)*dim_idx;
  });
  flat_index += index.front();
  return(flat_index);
}

// Converts a flat index to its multi-index
// Based on: https://github.com/cran/rje/blob/master/src/rje.c
vector< uint_t > MultiScale::flat_to_multi(size_t idx){
  size_t sz = std::accumulate(begin(num_intervals), end(num_intervals), 1, std::multiplies< uint_t >());
  if (idx < 0 || idx >= sz){ stop("Given index too large, must be between [0, n), where n is the total number of sets."); }
  size_t rem;
  vector< size_t > res(d);
  for (auto d_i: d_range) {
    rem = idx % num_intervals[d_i];
    res[d_i] = rem;
    idx -= rem;
    idx /= num_intervals[d_i];
  }
  vector< uint_t > midx(begin(res), end(res));
  return midx; 
}


// Given a flat index 'lsfi' in [0, 1, ..., n), uses the current segment map to 
// collects the segments whose union forms the level set of the pullback f^{-1}(U_{lsfi}) 
vector< int_t > MultiScale::extract_level_set(const size_t lsfi){
  
  // Step 1. Convert the flat index into a multi index 
  v_uint_t lsmi = flat_to_multi(lsfi);
  // Rcout << "LSFI: " << lsfi << ", LSMI: " << to_ivec<uint_t>(lsmi) << std::endl; 
  
  // Step 2. Expand the segment indices in the range given by the LS mapping in each direction
  vector< v_uint_t > seg_expansions(d);
  size_t cc = 0;
  for (auto& d_i: d_range){
    v_uint_t c_ls_idx = canonical_cover.at(d_i);
    IntegerVector tmp = IntegerVector(c_ls_idx.size());
    std::transform(begin(c_ls_idx), end(c_ls_idx), tmp.begin(), [](uint_t el){
      return static_cast< int >(el); 
    });
    // Rcout << "d_i: " << int(d_i) << ", LS Segment idx: " << tmp << std::endl;
    uint_t begin_segment = c_ls_idx.at(lsmi.at(d_i)*2);
    uint_t end_segment = c_ls_idx.at(lsmi.at(d_i)*2 + 1) - 1; // exclusive outer segment 
    v_uint_t segment_expansion = seq_ij< uint_t >(begin_segment, end_segment);
    // Rcout << "segments: " << to_ivec< uint_t >(segment_expansion) << std::endl; 
    seg_expansions.at(d_i) = segment_expansion;
    cc += seg_expansions.size();
  }
  
  // Step 3. The cartesian product of these expansions gives the segments that intersect the given level set 
  vector< v_int_t* > c_segments = vector< v_int_t* >();  
  CartesianProduct(seg_expansions, [&](const v_uint_t segment){
    if (segments.find(segment) != segments.end()){
      v_int_t& c_segment = segments.at(segment);
      if (c_segment.size() > 0){
        // Rcout << "adding point from segment: " << to_ivec< uint_t >(segment) << std::endl; 
        c_segments.push_back(&c_segment);
      }
    }
  });
  
  // Step 4. Merge the points at each segment into the level set, and return
  vector< int_t > res = merge_vectors< int_t >(c_segments);
  return(res);
}
  
// Converts vector to string
template <typename T>
std::string MultiScale::vec_to_string(vector<T> index, std::string sep){
  std::stringstream result;
  for (const T& item: index){ result << item << sep; }
  std::string key = result.str(); 
  return(key.substr(0, key.length() - 1));
}

// Returns a list of the segments. R-interface only.   
List MultiScale::get_segments(){
  using multi_kv = std::pair< v_uint_t, v_int_t >;
  List res = List();
  for (const multi_kv& kv: segments){
    vector< int > v = vector< int >(begin(kv.first), end(kv.first));
    std::string key = vec_to_string(v, ",");
    res[key] = wrap(kv.second);
  }
  return(res); // remove last sep
}
  
// Expects 1-based RLE built from the point swap indices. Lengths should be cumulative.
// void MultiScale::set_filtration_rle(const IntegerVector& ls_changes, const uint_t d_i){
//   cc_offsets.at(d_i) = v_int_t(ls_changes.begin(), ls_changes.end());
// }

  
// Expects a 0-based vector of integers representing the total order that point indices 
// along dimension 'd_i' change level sets.
// void MultiScale::create_filtration(const IntegerVector& f_idx, const NumericVector& intervals, const uint_t d_i){
//   indices.at(d_i) = v_int_t(f_idx.begin(), f_idx.end());
//   eps.at(d_i) = intervals;
// }

// Given a set of interval sizes (per dimension), retrieves the filtration index corresponding to 
// the largest interval length that is less than the given interval size
IntegerVector MultiScale::get_nearest_index(NumericVector c_eps){
  if (c_eps.size() != d){ stop("Input vector size does not match dimension."); }
  IntegerVector res(d);
  for (auto& d_i: d_range){
    const NumericVector& R = eps.at(d_i);
    auto it = std::upper_bound(R.begin(), R.end(), c_eps.at(d_i));
    res.at(d_i) = std::distance(R.begin(), it)-1; // res allowed to be -1
  }
  return(res);
}

// Returns a list of matrices showing the segment tables per dimension
List MultiScale::get_segment_table(){
  List res = List(d);
  for (auto& d_i: d_range){
    const size_t k = num_intervals.at(d_i);
    IntegerMatrix c_seg_table = IntegerMatrix(k, 2*k);
    vector< size_t > k_rng = seq_ij< size_t >(size_t(0), k-1);
    for (size_t i: k_rng){
      v_uint_t v = segment_cover_idx(i, d_i);
      c_seg_table(i, _) = IntegerVector(begin(v), end(v));
    }
    res.at(d_i) = c_seg_table;
  }
  return(res);
}
  
// Extracts the segment multi index the supplied point 'i' lies in
v_uint_t MultiScale::get_segment_midx(const int_t i){
  if (i < 0 || i >= n){ stop("Invalid i"); }
  v_uint_t res(d);
  for (auto d_i: d_range){ res.at(d_i) = pt_info.at(d_i).at(i).c_segment; }
  return(res);
}

// For a given point index 'idx', retrieves the current level set multi index.
v_uint_t MultiScale::get_current_lsmi(uint_t idx){
  v_uint_t current_lsmi = v_uint_t(d);
  std::transform(begin(d_range), end(d_range), begin(current_lsmi), [this, &idx](const uint_t d_i){
    const path_info& pt_meta = pt_info.at(d_i).at(idx);
    const vector< uint_t >& c_ls_path = ls_paths.at(d_i).at(pt_meta.k_idx);
    return c_ls_path.at(pt_meta.c_idx);
  }); 
  return(current_lsmi);
}

// For a given point index 'idx', retrieves the previous level set multi index.
// The previous LSMI depends on what state the point is in, e.g. whether its 
// expanding or contracting segments.
v_uint_t MultiScale::get_previous_lsmi(uint_t idx){
  v_uint_t prev_lsmi = v_uint_t(d);
  std::transform(begin(d_range), end(d_range), begin(prev_lsmi), [this, &idx](const uint_t d_i){
    const path_info& pt_meta = pt_info.at(d_i).at(idx);
    const vector< uint_t >& c_ls_path = ls_paths.at(d_i).at(pt_meta.k_idx);
    return c_ls_path.at(pt_meta.p_idx);
  }); 
  return(prev_lsmi);
}


// Inserts points into their initial segments
// pt_multi_idx := the point multi-indices 
// pt_path_idx := the index of unique path each pint follows 
void MultiScale::insert_pts(const IntegerMatrix& pt_multi_idx, const IntegerMatrix& pt_path_idx){
  if (pt_multi_idx.nrow() != n || pt_multi_idx.ncol() != d){ stop("'pt_multi_idx' must be a matrix of size equal to the number of points in the indices."); }
  if (pt_path_idx.nrow() != n || pt_path_idx.ncol() != d){ stop("Dimensionality of A must equal dimensionality of indices."); }
  
  // Clears any points 
  for (auto d_i: d_range){ pt_info.at(d_i).clear();  }
  
  // Add the points to the structure 
  for (size_t i = 0; i < n; ++i){
    
    // Assemble the initial path information to each dimension
    path_info new_pt = path_info(); 
    for (size_t d_i = 0; d_i < d; ++d_i){
      new_pt.k_idx = static_cast< uint_t >(pt_path_idx(i, d_i)); // the path the point follows along the current dimension
      new_pt.c_idx = new_pt.p_idx = 0; // every point starts at its originating level set
      new_pt.c_segment = (static_cast< uint_t >(pt_multi_idx(i, d_i))) * 2; // every point starts in its originating level sets lower segment
      pt_info.at(d_i).push_back(new_pt);
    }
    
    // Add point to the segment map 
    v_uint_t c_segment = get_segment_midx(i);
    segments[c_segment].push_back(i);
  }
}
  
// Stores a vector of the distinct ls paths taken by an individual point along dimension d_i. 
// void MultiScale::create_ls_paths(const IntegerMatrix& _ls_paths, const uint_t d_i){
//   const size_t _n = _ls_paths.nrow(), _k = _ls_paths.ncol(); 
//   if (_k != num_intervals.at(d_i)){ stop("Nope"); }
//   ls_paths.at(d_i).clear();
//   for (size_t i = 0; i < _n; ++i){
//     IntegerVector tmp = _ls_paths.row(i);
//     v_uint_t c_path = v_uint_t(tmp.begin(), tmp.end());
//     ls_paths.at(d_i).push_back(c_path);
//   }
// }
  
  
// There is a constant mapping between segments indices and overlapping intervals. 
// Specifically, each interval j is comprised of union of segments in the range [2j, 2j+1). 
// The number of segments each interval requires depends on the canonical cover. The collection 
// of all such maps is what is referred to as the segment table. 
// This function dynamically generates row 'i' in the segment table for dimension d_i. 
v_uint_t MultiScale::segment_cover_idx(int_t i, uint_t d_i){
  if (d_i >= num_intervals.size()){ stop("Invalid index: Must be 0-based index less than filter dimensionality"); }
  const uint_t n_level_sets = num_intervals.at(d_i); 
  const uint_t n_endpts = n_level_sets*2;
  v_uint_t ls_idx_res = v_uint_t(n_endpts);
  if (i < 0){ i = 0; }
  uint_t hc = n_level_sets - 1, lc = 0;
  for (uint_t j = 0; j < n_level_sets; ++j){
    uint_t j_tmp = n_endpts - j - 1;
    // Rprintf("j: %d, j_tmp: %d, i: %d, lc: %d, hc: %d\n");
    if (j % 2 == 0){
      ls_idx_res[j] = uint_t(i >= lc ? j - lc : j - i);
      ls_idx_res[j_tmp] = uint_t(i >= lc ? j_tmp + lc : j_tmp + i);
      ++lc;
    }
    else {
      ls_idx_res[j] = uint_t(i >= hc ? j + hc : j + i);
      ls_idx_res[j_tmp] = uint_t(i >= hc ? j_tmp - hc : j_tmp - i);
      --hc;
    }
  }
  return(ls_idx_res);
}  
  
// Given a global index 'i' \in { -1, 0, 1, ..., n*(k-1)-1 } and dimension 'd_i', 
// updates the 'canonical_cover' and 'cc_index' variables.
void MultiScale::update_canonical_cover(const int_t i, const uint_t d_i){
  size_t ii;
  if (i  == -1){ ii = 0; } 
  else {
    const auto& cc_off = cc_offsets.at(d_i); 
    auto ub = std::upper_bound(begin(cc_off), end(cc_off), i);
    ii = std::distance(begin(cc_off), ub);
  }
  // if (ii > 0){ ii -= 1; }
  // Rprintf("swap idx: %d, c_segment_idx: %d \n", ii, cc_index.at(d_i));
  // Updates the canonical cover based on the current index 'i', if needed
  if (ii != cc_index.at(d_i)){
    canonical_cover.at(d_i) = segment_cover_idx(ii, d_i);
    cc_index.at(d_i) = ii;
  }
  // Rcout << "here_extra" << std::endl;
}
  
List MultiScale::update_segments(const IntegerVector target_idx){
  
  // The information needed to update the indices
  // std::map< int_t, pt_update > pts_to_update = std::map< int_t, pt_update >(); 
  // 
  // std::unordered_set< size_t > pts_to_update2 = std::unordered_set< size_t >(); 
  std::unordered_map< size_t, vector< uint_t > > update_list;
  
  vector< bool > expansion_status(d, true); // by default, points are expanding 
    
  // Aggregate per-dimension information
  for (auto d_i: d_range){
    // Locals 
    const int_t c_idx = current_index.at(d_i); // current filtration index
    const int_t t_idx = target_idx.at(d_i); // target filtration index
    if (c_idx == t_idx){ continue; } // Don't update if we're at the current index
    const bool expanding = c_idx < t_idx; // are we expanding or contracting the current dimension?
    expansion_status.at(d_i) = expanding;
      
    // Rprintf("d=%d: expanding? %d. Current index: %d, target index: %d\n", d_i, expanding, c_idx, t_idx);
    // Adjust start and end indices
    int_t s_i = expanding ? c_idx+1 : c_idx;
    int_t e_i = expanding ? t_idx : t_idx+1;
    
    // Update lambda updates the start index, 
    // Updating lambda checks whether the start index as arrived at the end
    auto update = [&s_i, &expanding](){ return(expanding ? ++s_i : --s_i); };
    auto updating = [&s_i, &e_i, &expanding](){ return(expanding ? s_i <= e_i : s_i >= e_i); };

    // Update to the target state of the filtration
    for (; updating(); update()){
      int_t f_i = indices.at(d_i).at(s_i); // current filtration index
      size_t pt_idx = (f_i % n); // the point changing at this step
      
      // Extract the path information
      path_info& c_path = pt_info.at(d_i).at(pt_idx);
      v_uint_t ls_path = ls_paths.at(d_i).at(c_path.k_idx);
      
      // Debugging
      // IntegerVector tmp1 = IntegerVector(ls_path.begin(), ls_path.end());
      // Rprintf("pt_idx %d  path: k=%d, p=%d, c=%d, c_seg=%d\n", pt_idx, c_path.k_idx, c_path.p_idx, c_path.c_idx, c_path.c_segment);
      // Rcout << "LS Path: " << tmp1 << std::endl;
      
      // Get the source and target level sets for the current point
      uint_t next_idx = expanding ? c_path.c_idx+1 : c_path.c_idx-1;
      uint_t source_ls = ls_path.at(c_path.c_idx);
      uint_t target_ls = ls_path.at(next_idx);
      c_path.c_idx = next_idx;  // Update the points position 
     
      // Potentially updates the canonical cover
      update_canonical_cover(s_i, d_i);
      
      // IntegerVector tmp = IntegerVector(canonical_cover.at(d_i).begin(), canonical_cover.at(d_i).end());
      // Rcout << "Current segment indices: " << tmp << std::endl;
      
      // Compute the target segment 
      const bool intersecting_right = source_ls < target_ls;
      uint_t target_segment = 
        expanding ? 
          (intersecting_right ? canonical_cover.at(d_i).at(target_ls*2) : canonical_cover.at(d_i).at((target_ls*2)+1)-1) :
          (intersecting_right ? canonical_cover.at(d_i).at((source_ls*2)+1) : canonical_cover.at(d_i).at(source_ls*2)-1)
      ;
      
      // If the point isn't in the update list, put it in, otherwise update 
      // the target segment for the current point
      auto element = update_list.find(pt_idx);
      if (element != update_list.end()){
        (*element).second.at(d_i) = target_segment;
      } else {
        // Add the current points segment to the map, then update the target segment index for the current dimension
        update_list[pt_idx] = get_segment_midx(pt_idx); //vector< uint_t >(d, 255); // initialize to 255 to signal the segment hasn't changed
        update_list[pt_idx].at(d_i) = target_segment;
      }
      
      // Once the point has moved, check if the canonical cover needs to be updated
      // update_canonical_cover(s_i, d_i); 
      
      // Debug
      // Rprintf("pt id %d is going from ls %d to ls %d (going right? %d via segments f=%d, t=%d)\n",
      //         int(pt_idx)+1, int(source_ls), int(target_ls), int(intersecting_right), int(c_path.c_segment), int(target_segment));
    } // for(; updating(); update())
    
    current_index.at(d_i) = t_idx; // Update current filtration index
  } // for (auto d_i: d_range)
  
  // If any dimensions  are contracting, keep p_idx --> c_idx expansions conservative 
  const bool all_expanding = std::all_of(begin(expansion_status), end(expansion_status), identity());
  // Rcout << "all expanding: " << all_expanding << std::endl;
    
  // Move the points in the updated lists to their new respective segments
  for(auto& pt_seg: update_list){
    size_t pt_idx = pt_seg.first; 
    v_uint_t pt_target_segment = pt_seg.second;
    
    // Remove point from first segment
    v_uint_t pt_source_segment = get_segment_midx(pt_idx);
    v_int_t& from_pts = segments.at(pt_source_segment);
    erase_remove(from_pts, from_pts.begin(), from_pts.end(), int_t(pt_idx));
    
    // Push point into target segment. Save cached segment index 
    segments[pt_target_segment].push_back(pt_idx);
    for (auto& d_i: d_range){ pt_info.at(d_i).at(pt_idx).c_segment = pt_target_segment.at(d_i); }
      
    // Info
    // vector< size_t > src_seg = vector< size_t >(begin(pt_source_segment), end(pt_source_segment)); 
    // vector< size_t > tgt_seg = vector< size_t >(begin(pt_target_segment), end(pt_target_segment)); 
    // Rprintf("Pt idx: %d, source segment: %s, target_segment: %s\n", 
    //         pt_idx,
    //         vec_to_string(src_seg, ",").c_str(),
    //         vec_to_string(tgt_seg, ",").c_str()
    // );
  }
   
  // Collect which pullback sets need updating
  std::unordered_set< size_t > flat_ls;
  std::set< std::pair< size_t, size_t > > flat_ls_pairs;
  
  const size_t n_sets = std::accumulate(begin(num_intervals), end(num_intervals), 1, std::multiplies< size_t >());
  const size_t n_ls_pairs = choose< size_t >(n_sets, 2);
  
  for(auto& pt_seg: update_list){
    size_t pt_idx = pt_seg.first; 
    v_uint_t c_target_segment = pt_seg.second;
      
    // Step 1. Expand the level set indices in each direction
    vector< v_uint_t > path_rng(d);
    size_t cc = 0;
    
    // Collect the source set the point is coming from to remove from the cartesian product
    v_uint_t previous = get_previous_lsmi(pt_idx);
    v_uint_t target = get_current_lsmi(pt_idx);
    IntegerVector tmp1 = IntegerVector(begin(previous), end(previous));
    IntegerVector tmp2 = IntegerVector(begin(target), end(target));
    // Rprintf("Source/target lsmi for pt %d: ", pt_idx);
    // Rcout << tmp1 << ", ";
    // Rcout << tmp2 << std::endl;
    
    // The 'path_info' uses the cartesian product of the delta between the p_idx and c_idx to  
    // determine which sets were intersected. 
    for (auto& d_i: d_range){
      path_info& pt_meta = pt_info.at(d_i).at(pt_idx);
      const vector< uint_t >& c_ls_path = ls_paths.at(d_i).at(pt_meta.k_idx);
      
      // Rprintf("(d_i=%d) pt_idx %d path_info: k=%d, p=%d, c=%d, c_seg=%d, path: %s\n", d_i, pt_idx, pt_meta.k_idx, pt_meta.p_idx, pt_meta.c_idx, pt_meta.c_segment, vec_to_string(c_ls_path, ", ").c_str());
      
      // TODO: if contracting, maybe the begin should be the c_idx, and the end 
      // should be the p_idx, since p_idx > c_idx
      // auto path_begin = c_ls_path.begin();
      // auto path_end = c_ls_path.begin() + (expansion_status.at(d_i) ? int(pt_meta.c_idx)+1 : int(pt_meta.p_idx)+1);
      // 
      // auto idx_range = std::minmax_element(path_begin, path_end);
      // const int_t s = (int_t) *idx_range.first, e = (int_t) *idx_range.second;
      // path_rng.at(d_i) = seq_ij< uint_t >(s, e);
      // cc += path_rng.at(d_i).size();
      
      // Point hasn't moved in this direction 
      if (pt_meta.p_idx == pt_meta.c_idx){
        if (d_range.size() == 1){
          cc++;
          path_rng.at(d_i).push_back(c_ls_path.at(pt_meta.c_idx));
        } else {
          auto c_path_idx = seq_ij< uint_t >(0, pt_meta.c_idx); // add entire range
          for (uint_t idx: c_path_idx){
            path_rng.at(d_i).push_back(c_ls_path.at(idx));
          }
          cc += c_path_idx.size();
        }

      }
      // Contracting, add everything between
      else if (pt_meta.p_idx > pt_meta.c_idx){
        auto c_path_idx = seq_ij< uint_t >(pt_meta.c_idx, pt_meta.p_idx);
        for (uint_t idx: c_path_idx){
          path_rng.at(d_i).push_back(c_ls_path.at(idx));
        }
        cc += c_path_idx.size();
      }
      // Expanding, depends on dimensionality 
      else {
        // If 1-dimensional, just need to add up to and including the current index, but not the previous index
        // if (d_range.size() == 1){
          auto c_path_idx = seq_ij< uint_t >(pt_meta.p_idx+1, pt_meta.c_idx);
          for (uint_t idx: c_path_idx){
            path_rng.at(d_i).push_back(c_ls_path.at(idx));
          }
          cc += c_path_idx.size();
        // } 
        // Otherwise, need to add the previous index to collect the pullbacks in the product 
        // else {
        //   auto c_path_idx = seq_ij< uint_t >(pt_meta.p_idx, pt_meta.c_idx);
        //   for (uint_t idx: c_path_idx){
        //     path_rng.at(d_i).push_back(c_ls_path.at(idx));
        //   }
        //   cc += c_path_idx.size();
        // }
      }
      // Reset the previous idx to the current index
      pt_meta.p_idx = pt_meta.c_idx; 
    } 
    
    // Only add the updated level sets if the point moved level sets
    // if (cc > d){ 
      
    // Step 2. The cartesian product of the expanded indices comprise the level sets that
    // need to have their pullbacks recomputed. Save their corresponding flat indices. 
    vector< uint_t > local_flat_ls = vector< uint_t >();
    // TODO: only add target level set! This is not
    CartesianProduct(path_rng, [this, &previous, &target, &flat_ls, &local_flat_ls, &all_expanding](const v_uint_t lsmi){
      // Rcout << "lsmi: " << int(multi_to_flat(lsmi)) << std::endl;
      // If any dimensions are contracting, insert the index. If they're all expanding, only insert
      // indices not originating from the points source pullback
      if (lsmi != previous || !all_expanding){
        size_t lsfi = multi_to_flat(lsmi);
        flat_ls.insert(lsfi); // global 
        local_flat_ls.push_back(lsfi); // local to a point
      }
    });
      
    // Step 3. The pairwise combinations of the new level sets to update comprise the LS pairs that need to be recomputed. 
    // Save their corresponding (lower-triangular) flat indices. 
    
    // TODO: Change this structure, report instead 
    // 1) The source set whose pullback changed, as a key in a map 
    // 2) The newly intersecting target sets, as values in a std::set mapped by the source set 
    // Then in R, take the (k-1)-fold combinations of the values w/ the source set. 
    // Add a pair between the source set and any other sets the point intersected
    vector< uint_t > src_vec = { uint_t(multi_to_flat(previous)) };
    vector< v_uint_t > pair_sets = { src_vec, local_flat_ls};
    // for (size_t k = 0; k < max_dim; ++k){ pair_sets.push_back(local_flat_ls); }
    CartesianProduct(pair_sets, [&flat_ls_pairs, &n_sets](const v_uint_t pair_lsfi){
      // size_t ij_flat = index_lt(, n_sets);
      if (pair_lsfi.at(0) != pair_lsfi.at(1)){
        size_t i = (size_t) std::min(pair_lsfi.at(0), pair_lsfi.at(1));
        size_t j = (size_t) std::max(pair_lsfi.at(0), pair_lsfi.at(1));
        flat_ls_pairs.insert(std::make_pair(i, j));
      }
    });
        
      // combine_pairwise(flat_ls.begin(), flat_ls.end(), [&n_sets, &flat_ls_pairs](const size_t ls_i, const size_t ls_j){
      //   size_t ij_flat = index_lt(ls_i, ls_j, n_sets);
      //   flat_ls_pairs.insert(ij_flat);
      // });
    // } // if (cc > d){ 
   } // for(auto& pt_idx: update_list)
  
  // Convert the LS to update to an integer vector 
  // IntegerVector ls_res = IntegerVector(ls_to_update.begin(), ls_to_update.end());
  
  // // Unexpand the the LS pairs to update to an integer matrix 
  // const size_t n_sets = std::accumulate(begin(num_intervals), end(num_intervals), 1, std::multiplies< size_t >());
  IntegerMatrix ls_pairs = no_init_matrix(flat_ls_pairs.size(), 2);
  size_t i = 0;
  for (auto idx: flat_ls_pairs){
    // size_t to = index_to(idx, n_ls_pairs);
    // size_t from = index_from(idx, n_ls_pairs, to);
    ls_pairs(i++, _) = IntegerVector::create(idx.first, idx.second);
  }

  
  // Return the results 
  //, _["ls_pairs_to_update"] = ls_pairs
  std::set< size_t > pts_changed; 
  for(auto& pt_seg: update_list){ pts_changed.insert(pt_seg.first); } 
  return(List::create(_["indices_changed"] = flat_ls, _["pairs_changed"] = ls_pairs, _["points_updated"] = wrap(pts_changed))); 
} // update_segments

// Updates the simplex tree, and the corresponding vertices
// TODO: Rewrite as *update_pullback*
// void MultiScale::update_vertices(const IntegerVector which_levels, const NumericMatrix& X, const Function f, List& ls_vertex_map, SEXP stree){
//   // Update the vertices in the level sets dynamically
//   for (const int index: which_levels){
//     const IntegerVector level_set = extract_level_set(index) + 1; // convert to 1-based for R
//     update_level_set(index, level_set, X, f, vertices, ls_vertex_map, stree);
//   }
// }

// Converts a given list of vertices to a map and assigns to the structure
// void MultiScale::initialize_vertices(List& ls_vertex_map, List& vertices){
//   this->vertices = vertices_to_map(ls_vertex_map, vertices);
// }

// Given a pullback cover consisting of level sets, 
// void MultiScale::initialize(const List& pullback){
//   
// }

// void MultiScale::persistence(){
//   
// }


// This
// void dist_to_boundary(const vector< double >& x, const vector< double >& lb, const vector< double >& ub){
//   
// }


RCPP_MODULE(multiscale_module) {
  Rcpp::class_<MultiScale>("MultiScale")
  .constructor<const List>()
  .field_readonly( "indices", &MultiScale::indices )
  .field_readonly( "current_index", &MultiScale::current_index ) // The index of the cover
  .field_readonly( "canonical_cover", &MultiScale::canonical_cover)
  .field_readonly( "canonical_offsets", &MultiScale::cc_offsets)
  .field_readonly( "eps", &MultiScale::eps ) // the interval lengths which cause distinct mappers

  .field_readonly( "vertices", &MultiScale::vertices )
  .property( "segments", &MultiScale::get_segments ) // Read-only property to inspect the current segments 
  .property( "segment_table", &MultiScale::get_segment_table )
  .method( "as_XPtr", &MultiScale::as_XPtr ) 
  .method( "point_info", &MultiScale::point_info ) 
  .method( "uniq_paths", &MultiScale::uniq_paths ) 
  
  // .method( "create_filtration", &MultiScale::create_filtration )
  // .method( "create_ls_paths", &MultiScale::create_ls_paths )
  // .method( "set_filtration_rle", &MultiScale::set_filtration_rle )
  
  .method( "insert_pts", &MultiScale::insert_pts )
  
  .method( "flat_to_multi", &MultiScale::flat_to_multi )
  .method( "multi_to_flat", &MultiScale::multi_to_flat )

  .method( "segment_cover_idx", &MultiScale::segment_cover_idx )
  .method( "update_segments", &MultiScale::update_segments )

  .method( "get_nearest_index", &MultiScale::get_nearest_index )
  .method( "extract_level_set", &MultiScale::extract_level_set )
  // .method( "update_vertices", &MultiScale::update_vertices )
  // .method( "initialize_vertices", &MultiScale::initialize_vertices )
  ;
}