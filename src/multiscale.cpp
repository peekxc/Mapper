// MultiScale.cpp
// Implementation file for the indexed cover structure. 
// The goal of this structure is to:
// 1) Decompose a symmetric, interval-like cover in disjoint segments 
// 2) Speed up the construction of a cover from these segments, for varying interval sizes
// 3) Given two covers (U_r, U_r'), where r != r' and U_r is already constructed, 
//    return the indices of the preimages in U_r' that changed going from r -> r'. 
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
  // Rprintf("swap idx: %d, c_segment_idx: %d \n", ii, cc_index.at(d_i));
  // Updates the canonical cover based on the current index 'i', if needed
  if (ii != cc_index.at(d_i)){
    canonical_cover.at(d_i) = segment_cover_idx(ii, d_i);
    cc_index.at(d_i) = ii;
  }
}

// Adjusts the state of the fitration to match the given target index. 
// Returns a map between points that have changed and their (multi-index) target segments
std::unordered_map< size_t, vector< uint_t > > MultiScale::update_paths(const IntegerVector target_idx){
  std::unordered_map< size_t, vector< uint_t > > update_list;

  // Adjusts every points path state to match the target index
  for (auto d_i: d_range){
    const int_t c_idx = current_index[d_i]; // current filtration index
    const int_t t_idx = target_idx[d_i]; // target filtration index
    if (c_idx == t_idx){ continue; } // Don't update if we're at the current index
    
    // Identify start and end indices 
    bool expanding = c_idx < t_idx; 
    int_t s_i = expanding ? c_idx+1 : c_idx;
    int_t e_i = expanding ? t_idx : t_idx+1;
    const auto update = [&s_i, &expanding](){ return(expanding ? ++s_i : --s_i); };
    const auto updating = [&s_i, &e_i, &expanding](){ return(expanding ? s_i <= e_i : s_i >= e_i); };
    
    
    // Update to the target state of the filtration
    for (; updating(); update()){
      const int_t f_i = indices[d_i].at(s_i); // current filtration index
      const size_t pt_idx = (f_i % n); // the point changing at this step
      
      // Extract the path information
      path_info& c_path = pt_info[d_i].at(pt_idx); 
      const vector< uint_t > ls_path = ls_paths[d_i].at(c_path.k_idx);
      
      // Get the source and target level sets for the current point
      uint_t next_idx = expanding ? c_path.c_idx+1 : c_path.c_idx-1;
      uint_t source_ls = ls_path.at(c_path.c_idx);
      uint_t target_ls = ls_path.at(next_idx);
      c_path.c_idx = next_idx;  // Update the points position 
      
      // Potentially updates the canonical cover
      update_canonical_cover(s_i, d_i);
      
      // Compute the target segment 
      const bool intersecting_right = source_ls < target_ls;
      uint_t target_segment = 
        expanding ? 
        (intersecting_right ? canonical_cover[d_i].at(target_ls*2) : canonical_cover[d_i].at((target_ls*2)+1)-1) :
        (intersecting_right ? canonical_cover[d_i].at((source_ls*2)+1) : canonical_cover[d_i].at(source_ls*2)-1);
      
      // If the point isn't in the update list, put it in, otherwise update the target segment 
      // for the current point
      auto element = update_list.find(pt_idx);
      if (element != update_list.end()){
        (*element).second[d_i] = target_segment;
      } else {
        // Add the current points segment to the map, then update the target segment index for the current dimension
        update_list[pt_idx] = get_segment_midx(pt_idx); //vector< uint_t >(d, 255); // initialize to 255 to signal the segment hasn't changed
        update_list[pt_idx][d_i] = target_segment;
      }
    } // for(; updating(); update())
    
    current_index[d_i] = t_idx; // Update current index to the target 
  } // for (auto d_i: d_range)
  return update_list; 
}
  
// Given a target index, transfers points to the segments to match the state given by the target index
void MultiScale::update_segments(std::unordered_map< size_t, vector< uint_t > >& update_list){
  
  // Move the points in the updated lists to their new respective segments
  for(auto& pt_seg: update_list){
    const size_t pt_idx = pt_seg.first; 
    const vector< uint_t > pt_target_segment = pt_seg.second;
    
    // Remove point from first segment
    v_int_t& from_pts = segments.at(get_segment_midx(pt_idx));
    erase_remove(from_pts, from_pts.begin(), from_pts.end(), int_t(pt_idx));
    
    // Push point into target segment. Save cached segment index 
    segments[pt_target_segment].push_back(pt_idx);
    for (auto& d_i: d_range){ pt_info[d_i].at(pt_idx).c_segment = pt_target_segment[d_i]; }
  }
}

// Computes the range of relative indices the given point intersects
vector< v_uint_t > MultiScale::resolve_paths(const size_t pt_idx){
  // Lambda which enables adding sequential ranges of LS indices to path_rng
  vector< v_uint_t > path_rng(d);
  size_t cc = 0;
  const auto add_range = [&path_rng, &cc](const vector< uint_t >& ref_rng, const uint_t d_i, const uint_t a, const uint_t b){
    auto rng = seq_ij< uint_t >(a, b); // add entire range
    for (uint_t idx: rng){ path_rng[d_i].push_back(ref_rng.at(idx)); }
    cc += rng.size();
  };
  
  // Use delta between the p_idx and c_idx to determine which open sets each point intersected. 
  for (auto& d_i: d_range){
    path_info& pt_meta = pt_info[d_i].at(pt_idx);
    const vector< uint_t >& c_ls_path = ls_paths[d_i].at(pt_meta.k_idx);
    // 1. If point hasn't moved in this direction, depends on the dimensionality of the space 
    if (pt_meta.p_idx == pt_meta.c_idx){
      if (d_range.size() == 1){
        cc++;
        path_rng[d_i].push_back(c_ls_path.at(pt_meta.c_idx));
      } else {
        add_range(c_ls_path, d_i, 0, pt_meta.c_idx);
      }
    } // 2. Contracting, add [p_idx, c_idx]
    else if (pt_meta.p_idx > pt_meta.c_idx) {
      add_range(c_ls_path, d_i, pt_meta.c_idx, pt_meta.p_idx); 
    } // 3. Expanding, add (p_idx, c_idx]
    else { 
      add_range(c_ls_path, d_i, pt_meta.p_idx+1, pt_meta.c_idx);  
    }
    // Reset the previous idx to the current index
    pt_meta.p_idx = pt_meta.c_idx; 
  } 
  return(path_rng);
}


// Collect the subset of indices in the covers index set have had their decompositions modified 
void MultiScale::modified_indices(
    const std::unordered_map< size_t, vector< uint_t > >& update_list, 
    const size_t max_dim,
    std::unordered_set< size_t >& pullback_idx,
    std::unordered_set< vector< uint_t >, VectorHash< uint_t > >& nerve_idx, 
    const bool all_expanding)
{
  // Number of open sets 
  const size_t n_sets = std::accumulate(begin(num_intervals), end(num_intervals), 1, std::multiplies< size_t >());
  
  // Resolve the changes
  for(auto& pt_seg: update_list){
    const size_t pt_idx = pt_seg.first; 
    const vector< uint_t > c_target_segment = pt_seg.second;
      
    // Get previous and current LS multi-indices
    const vector< uint_t > previous = get_previous_lsmi(pt_idx);
    // const vector< uint_t > target = get_current_lsmi(pt_idx);
  
    // Step 1. Resolve point paths 
    vector< v_uint_t > path_rng = resolve_paths(pt_idx);
  
    // Step 2. The cartesian product of the expanded indices comprise the level sets that
    // need to have their pullbacks recomputed. Save their corresponding flat indices. 
    vector< uint_t > local_pullbacks = vector< uint_t >();
    CartesianProduct(path_rng, [this, &previous, &pullback_idx, &local_pullbacks, &all_expanding](const v_uint_t lsmi){
      if (lsmi != previous || !all_expanding){
        const size_t lsfi = multi_to_flat(lsmi);
        pullback_idx.insert(lsfi);
        local_pullbacks.push_back(lsfi); // local to a point
      }
    });
      
    // Step 3. The pairwise combinations of the new level sets to update comprise the LS pairs that need to be recomputed. 
    v_uint_t init = { uint_t(multi_to_flat(previous)) };
    vector< v_uint_t > local_nerve = { init };
    for (size_t k = 0; k < max_dim; ++k){ 
      local_nerve.push_back(local_pullbacks);
      CartesianProduct(local_nerve, [&nerve_idx](vector< uint_t > k_nerve){
        std::sort(begin(k_nerve), end(k_nerve));
        if (std::unique(begin(k_nerve), end(k_nerve)) == end(k_nerve)){// if all unique elements
          nerve_idx.insert(k_nerve); 
        }
      });
    }
  } // for(auto& pt_idx: update_list)
}

// Updates the index of the cover. Also records changes that occurred. 
// target_idx := target index for the sequence
// max_dim := maximum dimension of simplexes to consider 
List MultiScale::update_index(const IntegerVector target_idx, const size_t max_dim){
  
  // Which dimensions are expanding or contracting?
  vector< bool > expansion_status(d, true); 
  for (auto d_i: d_range){
    expansion_status[d_i] = current_index[d_i] < target_idx[d_i];
  }
  const bool all_expanding = std::all_of(begin(expansion_status), end(expansion_status), identity());
  
  // Move the points into states
  //std::unordered_map< size_t, vector< uint_t > >
  auto update_list = update_paths(target_idx);
  
  // Update the segments
  update_segments(update_list);
  
  // Collect whats changed
  std::unordered_set< size_t > pullback_idx;
  std::unordered_set< vector< uint_t >, VectorHash< uint_t > > nerve_idx;
  modified_indices(update_list, max_dim, pullback_idx, nerve_idx, all_expanding);

  // Return information back to R 
  IntegerVector pullbacks_to_update = IntegerVector(begin(pullback_idx), end(pullback_idx));
  List nerve_to_update = wrap(nerve_idx);
  // ntegerMatrix idx_to_check = no_init_matrix(nerve_idx.size(), max_dim+1);
  // size_t i = 0;
  // for (auto& local_pids: nerve_idx){
  //   IntegerVector tmp(begin(local_pids), end(local_pids));
  //   idx_to_check(i++, _) = tmp;
  // }
  
  // Return the results 
  std::set< size_t > pts_changed; 
  for(auto& pt_seg: update_list){ pts_changed.insert(pt_seg.first); } 
  return(List::create(_["indices_changed"] = pullbacks_to_update, _["pullback_indices_to_check"] = nerve_to_update, _["points_updated"] = wrap(pts_changed))); 
} // update_filtration


RCPP_MODULE(multiscale_module) {
  Rcpp::class_<MultiScale>("MultiScale")
  .constructor<const List>()
  .field_readonly( "indices", &MultiScale::indices )
  .field_readonly( "current_index", &MultiScale::current_index ) // The index of the cover
  .field_readonly( "canonical_cover", &MultiScale::canonical_cover)
  .field_readonly( "canonical_offsets", &MultiScale::cc_offsets)
  .field_readonly( "eps", &MultiScale::eps ) // the interval lengths which cause distinct mappers
  .property( "segments", &MultiScale::get_segments ) // Read-only property to inspect the current segments 
  .property( "segment_table", &MultiScale::get_segment_table )
  .method( "as_XPtr", &MultiScale::as_XPtr ) 
  .method( "point_info", &MultiScale::point_info ) 
  .method( "uniq_paths", &MultiScale::uniq_paths ) 
  .method( "insert_pts", &MultiScale::insert_pts )
  .method( "flat_to_multi", &MultiScale::flat_to_multi )
  .method( "multi_to_flat", &MultiScale::multi_to_flat )
  .method( "update_index", &MultiScale::update_index )
  .method( "segment_cover_idx", &MultiScale::segment_cover_idx )
  .method( "get_nearest_index", &MultiScale::get_nearest_index )
  .method( "extract_level_set", &MultiScale::extract_level_set )
  ;
}
// Info
// vector< size_t > src_seg = vector< size_t >(begin(pt_source_segment), end(pt_source_segment)); 
// vector< size_t > tgt_seg = vector< size_t >(begin(pt_target_segment), end(pt_target_segment)); 
// Rprintf("Pt idx: %d, source segment: %s, target_segment: %s\n", 
//         pt_idx,
//         vec_to_string(src_seg, ",").c_str(),
//         vec_to_string(tgt_seg, ",").c_str()
// );