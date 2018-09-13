#include <Rcpp.h>
using namespace Rcpp;

#include <array>
#include <cstdint>
#include <stack>
#include <set>
#include "MultiSegmentTree.h"

typedef std::vector< uint_fast8_t > index_t;
typedef std::vector< std::shared_ptr<uint_fast8_t> > index_ptr_t;
typedef std::ptrdiff_t s_size_t;

std::string index_to_str(index_t idx){
  std::stringstream result;
  std::copy(idx.begin(), idx.end(), std::ostream_iterator<int>(result, " "));
  return(result.str());
}

uint_fast8_t make_index_t(const int idx){
  return(static_cast<uint_fast8_t>(idx));
} 

std::shared_ptr<uint_fast8_t> make_index_ptr_t(const int idx){
  uint_fast8_t idx_val = static_cast<uint_fast8_t>(idx);
  std::shared_ptr<uint_fast8_t> idx_ptr = std::shared_ptr<uint_fast8_t>(new uint_fast8_t { idx_val });
  return(idx_ptr);
}

#define SWAP(x, y, T) do { T SWAP = x; x = y; y = SWAP; } while (0)

template<typename ForwardIterator>
std::map<int, int> get_unique_indices(ForwardIterator first, ForwardIterator last){
  std::map<int, int> pt_to_unique_idx;
  for(std::size_t i = 0; first != last; ++i, ++first){
    auto it = pt_to_unique_idx.find(*first);
    if (it == pt_to_unique_idx.end()) { // value doesn't exist
      pt_to_unique_idx.emplace(*first, i);
    }
  }
  return pt_to_unique_idx;
}

// Computes the absolute distances from a given point to the closest endpoint of each level set. This distance should 
// represent 1/2 the smallest interval length the target level set would have to be (via expansion) to intersect the given point.
// [[Rcpp::export]]
List dist_to_boxes(const IntegerVector& positions, const double interval_length, const int num_intervals, const NumericVector& dist_to_lower, const NumericVector& dist_to_upper) {
  
  // Sequence from 1 - < number of intervals >  
  std::vector<int> all_positions = std::vector<int>(num_intervals);
  std::iota(std::begin(all_positions), std::end(all_positions), 1);
  for (int j = 0; j < all_positions.size(); ++j){
    Rcout << all_positions.at(j) << ", ";
  }
  Rcout << std::endl; 
  
  // Iterators 
  IntegerVector::const_iterator pos_it = positions.begin();
  NumericVector::const_iterator dtl_it = dist_to_lower.begin(), dtu_it = dist_to_upper.begin();
  //IntegerVector::iterator k; 

  // Variables needed 
  const int n = positions.size();
  // int current_position[] = { 0 };
  std::array<int, 1> current_position = { {0} };
  
  // To fill each iteration
  IntegerVector target_positions = no_init(num_intervals - 1);
  //std::vector<int> target_positions = std::vector<int>(num_intervals - 1);
  NumericVector target_distances = no_init(num_intervals - 1);
  
  // Outputs 
  IntegerMatrix res_pos = no_init_matrix(n, num_intervals - 1);
  NumericMatrix res_dist = no_init_matrix(n, num_intervals - 1);
  
  double dtl = 0.0, dtu = 0.0;
  for (int i = 0, pos = 0; i < n; ++i, ++pos_it, ++dtl_it, ++dtu_it){
    pos = *pos_it, dtl = *dtl_it, dtu = *dtu_it;
    current_position[0] = pos;
    Rcout << current_position[0] << std::endl; 
    std::set_difference(all_positions.begin(), all_positions.end(), current_position.begin(), current_position.end(), target_positions.begin());
    
    for (int j = 0; j < target_positions.size(); ++j){
      Rcout << target_positions.at(j) << ", ";
    }
    Rcout << std::endl;
    Rcout << pos << ", " << dtl << ", " << dtu << std::endl; 
    // if (target_positions.size() != (num_intervals - 1)){ 
    //   stop("Something went wrong. ");
    // }
    // auto dist_to_box_f = [pos, dtl, dtu]  () {};
    std::transform(target_positions.begin(), target_positions.end(), target_distances.begin(), 
     [interval_length, num_intervals, pos, dtl, dtu](int target_position){
      if (target_position < pos){ return(dtl + (pos - target_position - 1) * interval_length); }
      else { return(dtu + (target_position - pos - 1) * interval_length); }
     });
    res_dist.row(i) = clone(target_distances);
    res_pos.row(i) = clone(target_positions);
  }
  return(List::create(_["target_pos"] = res_pos, _["target_dist"] = res_dist));
  // return(List::create());
}

struct MultiScale {
  const NumericMatrix& data;
  const Function clustering_function;
  const int d;
  int n;
  // Rcpp::XPtr<SegmentTree> tr; // the segment tree contains the level sets 
  // std::vector<LS_idx> level_set_idx;
  // std::vector< std::vector<uint> > level_sets;
  std::vector< IntegerVector > pt_idx, from_ls_idx, to_ls_idx, swap_idx; 
  std::vector< NumericVector > dist_to_box, swap_dist;
  // std::vector< uint_fast8_t > ls_idx;
  // std::stack< uint_fast8_t, std::vector< uint_fast8_t > > ls_idx_swaps;
  std::vector< s_size_t > current_index, current_cover; 
  IntegerVector resolution; 
  std::shared_ptr<MultiSegmentTree> mtree; 
  std::vector< index_ptr_t > ls_idx; // tracks the per dimension endpt configuration
  std::map< index_t, uint_fast8_t > ls_map;
  std::map< index_ptr_t, uint_fast8_t > endpt_to_ls_map;
  std::vector< std::vector< uint_fast8_t > > pt_segment_idx; 
  
  std::map< index_t, std::vector< int > > segment_map; 
  // std::vector< IntegerVector > segment_pts;
  
  MultiScale(const NumericMatrix& X, const Function f, const IntegerVector& k) : data(X), clustering_function(f), d(k.size()) { //, tr(stree) {
    bool res_check = std::any_of(resolution.begin(), resolution.end(), [](const int k_i){ return(k_i > 255); });
    if (res_check){ stop("Multiscale version of mapper only supports up to 255 intervals per dimension."); }
    //level_set_idx = std::vector<LS_idx>();
    current_index = std::vector< s_size_t >(d, -1); // points to the index of the filtration that has not been computed. 
    current_cover = std::vector< s_size_t >(d, 0);
    resolution = k;
    pt_idx = std::vector< IntegerVector >(d);
    from_ls_idx = std::vector< IntegerVector >(d);
    to_ls_idx = std::vector< IntegerVector >(d);
    swap_idx = std::vector< IntegerVector >(d);
    ls_idx = std::vector< index_ptr_t >(d);
    dist_to_box = std::vector< NumericVector >(d);
    swap_dist = std::vector< NumericVector >(d);
    pt_segment_idx = std::vector< std::vector< uint_fast8_t > >(d);
    
    // Multi-dimensional testing
    segment_map = std::map< index_t, std::vector< int > >(); 
    // Rcout << resolution << std::endl;
    
    
    // Prepare vector of dimension indices to use for each transform
    std::vector<uint_fast8_t> dim_idx(d);
    std::iota(dim_idx.begin(), dim_idx.end(), 0);
    
    // First, make the base level map, where each level set is disjoint. This populates the initial set of 
    // shared pointers to the endpoint vectors. 
    std::for_each(dim_idx.begin(), dim_idx.end(), [&](const int_fast8_t d_i){
      index_t base_ls_idx = compute_ls_idx(0, d_i);
      index_ptr_t idx_key = index_ptr_t(base_ls_idx.size());
      std::transform(base_ls_idx.begin(), base_ls_idx.end(), idx_key.begin(), [](const uint_fast8_t idx){
        return(make_index_ptr_t(idx));
      });
      ls_idx.at(d_i) = idx_key;
    });
    
    // TODO: inspect these
    make_ls_map();
    make_endpt_ls_map();
  }
  
  void set_everything(const IntegerVector& _pt_idx, 
                      const IntegerVector& _from_ls_idx, 
                      const IntegerVector& _to_ls_idx, 
                      const IntegerVector& _swap_idx, 
                      const NumericVector& _dist_to_box, 
                      const NumericVector& _swap_dist, 
                      const int d_i) {
    pt_idx.at(d_i) = _pt_idx;
    from_ls_idx.at(d_i) = _from_ls_idx; 
    to_ls_idx.at(d_i) = _to_ls_idx; 
    swap_idx.at(d_i) = _swap_idx;
    dist_to_box.at(d_i) = _dist_to_box;
    swap_dist.at(d_i) = _swap_dist;
  }
  
  // void create_flat_sets(const NumericMatrix& filter_values, const List& endpts){
  //   if (filter_values.ncol() != d){ stop("Filter points dimensionality does not match resolution dimension."); }
  //   
  //   // Load R functions
  //   Function findInterval = getFunctionR("findInterval"), rle = getFunctionR("rle"), order = getFunctionR("order");
  //   
  //   // Preprocessing
  //   List results = List(d);
  //   for (int d_i = 0; d_i < d; ++d_i){
  //     NumericVector pts = filter_values.column(d_i);
  //     NumericVector endpoints = endpts.at(d_i);
  //       
  //     // Extratc the order of the sorted endpoints
  //     NumericVector s_endpts = clone(endpoints).sort(); // store the (sorted) endpoints 
  //     IntegerVector o_endpts = match(s_endpts, endpoints); // save the original ordering for interval queries
  //     
  //     // Sort the point / retrieve the original ordering
  //     IntegerVector o_pts = order(pts); // sorted ordering
  //     NumericVector s_pts = pts[o_pts - 1]; // sort coordinates
  //     
  //     // Partition the points into the given intervals
  //     IntegerVector pt_int = findInterval(_["x"] = s_pts, _["vec"] = s_endpts, _["rightmost.closed"] = false, _["all.inside"] = true);
  //     
  //     // Use run-length encoding to get the interval sizes
  //     const List rle_pt_int = rle(pt_int);
  //     const IntegerVector int_len = rle_pt_int["lengths"];
  //     const IntegerVector int_val = rle_pt_int["values"];
  //   
  //     // Save everything
  //     results.at(d_i) = List::create(
  //       _["int_val"] = int_val, 
  //       _["int_len"] = int_len, 
  //       _["s_pts"] = s_pts, 
  //       _["o_pts"] = o_pts, 
  //       _["s_endpts"] = s_endpts,
  //       _["o_endpts"] = o_endpts, 
  //       _["endpoints"] = endpoints, 
  //       _["pts"] = pts
  //     );
  //   }
  //   
  //   segment_pts = std::vector< IntegerVector >();
  //   const std::size_t n_segments = std::accumulate(resolution.begin(), resolution.end(), 1, std::multiplies<int>())*2; 
  //   segment_pts.reserve(n_segments);
  //   for (int i = 0; i < n_segments; ++i){
  //     ls_pts.at(i);
  //   }
  //   
  //   // Assign the point indices to the leaves. The indices are sorted on copy.
  //   const int ne_int = int_len.size(); // number of non-empty intervals
  //   for (int i = 0, ci = 0; i < ne_int; ++i){
  //     const int int_index = int_val.at(i) - 1, interval_len = int_len.at(i);
  //     leaves.at(int_index).resize(interval_len);
  //     std::partial_sort_copy(o_pts.begin() + ci, o_pts.begin() + ci + interval_len,
  //                            leaves.at(int_index).begin(), leaves.at(int_index).end());
  //     ci += interval_len;
  //   }
  //   // std::vector<uint_fast8_t> dim_idx(d); 
  //   // std::iota(dim_idx.begin(), dim_idx.end(), 0);
  //   // for (const auto &key_pair : endpt_to_ls_map){
  //   //   IntegerVector query_l, query_r;
  //   //   std::for_each(dim_idx.begin(), dim_idx.end(), [&](const int_fast8_t d_i){
  //   //     const index_ptr_t& idx =  key_pair.first;
  //   //     query_l.push_back(*idx.at(2*d_i));  
  //   //     query_r.push_back((*idx.at(2*d_i + 1)) + 1);  
  //   //   });
  //   //   ls_pts.at(static_cast<int>(key_pair.second)) = mtree->query(query_l, query_r);
  //   // }
  // }
  
  // List get_flat_sets(){
  //   return(wrap(ls_pts));
  // }
  
  // Generates the level set multi-indexes as a vector of integer vectors
  std::vector< IntegerVector > get_cart_prod(std::vector< IntegerVector > elems){
    Function expand_grid = getFunctionR("expand.grid");
    List cp_args = List();
    for (int i = 0; i < elems.size(); ++i){
      std::string key("d"); 
      key += std::to_string(i);
      cp_args[key] = elems.at(i);
    }
    DataFrame ls_multi_idx = expand_grid(cp_args);
    
    // Convert the indices to vector of indices 
    std::vector< IntegerVector > indices =  std::vector< IntegerVector >();
    for (int i = 0; i < elems.size(); ++i){
      std::string key = "d" + std::to_string(i);
      indices.push_back( ls_multi_idx[key] );
    }
    return(indices);
  }
  
  // Generates the level set multi-indexes as a vector of integer vectors
  std::vector< IntegerVector > get_multi_indexes(){
    Function expand_grid = getFunctionR("expand.grid");
    List cp_args = List();
    for (int d_i = 0; d_i < d; ++d_i){
      std::vector<int> current_ls_idx = std::vector<int>(resolution.at(d_i));
      std::iota(current_ls_idx.begin(), current_ls_idx.end(), 1);
      std::string key("d"); 
      key += std::to_string(d_i);
      cp_args[key] = current_ls_idx;
    }
    DataFrame ls_multi_idx = expand_grid(cp_args);
    
    // Convert the indices to vector of indices 
    std::vector< IntegerVector > indices =  std::vector< IntegerVector >();
    for (int d_i = 0; d_i < d; ++d_i){
      std::string key = "d" + std::to_string(d_i);
      indices.push_back( ls_multi_idx[key] );
    }
    return(indices);
  }
  
  // Generates a mapping from an integer vector key to a flat level set index. The keys in the map are assumed to follow
  // cartesian product indexing scheme, and the values are assumed to be sequential from {0, 1, ..., n} where n is the number 
  // of level sets in the map. 
  void make_ls_map(){
    std::vector< IntegerVector > ls_multi_idx = get_multi_indexes();
    ls_map = std::map< index_t, uint_fast8_t >();
    
    // Prepare vector of dimension indices to use for each transform
    std::vector<uint_fast8_t> dim_idx(d); 
    std::iota(dim_idx.begin(), dim_idx.end(), 0);
  
    // Collect the endpoint indices for each level set 
    const uint_fast8_t n_level_sets = ls_multi_idx.at(0).size();
    for (uint_fast8_t i = 0; i < n_level_sets; ++i){
      index_t idx_key = index_t(d);
      std::for_each(dim_idx.begin(), dim_idx.end(), [&, i](const int_fast8_t d_i){
        const uint_fast8_t idx = (uint_fast8_t) ls_multi_idx.at(d_i).at(i) - 1;
        idx_key.at(d_i) = idx;
      });
      ls_map.emplace(idx_key, i);
    }
  }
  
  // Generates a mapping from an integer vector key to a flat index. The keys in the map are assumed to track a 
  // cartesian product indexing schems, and the values are assumed to be sequential from {0, 1, ..., n}^d where n is the number 
  // of end points in the map per dimension d_i in {1, ..., d}. 
  void make_endpt_ls_map(){
    std::vector< IntegerVector > ls_multi_idx = get_multi_indexes();
    endpt_to_ls_map = std::map< index_ptr_t, uint_fast8_t >();
    
    // Prepare vector of dimension indices to use for each transform
    std::vector<uint_fast8_t> dim_idx(d);
    std::iota(dim_idx.begin(), dim_idx.end(), 0);
    
    // Collect the endpoint indices for each level set
    const uint_fast8_t n_level_sets = ls_multi_idx.at(0).size();
    for (uint_fast8_t i = 0; i < n_level_sets; ++i){
      index_ptr_t segment_idx = index_ptr_t();
      std::for_each(dim_idx.begin(), dim_idx.end(), [&, i](const int d_i){
        const uint_fast8_t idx = (uint_fast8_t) ls_multi_idx.at(d_i).at(i) - 1;
        const std::shared_ptr<uint_fast8_t> seg_start = ls_idx.at(d_i).at(idx*2);
        const std::shared_ptr<uint_fast8_t> seg_end = ls_idx.at(d_i).at((idx*2) + 1);
        segment_idx.push_back(seg_start);
        segment_idx.push_back(seg_end);
      });
      endpt_to_ls_map.emplace(segment_idx, i);
    }
  }
  
  List get_endpt_ls_map(){
    std::vector< IntegerVector > ls_multi_idx = get_multi_indexes();
    endpt_to_ls_map = std::map< index_ptr_t, uint_fast8_t >();
    
    // Prepare vector of dimension indices to use for each transform
    std::vector<uint_fast8_t> dim_idx(d);
    std::iota(dim_idx.begin(), dim_idx.end(), 0);
    
    // Collect the endpoint indices for each level set
    const uint_fast8_t n_level_sets = ls_multi_idx.at(0).size();
    List res = List();
    for (uint_fast8_t i = 0; i < n_level_sets; ++i){
      IntegerVector segment_idx = IntegerVector();
      std::string key = "( ";
      std::for_each(dim_idx.begin(), dim_idx.end(), [&, i](const int d_i){
        const uint_fast8_t idx = (uint_fast8_t) ls_multi_idx.at(d_i).at(i) - 1;
        const uint_fast8_t seg_start = *ls_idx.at(d_i).at(idx*2);
        const uint_fast8_t seg_end = *ls_idx.at(d_i).at((idx*2) + 1);
        segment_idx.push_back(static_cast<int>(seg_start));
        segment_idx.push_back(static_cast<int>(seg_end));
        key += std::to_string(idx + 1) + (d_i == (d - 1) ? " " : ", ");
      });
      key += ")";
      res[key] = segment_idx;
    }
    return(res);
  }
  
  int query_ls_map(const IntegerVector& query_idx){
    if (query_idx.size() != d){ stop("Illegal query dimension."); }
    bool res_check = std::any_of(query_idx.begin(), query_idx.end(), [](const int k_i){ return(k_i > 255); });
    if (res_check){ stop("Illegal query."); }
    index_t dim_idx(d); 
    std::copy(query_idx.begin(), query_idx.end(), dim_idx.begin());
    return(static_cast<int>(ls_map.at(dim_idx)));
  }
  
  void build_segment_trees(const List& endpts, const NumericMatrix& filter_pts){
    mtree = std::shared_ptr<MultiSegmentTree>(new MultiSegmentTree(endpts));
    mtree->insert_pts(filter_pts);
    
    // Setup up cached index to store which segment the points fall in
    const std::size_t n = filter_pts.nrow();
    Rcout << "n: " << n << std::endl; 
    for (int d_i = 0; d_i < d; ++d_i){
      std::vector< uint_fast8_t > tmp = std::vector< uint_fast8_t >(n, -1);
      uint_fast8_t i = 0; 
      const std::shared_ptr<SegmentTree>& ctree = mtree->tr.at(d_i);
      std::vector< std::vector<int> >::iterator leaf_it;
      for (leaf_it = ctree->leaves.begin(); leaf_it != ctree->leaves.end(); ++leaf_it, ++i){
        // For each leaf, assign to values of tmp at the corresponding point indices the current segment index i
        std::for_each((*leaf_it).begin(), (*leaf_it).end(), [&tmp, &i](const int pt){
          tmp.at(pt - 1) = i;
        });
      }
      pt_segment_idx.at(d_i) = tmp;
    }
  }
  
  List get_segment_map(){
    List res = List();
    for (auto& key_pair: segment_map){
      std::string key = index_to_str(key_pair.first); 
      res[key] = wrap(key_pair.second);
    }
    return(res);
  }
  
  void build_multiscale_configuration(const IntegerMatrix& A, const NumericMatrix& filter_pts){
    if (A.nrow() != filter_pts.nrow() || A.ncol() != filter_pts.ncol()){
      stop("Dimensionality of inputs do not match.");
    }
    const int n = A.nrow();
    index_t key = index_t(d);
    for (std::size_t i = 0; i < n; ++i){
      IntegerMatrix::ConstRow pt_idx = A.row(i);
      std::transform(pt_idx.begin(), pt_idx.end(), key.begin(), make_index_t);
      auto it = segment_map.lower_bound(key);
      if (it != segment_map.end() && it->first == key) {
        it->second.push_back(static_cast<int>(i));
      } else {
        std::vector< int > v = { static_cast<int>(i) };
        segment_map.emplace_hint(it, key, v);
      }
    }
    
    
    // if (filter_pts.ncol() != d){ stop("Filter points dimensionality does not match resolution dimension."); }
    // Function findInterval = getFunctionR("findInterval"), rle = getFunctionR("rle"), order = getFunctionR("order");
    // List tmp = List();
    // std::vector< IntegerVector > int_len_v = std::vector< IntegerVector >();
    // std::vector< IntegerVector > int_val_v = std::vector< IntegerVector >();
    // std::vector< IntegerVector > o_pts_v = std::vector< IntegerVector >();
    // for (int d_i = 0; d_i < d; ++d_i){
    //   NumericVector pts = filter_pts.column(d_i);
    //   NumericVector endpoints = endpts.at(d_i);
    // 
    //   // Extratc the order of the sorted endpoints
    //   NumericVector s_endpts = clone(endpoints).sort(); // store the (sorted) endpoints
    //   IntegerVector o_endpts = match(s_endpts, endpoints); // save the original ordering for interval queries
    //   
    //   // Sort the point / retrieve the original ordering
    //   IntegerVector o_pts = order(pts); // sorted ordering
    //   NumericVector s_pts = pts[o_pts - 1]; // sort coordinates
    //   
    //   // Partition the points into the given intervals
    //   IntegerVector pt_int = findInterval(_["x"] = s_pts, _["vec"] = s_endpts, _["rightmost.closed"] = false, _["all.inside"] = true);
    //   
    //   // Use run-length encoding to get the interval sizes
    //   const List rle_pt_int = rle(pt_int);
    //   const IntegerVector int_len = rle_pt_int["lengths"];
    //   const IntegerVector int_val = rle_pt_int["values"];
    //   
    //   // Store the needed intermediate result prior to building the map
    //   // tmp.at(d_i) = List::create(_["int_len"] = int_len, _["int_val"] = int_val, _["o_pts"] = o_pts);
    //   int_len_v.push_back(int_len);
    //   int_val_v.push_back(int_val);
    //   o_pts_v.push_back(o_pts);
    // }
    
    // std::vector< IntegerVector > segment_idx = std::vector< IntegerVector >();
    // std::transform(tmp.begin(), tmp.end(), segment_idx.begin(), [](const List& l_i){
    //   as<IntegerVector>(l_i["int_val"]);
    // });
    // Assign the point indices to the leaves. The indices are sorted on copy.
    // const int ne_int = int_len.size(); // number of non-empty intervals
    // for (int i = 0, ci = 0; i < ne_int; ++i){
    //   const int int_index = int_val.at(i) - 1, interval_len = int_len.at(i);
    //   leaves.at(int_index).resize(interval_len);
    //   std::partial_sort_copy(o_pts.begin() + ci, o_pts.begin() + ci + interval_len,
    //                          leaves.at(int_index).begin(), leaves.at(int_index).end());
    //   ci += interval_len;
    // }
    // LogicalVector
    
    // std::vector< IntegerVector > segment_multi_idx = get_cart_prod(int_val_v);
    // const int n_idx = segment_multi_idx.at(0).size();
    // std::vector<uint_fast8_t> dim_idx(d), multi_key(d);
    // std::iota(dim_idx.begin(), dim_idx.end(), 0);
    // for (int i = 0; i < n_idx; ++i){
    //   std::transform(dim_idx.begin(), dim_idx.end(), multi_key.begin(), [&](const uint_fast8_t d_i){
    //     segment_multi_idx.at(i).at(d_i);
    //   });
    //   std::vector<int> pt_idx = std::vector<int>();
    //   std::for_each(dim_idx.begin(), dim_idx.end(), [&](const uint_fast8_t d_i){
    //     int_
    //   });
    //   // segment_map.emplace(multi_key, );
    // }
    // 
    // 
    //   for (const auto &key_pair : endpt_to_ls_map){
    //     IntegerVector query_l, query_r;
    //     std::for_each(dim_idx.begin(), dim_idx.end(), [&](const int_fast8_t d_i){
    //       const index_ptr_t& idx =  key_pair.first;
    //       query_l.push_back(*idx.at(2*d_i));
    //       query_r.push_back((*idx.at(2*d_i + 1)) + 1);
    //     });
    //     ls_pts.at(static_cast<int>(key_pair.second)) = mtree->query(query_l, query_r);
    //   }
  }
  
  // inline void swap_ls_index(std::size_t from, std::size_t to){
  //   if (from > level_set_idx.size() || to > level_set_idx.size()){ stop("Indexing error."); }
  //   if (from < to){
  //     uint_fast8_t tmp = level_set_idx[from].endpt2; 
  //     level_set_idx[from].endpt2 = level_set_idx[to].endpt1;
  //     level_set_idx[to].endpt1 = tmp;
  //     // level_set_idx[from].endpt2 += level_set_idx[to].endpt1;
  //     // level_set_idx[to].endpt1 = level_set_idx[from].endpt2 - level_set_idx[to].endpt1;
  //     // level_set_idx[from].endpt2 -= level_set_idx[to].endpt1;
  //   } else {
  //     uint_fast8_t tmp = level_set_idx[to].endpt1; 
  //     level_set_idx[to].endpt2 = level_set_idx[from].endpt1;
  //     level_set_idx[from].endpt1 = tmp;
  //   }
  // }

  
  // inline void swap_ls(std::size_t i){
  //   if (i == 0 || i >= (ls_idx.size()/2)) { stop("Invalid swap configuration."); }
  //   const std::size_t n_swaps = (ls_idx.size() - i*2)/2;
  //   for (std::size_t c = 0; c < n_swaps; c++, i+=2){
  //     std::swap(ls_idx[i], ls_idx[i + 1]);
  //   }
  // }
  
  
  // Returns the level sets at a given cover parameterization
  // Available data: pt_idx, from_ls_idx, to_ls_idx, swap_idx;
  // dist_to_box, swap_dist;
  void update_cover(const s_size_t target_index, const s_size_t d_i){
    if (target_index < -1 || target_index >= pt_idx.at(d_i).size()){ stop("'next_index' must be non-negative and less than the size of the filtration."); }
    // Rcout << "current index: " << current_index.at(d_i) << std::endl;
    // If we're already at the solution, return
    if (current_index.at(d_i) == target_index){ return; }
    
    // Otherwise, check which direction we're going
    const bool is_increasing = current_index.at(d_i) < target_index; 
    s_size_t i; 
    auto update_index = [&i, is_increasing](){ is_increasing ? ++i : --i; };
    
    const IntegerVector& c_pt_idx = pt_idx.at(d_i);
    const IntegerVector& c_swap_idx = swap_idx.at(d_i);
    const IntegerVector& c_from_ls = from_ls_idx.at(d_i);
    const IntegerVector& c_to_ls = to_ls_idx.at(d_i);
    std::shared_ptr<SegmentTree>& c_tree = mtree->tr.at(d_i);
    
    // Loop through the parameterization
    if (is_increasing){
      for (i = current_index.at(d_i) + 1; i <= target_index; update_index()){
        // Rcout << "index: " << i << std::endl;
        if (current_cover.at(d_i) != c_swap_idx.at(i)){ // Update the level set segment index
          current_cover.at(d_i) = c_swap_idx.at(i);
          index_t new_ls_idx = compute_ls_idx(c_swap_idx.at(i), d_i);
          assign_ls_idx(new_ls_idx, ls_idx.at(d_i)); // dynamically update the indices
        }
        
        IntegerVector current_ls_config = get_ls_idx(c_swap_idx.at(i), d_i);
        // Rcout << "Current ls configuration: " << current_ls_config << std::endl; 
        
        // Transfer points
        const int from_ls = c_from_ls.at(i);
        const int to_ls = c_to_ls.at(i);
        // Rprintf("Searching for pt %d between leaf indices [%d, %d)\n", c_pt_idx.at(i), current_ls_config.at(from_ls*2), current_ls_config.at((from_ls*2) + 1));
        // 
        // TODO: fix this with either a vector mapping pt index -> segment id, or via 
        // some distance calculation that takes the current threshold into account.
        int from_leaf = mtree->tr.at(d_i)->find_leaf(c_pt_idx.at(i), current_ls_config.at(from_ls*2), current_ls_config.at((from_ls*2) + 1)); 
        int to_leaf = from_ls < to_ls ? current_ls_config.at(to_ls*2) : current_ls_config.at(to_ls*2 + 1) - 1; // exclusive outer
        // Rprintf("Swapping point %d from leaf %d (in ls %d) to leaf %d (in ls %d)\n", c_pt_idx.at(i), from_leaf, from_ls, to_leaf, to_ls);
        
        // Swap in the segemnt trees
        c_tree->swap(from_leaf, to_leaf, c_pt_idx.at(i));
      }
      // Final update, since the current index is exclusive 
      //update_index();
      current_index.at(d_i) = target_index; // Finish up 
      // Rcout << "final index: " << current_index.at(d_i) << std::endl; 
    } else {
      for(i = current_index.at(d_i); i > target_index; update_index()){
        // Rcout << "index: " << i << std::endl;
        if (current_cover.at(d_i) != c_swap_idx.at(i)){ // Update the level set segment index
          current_cover.at(d_i) = c_swap_idx.at(i);
          index_t new_ls_idx = compute_ls_idx(c_swap_idx.at(i), d_i);
          assign_ls_idx(new_ls_idx, ls_idx.at(d_i)); // dynamically update the indices
        }
        IntegerVector current_ls_config = get_ls_idx(c_swap_idx.at(i), d_i);
        
        // Transfer points
        const int from_ls = c_to_ls.at(i);
        const int to_ls = c_from_ls.at(i);
        //Rprintf("Searching for pt %d between leaf indices [%d, %d)\n", c_pt_idx.at(i), current_ls_config.at(from_ls*2), current_ls_config.at((from_ls*2) + 1));
        
        // TODO: fix this with either a vector mapping pt index -> segment id, or via 
        // some distance calculation that takes the current threshold into account.
        int from_leaf = mtree->tr.at(d_i)->find_leaf(c_pt_idx.at(i), current_ls_config.at(from_ls*2), current_ls_config.at((from_ls*2) + 1)); 
        int to_leaf = to_ls < from_ls ? from_leaf - 1 : from_leaf + 1;
        // Rprintf("Swapping point %d from leaf %d (in ls %d) to leaf %d (in ls %d)\n", c_pt_idx.at(i), from_leaf, from_ls, to_leaf, to_ls);
        
        // TODO: improve this! remove segment tree and replace w/ index mapped list updated by reference! 
        c_tree->swap(from_leaf, to_leaf, c_pt_idx.at(i));
      }
      // Final update, since the current index is exclusive 
      //update_index();
      current_index.at(d_i) = target_index; // Finish up 
      // Rcout << "final index: " << current_index.at(d_i) << std::endl; 
    }
    
    // If we updated all the way to start, set cover to original disjoint cover
    if (current_index.at(d_i) == -1){
      current_cover.at(d_i) = 0;
      index_t new_ls_idx = compute_ls_idx(0, d_i);
      assign_ls_idx(new_ls_idx, ls_idx.at(d_i)); // dynamically update the indices
    }
  }
  
  void update_multi_cover(const IntegerVector target_index){
    if (target_index.size() != d){ stop("Invalid query, doesn't match dimensionality of filter space."); }
    int d_i = 0; 
    bool all_equal = std::all_of(target_index.begin(), target_index.end(), [&](const int ti){ return(ti == current_index.at(d_i++)); });
    if (all_equal) { return; }
    std::vector< std::vector<int> > pts_to_change_per_d(d); 
    std::vector< std::vector<uint_fast8_t> > pt_seg_from(d); 
    std::vector< std::vector<uint_fast8_t> > pt_seg_to(d); 
    std::vector< std::vector<std::size_t> > pt_absolute_idx(d); 
    std::unordered_set<int> pts_to_change;
    for (d_i = 0; d_i < d; ++d_i){
      const bool is_increasing = current_index.at(d_i) < target_index.at(d_i); // Otherwise, check which direction we're going
      const int start_idx = is_increasing ? current_index.at(d_i)+1 : current_index.at(d_i);
      const int end_idx = is_increasing ? target_index.at(d_i) : target_index.at(d_i) + 1;
      const IntegerVector& c_swap_idx = swap_idx.at(d_i);
      const IntegerVector& c_pt_idx = pt_idx.at(d_i);
      const IntegerVector& c_from_ls = from_ls_idx.at(d_i);
      const IntegerVector& c_to_ls = to_ls_idx.at(d_i);
      std::vector< int > pt_idx_vec = std::vector<int>(c_pt_idx.begin() + start_idx, c_pt_idx.begin() + end_idx);
      // IntegerVector pt_ids = wrap(pt_idx_vec);
      // Rcout << pt_ids << std::endl; 
      std::map<int, int> pt_uniq_idx = get_unique_indices(pt_idx_vec.rbegin(), pt_idx_vec.rend());
      for (auto& pidx: pt_uniq_idx){
        const int idx = pt_idx_vec.size() - pidx.second - 1; // relative index
        // Rprintf("pt id: %d, relative idx: %d, start idx: %d, end_idx: %d\n", pidx.first, idx, start_idx, end_idx);
        index_t new_ls_idx = compute_ls_idx(c_swap_idx.at(start_idx + idx), d_i);
        const int from_ls = c_from_ls.at(idx), to_ls = c_to_ls.at(idx);
        const int pt = pidx.first;
        const int from_segment = pt_segment_idx.at(d_i).at(pt - 1);
        const int to_segment = from_ls < to_ls ? new_ls_idx.at(to_ls*2) : new_ls_idx.at(to_ls*2 + 1) - 1; // exclusive outer
        Rprintf("pt id: %d (rel idx=%d) going from ls %d to %d (segment %d to %d)\n", pidx.first, idx, from_ls, to_ls, from_segment, to_segment);
        pts_to_change_per_d.at(d_i).push_back(pt);
        pt_seg_from.at(d_i).push_back(from_segment);
        pt_seg_to.at(d_i).push_back(to_segment);
        pts_to_change.insert(pt);
        pt_absolute_idx.at(d_i).push_back(start_idx + idx);
      } // given the unique pts and their indices, get their corresponding segment locations
    }
    Rcout << "finsihg preprocessing " << std::endl; 
    // Prepare vector of dimension indices
    std::vector<int> dim_idx(d);
    std::iota(dim_idx.begin(), dim_idx.end(), 0);
    
    // Swap the points the between the appropriate segments!  
    index_t key_from(d), key_to(d); 
    std::for_each(pts_to_change.begin(), pts_to_change.end(), [&](const int pt){
      std::for_each(dim_idx.begin(), dim_idx.end(), [&](const int d_i){
        std::vector<int> candidate_pts = pts_to_change_per_d.at(d_i);
        ptrdiff_t pos = std::distance(candidate_pts.begin(), std::find(candidate_pts.begin(), candidate_pts.end(), pt));
        if (pos == candidate_pts.size()) { // doesn't exist
          Rprintf("Using point %d's (dim=%d) default current positon %d\n", pt, d_i, pt-1);
          key_from.at(d_i) = pt_segment_idx.at(d_i).at(pt-1); // use wherever it's at 
          key_to.at(d_i) = pt_segment_idx.at(d_i).at(pt-1); // use wherever it's at 
        } else {
          key_from.at(d_i) = pt_seg_from.at(d_i).at(pos);
          key_to.at(d_i) = pt_seg_to.at(d_i).at(pos);
        }
      });
      std::string key_str_from = index_to_str(key_from);
      std::string key_str_to = index_to_str(key_to);
      Rprintf("pt: %d from %s to %s\n", pt, key_str_from.c_str(), key_str_to.c_str());
    });

    
    
    // s_size_t i; 
    // std::vector< index_t > from_segments(d), to_segments(d);
    // for (d_i = 0; d_i < d; ++d_i){
    //   const bool is_increasing = current_index.at(d_i) < target_index.at(d_i); // Otherwise, check which direction we're going
    //   auto update_index = [&i, is_increasing](){ is_increasing ? ++i : --i; };
    //   const IntegerVector& c_pt_idx = pt_idx.at(d_i);
    //   const IntegerVector& c_swap_idx = swap_idx.at(d_i);
    //   const IntegerVector& c_from_ls = from_ls_idx.at(d_i);
    //   const IntegerVector& c_to_ls = to_ls_idx.at(d_i);
    //   // std::shared_ptr<SegmentTree>& c_tree = mtree->tr.at(d_i);
    //   if (is_increasing){
    //     for (i = current_index.at(d_i) + 1; i <= target_index.at(d_i); update_index()){
    //       // Rcout << "index: " << i << std::endl;
    //       if (current_cover.at(d_i) != c_swap_idx.at(i)){ // Update the level set segment index
    //         current_cover.at(d_i) = c_swap_idx.at(i);
    //         index_t new_ls_idx = compute_ls_idx(c_swap_idx.at(i), d_i);
    //         assign_ls_idx(new_ls_idx, ls_idx.at(d_i)); // dynamically update the indices
    //       }
    //       IntegerVector current_ls_config = get_ls_idx(c_swap_idx.at(i), d_i);
    //       const int from_ls = c_from_ls.at(i), to_ls = c_to_ls.at(i);
    //       from_segments.at(d_i).push_back(pt_segment_idx.at(d_i).at(c_pt_idx.at(i)));
    //       to_segments.at(d_i).push_back(from_ls < to_ls ? current_ls_config.at(to_ls*2) : current_ls_config.at(to_ls*2 + 1) - 1); // exclusive outer
    //     }
    //   } else {
    //     for(i = current_index.at(d_i); i > target_index.at(d_i); update_index()){
    //       if (current_cover.at(d_i) != c_swap_idx.at(i)){ // Update the level set segment index
    //         current_cover.at(d_i) = c_swap_idx.at(i);
    //         index_t new_ls_idx = compute_ls_idx(c_swap_idx.at(i), d_i);
    //         assign_ls_idx(new_ls_idx, ls_idx.at(d_i)); // dynamically update the indices
    //       }
    //       IntegerVector current_ls_config = get_ls_idx(c_swap_idx.at(i), d_i);
    //       const int from_ls = c_to_ls.at(i), to_ls = c_from_ls.at(i);
    //       from_segments.at(d_i).push_back(pt_segment_idx.at(d_i).at(c_pt_idx.at(i)));
    //       to_segments.at(d_i).push_back(to_ls < from_ls ? from_segments.at(d_i).back() - 1 : from_segments.at(d_i).back() + 1);
    //     }
    //   }
    // }
    // 

    // 
    // // Update the segments appropriately
    // std::vector<int> from_update, to_update;
    // std::unique_copy(from_segments.at(d_i).rbegin(), from_segments.at(d_i).rend(), std::back_inserter(from_update));
    // std::unique_copy(to_segments.at(d_i).rbegin(), to_segments.at(d_i).rend(), std::back_inserter(to_update));
    // // std::map<int, int> sti_from, sti_to;
    // // std::for_each(dim_idx.begin(), dim_idx.end(), [&](const int d_i){
    // //   sti_from = get_unique_indices(from_segments.at(d_i).rbegin(), from_segments.at(d_i).rend());
    // //   sti_to = get_unique_indices(to_segments.at(d_i).rbegin(), to_segments.at(d_i).rend());
    // //   for (auto& key_pair: sti_from){
    // //     int segment_from = key_pair.first
    // //     int pt_idx = 
    // //   }
    // // });
    // 
    // // std::find_first_of(pt_idx.)
    // // c_tree->swap(from_leaf, to_leaf, c_pt_idx.at(i));
    // for (auto& key_pair: segment_map){
    //   std::string key = index_to_str(key_pair.first); 
    //   res[key] = wrap(key_pair.second);
    // }
    // segment_map.
    // 
    // // Finish up
    // std::copy(target_index.begin(), target_index.end(), current_index.begin());
    // 
    // // If we updated all the way to start, set cover to original disjoint cover
    // for (int d_i = 0; d_i < d; ++d_i){
    //   if (current_index.at(d_i) == -1){
    //     current_cover.at(d_i) = 0;
    //     index_t new_ls_idx = compute_ls_idx(0, d_i);
    //     assign_ls_idx(new_ls_idx, ls_idx.at(d_i)); // dynamically update the indices
    //   }
    // }
  } // update_multi_cover 
  
  IntegerVector get_segment_idx(const int d_i){
    return wrap(pt_segment_idx.at(d_i));
  }
  
  // Generates the level sets given the current cover parameterization, i.e. a list of which points (by index) intersect each level set. 
  List get_level_sets(){
    Function expand_grid = getFunctionR("expand.grid");
    
    // Generate the level set multi-indexes
    List cp_args = List();
    for (int d_i = 0; d_i < d; ++d_i){
      std::vector<int> ls_idx = std::vector<int>(resolution.at(d_i));
      std::iota(ls_idx.begin(), ls_idx.end(), 1);
      std::string key("d"); 
      key += std::to_string(d_i);
      cp_args[key] = ls_idx;
    }
    DataFrame ls_multi_idx = expand_grid(cp_args);
    
    // Convert the indices to vector of indices 
    std::vector< IntegerVector > indices =  std::vector< IntegerVector >();
    for (int d_i = 0; d_i < d; ++d_i){
      std::string key = "d" + std::to_string(d_i);
      indices.push_back( ls_multi_idx[key] );
    }
    
    // Prepare vector of dimension indices
    std::vector<int> dim_idx(d); 
    std::iota(dim_idx.begin(), dim_idx.end(), 0);
    
    // Generate the segment indices 
    std::vector< IntegerVector > ls_segment_idx = std::vector< IntegerVector >(d);
    std::transform(dim_idx.begin(), dim_idx.end(), ls_segment_idx.begin(), [&](const int d_i){
      return(compute_ls_idx(current_cover.at(d_i), d_i));
    });
    
    // Collect the endpoint indices for each level set 
    const int n_level_sets = ls_multi_idx.nrow();
    List res = List(); 
    IntegerVector segment_query_starts = IntegerVector(d);
    IntegerVector segment_query_ends = IntegerVector(d);
    for (int i = 0; i < n_level_sets; ++i){
      std::vector<int> idx = std::vector<int>(d);
      std::string key = "( ";
      std::for_each(dim_idx.begin(), dim_idx.end(), [&, i](const int d_i){
        const int idx = indices.at(d_i).at(i) - 1;
        key += std::to_string(idx + 1) + (d_i == (d - 1) ? " " : ", ");
        const int seg_start = ls_segment_idx.at(d_i).at(idx*2);
        const int seg_end = ls_segment_idx.at(d_i).at((idx*2) + 1);
        segment_query_starts.at(d_i) = seg_start;
        segment_query_ends.at(d_i) = seg_end;
      });
      key += ")";
      res[key] = mtree->query(segment_query_starts, segment_query_ends);
    }
    return(res);
    
    // TODO: go only by endpts indexing instead of segment tree
    // s_size_t n_level_sets = ls_map.size();
    // List res = List(n_level_sets);
    // for (s_size_t l_i = 0; l_i < n_level_sets; ++l_i){
    //   uint_fast8_t i = ls_idx.at(l_i * 2), j = ls_idx.at((l_i * 2) + 1);
    //   res.push_back(mtree->query((int(i), int(j)));
    // }
  }
  
  // List get_base_cover_ls(){
  //   const int n_level_sets = tr->n/2;
  //   std::vector< IntegerVector > res;
  //   res.reserve(n_level_sets);
  //   for (int i = 0; i < tr->n; i += 2){
  //     res.push_back(wrap(tr->queryInterval(i, i+1)));
  //   }
  //   return(wrap(res));
  // }
  
  IntegerVector query(const IntegerVector& l, const IntegerVector& r){
    if (l.size() != d || r.size() != d){ stop("Query must match dimensionality of filter space."); }
    IntegerVector res = mtree->query(l, r);
    return(res);
  }
  
  IntegerVector get_ls_idx(std::size_t i, std::size_t d_i){
    IntegerVector res = no_init(ls_idx.at(d_i).size());
    std::transform(ls_idx.at(d_i).begin(), ls_idx.at(d_i).end(), res.begin(), [](const std::shared_ptr<uint_fast8_t> idx_ptr){
      return(static_cast<int>(*idx_ptr));
    });
    return(res);
  }
  
  // Assigns the values of 'new_idx' to the values pointed to by the shared pointers in the 'current_idx' 
  void assign_ls_idx(const index_t new_idx, const index_ptr_t& current_idx){
    if (new_idx.size() != current_idx.size()){ stop("New indices size do not match current indices."); }
    std::size_t cc = 0; 
    std::for_each(new_idx.begin(), new_idx.end(), [&](const uint_fast8_t idx){
      (*current_idx.at(cc++)) = idx; // dynamically update the endpoint configuration
    });
  }
  
  // Computes an integer vector where each consecutive pair of integers represents the segment intervals the level set represents. 
  // Note that all intervals are stored in the form [s, e), where e > s, and thus the points composing a given level set w/ index 
  // [s, e) are stored in the segments {s, s+1, ..., e - 1}, i.e. not including e. The parameter 'i' controls how many times the 
  // endpoints have 'switched' position, i.e. as the intervals grow uniformly from their starting disjoint configuration, the 
  // endpoints representing the boundaries of each level set occasionally 'swap' with their neighboring sets.  
  index_t compute_ls_idx(std::size_t i, std::size_t d_i){
    if (d_i >= resolution.size()){ stop("Invalid index: Must be 0-based index less than filter dimensionality"); }
    const int n_endpts = resolution.at(d_i)*2;
    const int n_level_sets = n_endpts/2; 
    index_t ls_idx_res = index_t(n_endpts);
    int hc = n_level_sets - 1, lc = 0;
    for (std::size_t j = 0; j < n_level_sets; ++j){
      int j_tmp = n_endpts - j - 1;
      // Rprintf("j: %d, j_tmp: %d, i: %d, lc: %d, hc: %d\n");
      if (j % 2 == 0){
        ls_idx_res[j] = make_index_t(i >= lc ? j - lc : j - i);
        ls_idx_res[j_tmp] = make_index_t(i >= lc ? j_tmp + lc : j_tmp + i);
        ++lc;
      }
      else {
        ls_idx_res[j] = make_index_t(i >= hc ? j + hc : j + i);
        ls_idx_res[j_tmp] = make_index_t(i >= hc ? j_tmp - hc : j_tmp - i);
        --hc;
      }
    }
    return(ls_idx_res);
  }
  
  // Given an integer representing the current swap configuration, generate a 
  // mapping between flat level set indices to which segment endpoints each level set points to
  List get_ls_mapping(){
    
    // Convert the indices to vector of indices 
    std::vector< IntegerVector > indices = get_multi_indexes();
    
    // Prepare vector of dimension indices
    std::vector<int> dim_idx(d); 
    std::iota(dim_idx.begin(), dim_idx.end(), 0);
  
    // Collect the endpoint indices for each level set 
    const int n_level_sets = indices.at(0).size();
    List res = List(); 
    for (int i = 0; i < n_level_sets; ++i){
      std::vector<int> idx = std::vector<int>(d);
      IntegerMatrix segment_idx = IntegerMatrix(d, 2);
      std::string key = "( ";
      std::for_each(dim_idx.begin(), dim_idx.end(), [&, i](const int d_i){
        const int idx = indices.at(d_i).at(i) - 1;
        key += std::to_string(idx + 1) + (d_i == (d - 1) ? " " : ", ");
        const int seg_start = static_cast<int>(*ls_idx.at(d_i).at(idx*2));
        const int seg_end = static_cast<int>(*ls_idx.at(d_i).at((idx*2) + 1));
        segment_idx.row(d_i) = IntegerVector::create(seg_start, seg_end);
      });
      key += ")";
      res[key] = (segment_idx);
    }
    return(res);
  }
  
  // Given an integer representing the current swap configuration, generate a 
  // mapping between flat level set indices to which segment endpoints each level set points to
  List make_ls_mapping(IntegerVector swap_idx){
    if (swap_idx.size() != d){ stop("swap indices size != dimensionality of filter space. "); }
    Function expand_grid = getFunctionR("expand.grid");
    
    // Generate the level set multi-indexes
    List cp_args = List();
    for (int d_i = 0; d_i < d; ++d_i){
      std::vector<int> ls_idx = std::vector<int>(resolution.at(d_i));
      std::iota(ls_idx.begin(), ls_idx.end(), 1);
      std::string key("d"); 
      key += std::to_string(d_i);
      cp_args[key] = ls_idx;
    }
    DataFrame ls_multi_idx = expand_grid(cp_args);
    
    // Convert the indices to vector of indices 
    std::vector< IntegerVector > indices =  std::vector< IntegerVector >();
    for (int d_i = 0; d_i < d; ++d_i){
      std::string key = "d" + std::to_string(d_i);
      indices.push_back( ls_multi_idx[key] );
    }
    
    // Prepare vector of dimension indices
    std::vector<int> dim_idx(d); 
    std::iota(dim_idx.begin(), dim_idx.end(), 0);
    
    // Generate the segment indices 
    std::vector< IntegerVector > ls_segment_idx = std::vector< IntegerVector >(d);
    std::transform(dim_idx.begin(), dim_idx.end(), ls_segment_idx.begin(), [&](const int d_i){
      return(compute_ls_idx(swap_idx.at(d_i), d_i));
    });
    

    // return(List::create(_["ls_segment_idx"] = ls_segment_idx, _["dim_idx"] = dim_idx, _["indices"] = wrap(indices)));
    
    // Collect the endpoint indices for each level set 
    const int n_level_sets = ls_multi_idx.nrow();
    List res = List(); 
    for (int i = 0; i < n_level_sets; ++i){
      std::vector<int> idx = std::vector<int>(d);
      IntegerMatrix segment_idx = IntegerMatrix(d, 2);
      std::string key = "( ";
      std::for_each(dim_idx.begin(), dim_idx.end(), [&, i](const int d_i){
        const int idx = indices.at(d_i).at(i) - 1;
        key += std::to_string(idx + 1) + (d_i == (d - 1) ? " " : ", ");
        const int seg_start = ls_segment_idx.at(d_i).at(idx*2);
        const int seg_end = ls_segment_idx.at(d_i).at((idx*2) + 1);
        segment_idx.row(d_i) = IntegerVector::create(seg_start, seg_end);
      });
      key += ")";
      res[key] = (segment_idx);
    }
    
    // const IntegerVector& ls_idx = compute_ls_idx(swap_idx.at(d_i), d_i);
    //test.at(0);
    return(res);
  }
  
};

RCPP_MODULE(multiscale_module) {
  Rcpp::class_<MultiScale>("MultiScale")
  .constructor<NumericMatrix, Function, SEXP>()
  .method( "set_everything", &MultiScale::set_everything )
  .method( "build_segment_trees", &MultiScale::build_segment_trees )
  .method( "make_ls_map", &MultiScale::make_ls_map )
  .method( "query_ls_map", &MultiScale::query_ls_map )
  .method( "update_cover", &MultiScale::update_cover )
  .method( "query", &MultiScale::query )
  .method( "get_multi_indexes", &MultiScale::get_multi_indexes )
  .method( "get_ls_idx", &MultiScale::get_ls_idx )
  .method( "make_ls_map", &MultiScale::make_ls_map)
  .method( "make_endpt_ls_map", &MultiScale::make_endpt_ls_map)
  .method( "get_endpt_ls_map", &MultiScale::get_endpt_ls_map)
  .method( "make_ls_mapping", &MultiScale::make_ls_mapping)
  .method( "get_ls_mapping", &MultiScale::get_ls_mapping )
  // .method( "get_base_cover_ls", &MultiScale::get_base_cover_ls )
  .method( "compute_ls_idx", &MultiScale::compute_ls_idx )
  .method( "get_level_sets", &MultiScale::get_level_sets )
  .method( "get_segment_idx", &MultiScale::get_segment_idx )
  .method( "build_multiscale_configuration", &MultiScale::build_multiscale_configuration )
  .method( "get_segment_map", &MultiScale::get_segment_map )
  .method( "update_multi_cover", &MultiScale::update_multi_cover )
  ;
}


// Computes Mapper across multiple cover parameterizations. Requires as input:
//  pt_idx := The point indices to change per parameterization
//  from_ls := The level set index the point lies within at the current parameterization
//  to_ls := The level set index the point will intersect at the current parameterization
//  dist_thresh := 1/2th the interval size of each box at the current parameterization
//  clustering_function := the clustering function to apply per parameterization
// [[Rcpp::export]]
void multiscale(const IntegerVector& pt_idx){
  
}



/*** R
*/
