#include <Rcpp.h>
using namespace Rcpp;

#include <array>
#include <cstdint>
#include <stack>
#include <set>
#include <chrono>  // for high_resolution_clock
#include "MultiSegmentTree.h"

#include "cartesian_product.h"

typedef std::vector< uint_fast8_t > index_t;
typedef std::vector< std::shared_ptr<uint_fast8_t> > index_ptr_t;
typedef std::ptrdiff_t s_size_t;

inline std::size_t index_lower_triangular(std::size_t from, std::size_t to, const std::size_t N){
  if (from < to){ std::swap(from, to); }
  return((N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1));
}
#define INDEX_TO(k, n) n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5) // expects 0-based, returns 0-based
#define INDEX_FROM(k, n, i) k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2 // expects 0-based, returns 0-based


std::string index_to_str(index_t idx){
  std::string result;
  std::for_each(idx.begin(), idx.end()-1, [&result](const uint_fast8_t i){
    result.append(std::to_string(static_cast<int>(i)));
    result.append(" ");
  });
  result.append(std::to_string(static_cast<int>(idx.at(idx.size() - 1))));
  return(result);
}

std::string index_as_key(index_t idx, bool one_based = true){
  std::transform(idx.begin(), idx.end(), idx.begin(), std::bind2nd(std::plus<uint_fast8_t>(), static_cast<uint_fast8_t>(one_based)));
  std::string key = "("+index_to_str(idx)+")";
  return(key);
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

// TODO
// template<typename T> 
// std::vector<T> soa_to_aos(){
//   
// }

// Converts a vector of vectors (SoA memory structure) to an IntegerMatrix
template<typename T> 
IntegerMatrix vec_to_int_matrix(std::vector< std::vector<T> > vec){
  const int n = vec.at(0).size(), d = vec.size();
  IntegerMatrix res = Rcpp::no_init_matrix(n, d);
  for (int d_i = 0; d_i < d; ++d_i){
    IntegerVector tmp = wrap(vec.at(d_i));
    res(_, d_i) = tmp;
  }
  return(res);
}


template <typename T> 
std::vector<T> merge_vectors(const std::vector< std::vector<T>* >& vec){
  std::size_t total_vec_size = 0;
  std::for_each(vec.begin(), vec.end(), [&](const std::vector<T>* v){ total_vec_size += v->size(); });
  std::vector< T > final_res = std::vector< T >();
  final_res.reserve(total_vec_size);
  std::for_each(vec.begin(), vec.end(), [&](const std::vector<T>* v){ std::copy(v->begin(), v->end(), std::back_inserter(final_res)); });
  return(final_res);
}

// Applies the function Func to all pairwise combinations in the range [first, last)
template<typename Iter, typename Func>
void combine_pairwise(Iter first, Iter last, Func func)
{
  for(; first != last; ++first){
    for(Iter next = std::next(first); next != last; ++next){
      func(*first, *next);
    }
  }
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
  const int d, n_pts;
  // Rcpp::XPtr<SegmentTree> tr; // the segment tree contains the level sets 
  // std::vector<LS_idx> level_set_idx;
  // std::vector< std::vector<uint> > level_sets;
  std::vector< IntegerVector > pt_idx, from_ls_idx, to_ls_idx, swap_idx; 
  std::vector< NumericVector > dist_to_box, swap_dist;
  // std::vector< uint_fast8_t > ls_idx;
  // std::stack< uint_fast8_t, std::vector< uint_fast8_t > > ls_idx_swaps;
  std::vector< s_size_t > current_index; 
  index_t current_cover; 
  IntegerVector resolution; 
  std::shared_ptr<MultiSegmentTree> mtree; 
  std::vector< index_ptr_t > ls_idx; // tracks the per dimension endpt configuration
  std::map< index_t, std::size_t > ls_to_lsfi_map;
  std::map< index_ptr_t, uint_fast8_t > endpt_to_ls_map;
  std::vector< std::vector< uint_fast8_t > > pt_segment_idx; 
  
  std::map< index_t, std::vector< int > > segment_map; 
  
  std::vector< index_t> lsmi_AoS;
  StringVector index_set; 
  // std::vector< IntegerVector > segment_pts;
  std::vector<int> dim_idx;
  
  MultiScale(const NumericMatrix& X, const Function f, const IntegerVector& k) : data(X), n_pts(data.nrow()), clustering_function(f), d(k.size()) { //, tr(stree) {
    bool res_check = std::any_of(resolution.begin(), resolution.end(), [](const int k_i){ return(k_i > 255); });
    if (res_check){ stop("Multiscale version of mapper only supports up to 255 intervals per dimension."); }
    //level_set_idx = std::vector<LS_idx>();
    current_index = std::vector< s_size_t >(d, -1); // points to the index of the filtration that has not been computed. 
    current_cover = index_t(d, 0);
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
    
    // Useful vector to have
    dim_idx = std::vector<int>(d);
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
    
    std::vector< std::vector<uint_fast8_t> > lsmi_SoA(d);
    std::for_each(dim_idx.begin(), dim_idx.end(), [&lsmi_SoA, &k](const int d_i){
      std::vector<uint_fast8_t> idx = std::vector<uint_fast8_t>(k.at(d_i));
      std::iota(idx.begin(), idx.end(), 0);
      lsmi_SoA.at(d_i) = idx;
    });
    // std::for_each(dim_idx.begin(), dim_idx.end(), [&lsmi_SoA, &k](const int d_i){
    //   std::for_each(idx.begin(), idx.end(), [](const uint_fast8_t el){
    //     Rcout << static_cast<int>(el); 
    //   });
    //   Rcout << std::endl; 
    // });

    //const MultiScale& tmp = *this;
    
    
    
    lsmi_AoS = std::vector<index_t>(); 
    CartesianProduct<uint_fast8_t>(lsmi_SoA, [&](const index_t lsmi){
      lsmi_AoS.push_back(lsmi);
    });
    
    
    // TODO: inspect these
    make_ls_map();
    make_endpt_ls_map();
  }
  
  void set_index_set(const StringVector& _index_set){
    index_set = _index_set;
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
    ls_to_lsfi_map = std::map< index_t, std::size_t >();
    
    // Prepare vector of dimension indices to use for each transform
    std::vector<uint_fast8_t> dim_idx(d); 
    std::iota(dim_idx.begin(), dim_idx.end(), 0);
  
    // Fill in the mapping from index set --> counter 
    const std::size_t n_level_sets = ls_multi_idx.at(0).size();
    for (std::size_t i = 0; i < n_level_sets; ++i){
      index_t idx_key = index_t(d);
      std::for_each(dim_idx.begin(), dim_idx.end(), [&, i](const int_fast8_t d_i){
        const uint_fast8_t idx = (uint_fast8_t) ls_multi_idx.at(d_i).at(i) - 1;
        idx_key.at(d_i) = idx;
      });
      ls_to_lsfi_map.emplace(idx_key, i);
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
    index_t query(d); 
    std::copy(query_idx.begin(), query_idx.end(), query.begin());
    return(static_cast<int>(ls_to_lsfi_map.at(query)));
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
  
  // Builds the initial set of segment to point index mappings. Generates 1-based point indices, 0-based segment indices. 
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
        it->second.push_back(static_cast<int>(i+1));
      } else {
        std::vector< int > v = { static_cast<int>(i+1) };
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
  
  // Retrieves the last known target segment a point would've merged to, w/ increasing overlap,
  // up to and including the given target index of the filtration.
  std::map<int, int> last_known_target_segment(std::vector<int> pt_ids, const int target_idx, const int d_i){
    const IntegerVector& c_swap_idx = swap_idx.at(d_i);
    const IntegerVector& c_pt_idx = pt_idx.at(d_i);
    const IntegerVector& c_from_ls = from_ls_idx.at(d_i);
    const IntegerVector& c_to_ls = to_ls_idx.at(d_i);
    std::vector< int > pt_idx_vec = std::vector<int>(c_pt_idx.begin(), c_pt_idx.begin() + target_idx + 1);
    
    // Retrieve the unique point indices that last changed along the range, going backwards
    std::map<int, int> pt_uniq_idx = get_unique_indices(pt_idx_vec.rbegin(), pt_idx_vec.rend()); 
    
    // For each point, retrieve its last target segment
    std::map<int, int> pt_to_res; 
    for (auto& pidx: pt_uniq_idx){
      const int global_idx = pt_idx_vec.size() - pidx.second - 1; // relative index 
      index_t new_ls_idx = compute_ls_idx(c_swap_idx.at(global_idx), d_i);
      const int from_ls = c_from_ls.at(global_idx), to_ls = c_to_ls.at(global_idx);
      const int pt = pidx.first;
      const int to_segment = (from_ls < to_ls ? new_ls_idx.at(to_ls*2) : new_ls_idx.at(to_ls*2 + 1) - 1); 
      pt_to_res.emplace(pidx.first, to_segment);
    }
    
    // Any point not in the map is in its default position
    std::for_each(pt_ids.begin(), pt_ids.end(), [&](const int pt){
      std::map<int, int>::iterator it = pt_to_res.find(pt);
      if (it == pt_to_res.end()){ // pt wasn't found in the map
        std::size_t first_pt_idx = std::distance(pt_idx.at(d_i).begin(), std::find(pt_idx.at(d_i).begin(), pt_idx.at(d_i).end(), pt));
        int first_ls = from_ls_idx.at(d_i).at(first_pt_idx);
        pt_to_res.emplace_hint(it, pt, first_ls*2); // every point initially exists in the lower segment of its initial level set
      }
    });
    
    // Return the map
    return(pt_to_res);
  }
  
  List test(){
    // std::vector<int> test_pts = std::vector<int>(10);
    // std::iota(test_pts.begin(), test_pts.end(), 1);
    // std::map<int, int> last_from = last_known_ls(test_pts, 3, 10, 0);
    // List res = List();
    // for (auto& pair: last_from){
    //   res[std::to_string(pair.first).c_str()] = pair.second;
    // }
    // return(res);
  }
  
  // A point may exist in multiple level sets. It is necessary to know the full level set index a point 
  // last intersected in the filtration (from the expansion sense).
  std::map<int, int> last_known_ls(std::vector<int> pt_ids, const int start_idx, const int d_i){
    const IntegerVector& c_swap_idx = swap_idx.at(d_i);
    const IntegerVector& c_pt_idx = pt_idx.at(d_i);
    const IntegerVector& c_from_ls = from_ls_idx.at(d_i);
    const IntegerVector& c_to_ls = to_ls_idx.at(d_i);
    
    // Main result
    std::map<int, int> pt_to_res; // For each point, retrieve its last to ls index
    
    // If at the initial configuration where the level sets are disjoint 
    if (start_idx == -1){
      std::for_each(pt_ids.begin(), pt_ids.end(), [&](const int pt){
        // Extract first occurrence of point in point ids. It's from level set is where it came from originally.
        std::size_t first_pt_idx = std::distance(c_pt_idx.begin(), std::find(c_pt_idx.begin(), c_pt_idx.end(), pt));
        int first_ls = c_from_ls.at(first_pt_idx);
        pt_to_res.emplace(pt, first_ls); // every point initially exists in the lower segment of its initial level set
      });
      return(pt_to_res);
    }
    
    // Subrange of the filtration to focus on which is less or equal to the current starting index
    std::vector< int > pt_idx_vec = std::vector<int>(c_pt_idx.begin(), c_pt_idx.begin() + start_idx);
    
    // Retrieve the unique point indices that last changed along the range, going backwards
    std::map<int, int> pt_uniq_idx = get_unique_indices(pt_idx_vec.rbegin(), pt_idx_vec.rend()); 
    for (auto& pidx: pt_uniq_idx){
      const int global_idx = pt_idx_vec.size() - pidx.second - 1; // relative index 
      index_t new_ls_idx = compute_ls_idx(c_swap_idx.at(global_idx), d_i);
      const int to_ls = c_to_ls.at(global_idx);
      pt_to_res.emplace(pidx.first, to_ls);
    }
    
    // Any point not in the map is in its default position
    std::for_each(pt_ids.begin(), pt_ids.end(), [&](const int pt){
      std::map<int, int>::iterator it = pt_to_res.find(pt);
      if (it == pt_to_res.end()){ // pt wasn't found in the map
        // Extract first occurrence of point in point ids. It's from level set is where it came from originally.
        std::size_t first_pt_idx = std::distance(pt_idx.at(d_i).begin(), std::find(pt_idx.at(d_i).begin(), pt_idx.at(d_i).end(), pt));
        int first_ls = from_ls_idx.at(d_i).at(first_pt_idx);
        pt_to_res.emplace_hint(it, pt, first_ls); // every point initially exists in the lower segment of its initial level set
      }
    });
    
    // Return the map
    return(pt_to_res);
  }

  std::vector< std::vector<std::size_t> > unexpand_lsfi_pairs(const std::set<std::size_t>& lsfi_pairs) {
    // Unexpand the the pairs to their lsfis...
    // IntegerMatrix res = no_init_matrix(ls_pair_idx.size(), 2);
    std::vector< std::vector<std::size_t> > res = std::vector< std::vector<std::size_t> >(2);
    const std::size_t n_level_sets = ls_to_lsfi_map.size(); 
    std::size_t i = 0; 
    std::for_each(lsfi_pairs.begin(), lsfi_pairs.end(), [&i, &res, &n_level_sets](const std::size_t idx){
      std::size_t to = INDEX_TO(idx, n_level_sets);
      std::size_t from = INDEX_FROM(idx, n_level_sets, to);
      res.at(0).push_back(to);
      res.at(1).push_back(from);
      // res(i++, _) = IntegerVector::create(to, from);
    });
    return(res);
  }
  
  
  // Returns the level set pairs that have changed state given the current filtration state 'start_idx' 
  // and the end filtration state 'target_idx'. The output is a set of index-types which uniquely 
  // represent a pair of level set indices. The set must be unexpanded into a two column matrix of LSFIs 
  // to, for example, actually retrieve the level sets. 
  std::pair< std::set<std::size_t>, std::set<std::size_t> > get_ls_that_change(const IntegerVector start_idx, const IntegerVector target_idx){
    if (target_idx.size() != d){ stop("Invalid query, doesn't match dimensionality of filter space."); }
    int d_i = 0; 
    bool invalid = std::any_of(target_idx.begin(), target_idx.end(), [&](const int ti){ 
      return(ti < -1 || ti >= pt_idx.at(d_i++).size());
    });

    // Get all the point ids 
    std::vector<int> all_pts = std::vector<int>(n_pts);
    std::iota(all_pts.begin(), all_pts.end(), 1);
    
    // To get the level set pairs to update, need to know the range of the level set expansion for each point
    std::vector< std::vector< std::pair<int, int> > > min_max_ls(d);
    
    // Retrieve the last known level sets each point intersected before the current starting index 
    // std::vector< std::vector<int> > pt_ls_pos(d);
    for (int& d_i : dim_idx){
      std::map<int, int> ls_pos = last_known_ls(all_pts, start_idx.at(d_i), d_i);
      min_max_ls.at(d_i) = std::vector< std::pair<int, int> >(n_pts);
      for (const auto& kv : ls_pos){
        Rprintf("Default range (i=%d, d=%d): [%d, %d]\n", kv.first, d_i, kv.second, kv.second);
        min_max_ls.at(d_i).at(kv.first-1).first = kv.second;
        min_max_ls.at(d_i).at(kv.first-1).second = kv.second;
      } 
      // std::transform(ls_pos.begin(), ls_pos.end(), min_max_ls.begin(), [&](std::pair<int, int>& kv){
      //   // std::pair<int, int> ls_range = std::make_pair(kv.second, kv.second);
      // 
      // });
      // // pt_ls_pos.at(d_i) = std::vector<int>(n_pts);
      // for (auto& pt_to_ls: ls_pos){
      //   // pt_ls_pos.at(d_i).at(pt_to_ls.first-1) = pt_to_ls.second;
      //   
      // }
    }
    
    // Compute the min/max LS indexes in the current range for each point
    // std::vector< std::vector< std::pair<int, int> > > min_max_ls(d);
    // Rcout << "Computing the ls indices..." << std::endl;  
    for (int& d_i : dim_idx){
      
      // Vector mapping a point id to its min/max ls index within the current range
      // std::pair<int, int> default_range = std::make_pair<int, int>(std::numeric_limits<int>::max(), std::numeric_limits<int>::min()); 
      // min_max_ls.at(d_i) = std::vector< std::pair<int, int> >(n_pts, default_range));
      
      // Loop through the current filtration range, updating the min/max ls bounds for each point.
      for (int i = start_idx.at(d_i)+1; i <= target_idx.at(d_i); ++i){
        // Rprintf("pt index: %d\n", pt_idx.at(d_i).at(i)-1);
        std::pair<int, int>& c_range = min_max_ls.at(d_i).at(pt_idx.at(d_i).at(i)-1);
        const int from_ls = from_ls_idx.at(d_i).at(i);
        const int to_ls = to_ls_idx.at(d_i).at(i);
        if (from_ls < to_ls){
          if (from_ls < c_range.first){ c_range.first = from_ls; }
          if (to_ls > c_range.second){ c_range.second = to_ls; }
        } else {
          if (to_ls < c_range.first){ c_range.first = to_ls; }
          if (from_ls > c_range.second){ c_range.second = from_ls; }
        }
      }
    }
    
    // Debugging 
    // IntegerMatrix res = IntegerMatrix(n_pts, d*2);
    // IntegerVector tmp = IntegerVector(d*2);
    // for (int i = 0; i < n_pts; ++i){
    //   std::for_each(dim_idx.begin(), dim_idx.end(), [&](const int d_i){
    //     std::pair<int, int>& c_range = min_max_ls.at(d_i).at(i);
    //     tmp.at(d_i*2) = c_range.first;
    //     tmp.at(d_i*2 + 1) = c_range.second;
    //   });
    //   res(i, _) = clone(tmp);
    // }
    
    // Rcout << "Collecting the expansions..." << std::endl;  
    // return(res);
    std::set<std::size_t> lsfi_to_recompute = std::set<std::size_t>();
    std::set<std::size_t> ls_pair_idx = std::set<std::size_t>();
    std::vector< std::vector<uint_fast8_t> > expansions = std::vector< std::vector<uint_fast8_t> >(d);
    std::vector<uint_fast8_t> expanded = std::vector<uint_fast8_t>();
    for (int i = 0; i < n_pts; ++i){
      expanded.clear();
      std::for_each(dim_idx.begin(), dim_idx.end(), [&](const int d_i){
        std::pair<int, int>& c_range = min_max_ls.at(d_i).at(i);
        // tmp.at(d_i*2) = c_range.first;
        // tmp.at(d_i*2 + 1) = c_range.second;
        // Expand the index inclusion range
        expanded.resize((c_range.second - c_range.first) + 1);
        std::iota(expanded.begin(), expanded.end(), c_range.first);
        expansions.at(d_i) = expanded;
      });
      // Rcout << "Expansion done..." << std::endl;  
      
      // Convert the cartesian product of the expanded pairs to level set flat indices
      std::vector<std::size_t> changed_lsfis = std::vector<std::size_t>();
      std::size_t cc = 0;
      std::for_each(expansions.begin(), expansions.end(), [&cc](const std::vector<uint_fast8_t>& v){ 
        cc += v.size();
      });
      if (cc > d){ 
        CartesianProduct(expansions, [&](const index_t& ls_idx){
          changed_lsfis.push_back(ls_to_lsfi_map.at(ls_idx));
          Rprintf("LSFI %d changed due to pt %d\n", ls_to_lsfi_map.at(ls_idx), i+1);
        });
      }
      
      // Record the unique set of LSFIs that are changing
      std::for_each(changed_lsfis.begin(), changed_lsfis.end(), [&lsfi_to_recompute](const std::size_t lsfi){
        lsfi_to_recompute.insert(lsfi);
      });
      
      // Rcout << "Changed LSFIs collected..." << std::endl;  
      
      // The cartesian product of these LSFIs should be updated
      if (changed_lsfis.size() > 0){
        std::vector< std::vector<std::size_t> > lsfi_lst(2);
        lsfi_lst.at(0) = changed_lsfis;
        lsfi_lst.at(1) = changed_lsfis;
        const std::size_t n_level_sets = ls_to_lsfi_map.size(); 
        CartesianProduct(lsfi_lst, [&ls_pair_idx, &n_level_sets](const std::vector<std::size_t> ls_pair){
          std::size_t i = ls_pair.at(0);
          std::size_t j = ls_pair.at(1);
          if (i != j){
            ls_pair_idx.insert(index_lower_triangular(i, j, n_level_sets));
            Rprintf("Adding LS pair (%d, %d)\n", i, j);
          }
        });
      }
      
      
      // res(i, _) = clone(tmp);
    }
    
    
    // Rcout << "Returning..." << std::endl;  
    // IntegerVector res = IntegerVector(ls_pair_idx.begin(), ls_pair_idx.end());
    std::pair< std::set<std::size_t>, std::set<std::size_t> > res(lsfi_to_recompute, ls_pair_idx);
    return(res);
    // return(vectors_to_matrix(pt_ls_pos));
    
    // IntegerMatrix A_tmp = clone(A);
    // std::vector< std::vector<int> > all_from_ls, all_to_ls;
    // 
    // std::map<int, std::pair<int, int> > (); 
    // for (d_i = 0; d_i < d; ++d_i){
    //   
    //   // Which direction is the filtration going, expanding or contracting? 
    //   const bool is_increasing = current_index.at(d_i) < target_index.at(d_i);
    //   
    //   // Extract the range [start, end] of the filtration containing the information needed to change the mapper
    //   const int start_idx = is_increasing ? current_index.at(d_i)+1 : target_index.at(d_i)+1; 
    //   const int end_idx = is_increasing ? target_index.at(d_i) : current_index.at(d_i); // inclusive
    //   
    //   // Get the current information to extract 
    //   const IntegerVector& c_swap_idx = swap_idx.at(d_i), 
    //       &c_pt_idx = pt_idx.at(d_i), 
    //       &c_from_ls = from_ls_idx.at(d_i), 
    //       &c_to_ls = to_ls_idx.at(d_i);
    //   
    //   // Record the level sets to update, and the level set pairs 
    //   std::vector< index_t > ls_to_update;
    //   std::vector< std::pair<index_t, index_t> > ls_pairs_to_update;
    //   if (is_increasing){
    //     
    //     for (int i = start_idx; i < end_idx; ++i){
    //       int from_ls = c_from_ls.at(i), to_ls = c_to_ls.at(i);
    //       
    //       // Extract the index of the actual level set the pt is at/going to
    //       // But! A point exists in multiple level sets at one time
    //       // Here, extract the 
    //       // ls_to_update.push_back();
    //     }
    //   }
    // 
    //     
    //   std::vector< int > pt_idx_vec = std::vector<int>(c_pt_idx.begin() + start_idx, c_pt_idx.begin() + end_idx + 1);
    //   std::map<int, int> pt_uniq_idx_forwards, pt_uniq_idx_backwards; 
    //   std::vector< std::map<int, int> > start_ls(d), end_ls(d);
    //   if (is_increasing){
    //     
    //     pt_uniq_idx_forwards = get_unique_indices(pt_idx_vec.begin(), pt_idx_vec.end()); 
    //     pt_uniq_idx_backwards = get_unique_indices(pt_idx_vec.rbegin(), pt_idx_vec.rend()); 
    //     
    //     for (auto& pidx: pt_uniq_idx_backwards){
    //       const int pt = pidx.first;
    //       const int global_idx = start_idx + (pt_idx_vec.size() - pidx.second - 1); // global pt index
    //       const int from_ls = c_from_ls.at(global_idx), to_ls = c_to_ls.at(global_idx);
    //       // A_tmp.at(pt-1, d_i) = to_ls+1;
    //     }
    //     for (auto& pidx: pt_uniq_idx_forwards){
    //       const int pt = pidx.first;
    //       const int global_idx = start_idx + pidx.second;  // global pt index
    //       const int from_ls = c_from_ls.at(global_idx), to_ls = c_to_ls.at(global_idx);
    //     }
    //   }
    // }
    // return(A_tmp);
  }
  
  // Iterate through all pairwise combinations of level sets, accepting only combination pairs that are 
  // sufficiently close, as determined by ls_config.
  IntegerMatrix ls_to_change(const IntegerVector ls_config){
    using kv = std::map<index_t, std::size_t>::value_type;
    std::vector<uint_fast8_t>::iterator max_dif;
    std::vector<bool> max_diff(d);
    IntegerVector c1, c2; 
    combine_pairwise(ls_to_lsfi_map.begin(), ls_to_lsfi_map.end(), [&](kv& from, kv& to){
      std::size_t d_i = 0; 
      std::transform(from.first.begin(), from.first.end(), to.first.begin(), max_diff.begin(), [&d_i, &ls_config](const uint_fast8_t i, const uint_fast8_t j){
        return(std::abs(i - j) <= ls_config.at(d_i++));
      });
      bool valid_pair = std::all_of(max_diff.begin(), max_diff.end(), [](const bool res){ return(res); });
      if (valid_pair){
        c1.push_back(static_cast<int>(from.second));
        c2.push_back(static_cast<int>(to.second));
      }
    });
    IntegerMatrix res(c1.size(), 2);
    res(_, 0) = c1, res(_, 1) = c2;
    return(res);
  }
  
  void ls_change(const int d_i, index_t& key, index_t& ls_config){
    // if (d_i > d){ return; }
    // else if (d_i == d){
    //   for (int i = 0; i < resolution.at(d_i); ++i){
    //     key.at(d_i) = i;
    //     std::vector<int> tmp(ls_config.at(d_i)*2 + 1); 
    //     std::iota(tmp.begin(), tmp.end(), i - ls_config.at(d_i));
    //     
    //   }
    // }
    // else {
    //   for (int i = 0; i < resolution.at(d_i); ++i){ 
    //     key.at(d_i) = i;
    //     ls_change(d_i, key, ls_config); 
    //   }
    // }
    
  }
  
  // Assumes the level set configuration has been set using the current cover 
  std::vector< IntegerVector > run_clustering(std::set<std::size_t> lsfi_to_compute, const std::vector< index_t >& ls_segment_idx){
    
    // Convert the LSFIs to LSMIs
    std::vector< index_t > lsmi_vec = std::vector< index_t >(lsfi_to_compute.size());
    std::transform(lsfi_to_compute.begin(), lsfi_to_compute.end(), lsmi_vec.begin(), [&](const std::size_t lsfi){
      return(lsmi_AoS.at(lsfi));
    });
    
    
    // Retrieve the level set membership vectors
    // TODO: investigate peformance of doing one level set at a time vs. all 
    std::vector< IntegerVector > all_vertices = std::vector< IntegerVector >(); 
    return(all_vertices);
    std::size_t cc = 0; 
    std::for_each(lsmi_vec.begin(), lsmi_vec.end(), [&](const index_t lsmi){
      
      // Get the point indices we're working with
      std::vector<int> ls_pt_idx = get_level_set(lsmi, ls_segment_idx);
      
      // Apply the clustering 
      const IntegerVector cl_results = clustering_function(_["X"] = data, _["idx"] = wrap(ls_pt_idx));
      const IntegerVector cl_idx = self_match(cl_results) - 1;  
      const IntegerVector ids = unique(cl_idx);
      all_vertices.resize(all_vertices.size() + ids.size());
      std::vector<int>::iterator c_pt = ls_pt_idx.begin(); 
      
      // Fill in the point memberships
      for (IntegerVector::const_iterator it = cl_idx.begin(); it != cl_idx.end(); ++it, ++c_pt){
        all_vertices.at(cc + (*it)).push_back((*c_pt));
      }
      cc += ids.size();
    });
    
    // Return
    return(all_vertices);
  }
  
  // The primary API function: Updates the filtration indices to retrieve the information needed to 
  // to compute the mapper at the target filtration index. 
  List compute_mapper(const IntegerVector target_index){
    
    // Get the current index fo the filtration
    IntegerVector current_index_copy = IntegerVector(current_index.begin(), current_index.end());
  
    // Update the segment containers
    // Rcout << "here1" << std::endl; 
    update_multi_cover(target_index);
    
    Rcout << "Current index: " << current_index_copy << std::endl; 
    Rcout << "Target index: " << target_index << std::endl; 
    
    // Extract the LSFI indices and the level set pairs that changed during this update
    using index_set = std::set< std::size_t >;
    // Rcout << "here2" << std::endl; 
    std::pair< index_set, index_set > to_update = get_ls_that_change(current_index_copy, target_index);
    
    Rcout << "Vertices to update: " << std::endl; 
    std::for_each(to_update.first.begin(), to_update.first.end(), [](const std::size_t v){
      Rcout << static_cast<int>(v) << ", ";
    });
    Rcout << std::endl; 
    
    std::vector< std::vector<std::size_t> > lsfi_pairs = unexpand_lsfi_pairs(to_update.second);
    return(wrap(lsfi_pairs));
    // Rcout << "Level set pairs to update: " << std::endl; 
    // std::for_each(lsfi_pairs.begin(), lsfi_pairs.end(), [](const std::vector<std::size_t> pair){
    //   Rcout << static_cast<int>(pair.at(0)) << ", " << static_cast<int>(pair.at(1));
    // });
    // Rcout << std::endl; 
    
    
    
    // Apply the clustering function to update the 0-skeleton 
    Rcout << "here3" << std::endl; 
    std::vector<index_t> c_ls_config = get_current_ls_idx();
    std::vector<IntegerVector> updated_vertices = run_clustering(to_update.first, c_ls_config);
    
    // Retrieve the level set pairs to update
    Rcout << "here4" << std::endl; 
    // std::vector< std::vector<std::size_t> > lsfi_pairs = unexpand_lsfi_pairs(to_update.second);
    
    // Return 
    Rcout << "here5" << std::endl; 
    return(List::create(_["new_vertices"] = updated_vertices, _["lsfi_pairs"] = lsfi_pairs));
    
    // Use the level set pairs that were known to have changed to update the 1-skeleton  
    
  }
  
  // Main functionality of the class. Given a set of indices into the filtration, updates the segment mapping to reflect that 
  // progression state of the cover. 
  void update_multi_cover(const IntegerVector target_index){
    if (target_index.size() != d){ stop("Invalid query, doesn't match dimensionality of filter space."); }
    int d_i = 0; 
    bool invalid = std::any_of(target_index.begin(), target_index.end(), [&](const int ti){ 
      return(ti < -1 || ti >= pt_idx.at(d_i++).size());
    });
    d_i = 0;
    bool all_equal = std::all_of(target_index.begin(), target_index.end(), [&](const int ti){ return(ti == current_index.at(d_i++)); });
    if (all_equal) { return; }
    std::vector< std::vector<int> > pts_to_change_per_d(d); 
    std::vector< std::vector<uint_fast8_t> > pt_seg_from(d); 
    std::vector< std::vector<uint_fast8_t> > pt_seg_to(d); 
    std::unordered_set<int> pts_to_change;
    for (d_i = 0; d_i < d; ++d_i){
      const bool is_increasing = current_index.at(d_i) < target_index.at(d_i); // Otherwise, check which direction we're going
      const int start_idx = is_increasing ? current_index.at(d_i)+1 : target_index.at(d_i)+1; // inclusive
      const int end_idx = is_increasing ? target_index.at(d_i) : current_index.at(d_i); // inclusive
      // Rprintf("is increasing? %s  ==> start: %d, end: %d\n", is_increasing ? "YES" : "NO", start_idx, end_idx);
      const IntegerVector& c_swap_idx = swap_idx.at(d_i);
      const IntegerVector& c_pt_idx = pt_idx.at(d_i);
      const IntegerVector& c_from_ls = from_ls_idx.at(d_i);
      const IntegerVector& c_to_ls = to_ls_idx.at(d_i);
      std::vector< int > pt_idx_vec = std::vector<int>(c_pt_idx.begin() + start_idx, c_pt_idx.begin() + end_idx + 1);
      
      // IntegerVector pt_ids = wrap(pt_idx_vec);
      // Rcout << pt_ids << std::endl; 
      
      std::map<int, int> pt_uniq_idx; 
      if (is_increasing){
        pt_uniq_idx = get_unique_indices(pt_idx_vec.rbegin(), pt_idx_vec.rend()); 
        for (auto& pidx: pt_uniq_idx){
          const int pt = pidx.first;
          const int global_idx = start_idx + (pt_idx_vec.size() - pidx.second - 1); // global pt index
          const int from_ls = c_from_ls.at(global_idx), to_ls = c_to_ls.at(global_idx);
          
          // Update the level set configuration, retrieve which segment(s) the point is going to/coming from
          index_t new_ls_idx = compute_ls_idx(c_swap_idx.at(global_idx), d_i);
          const int from_segment = pt_segment_idx.at(d_i).at(pt - 1);
          const int to_segment = (from_ls < to_ls ? new_ls_idx.at(to_ls*2) : new_ls_idx.at(to_ls*2 + 1) - 1);
          
          // Record the changes that need to happen
          Rprintf("pt id: %d (global=%d) going from ls %d to %d (segment %d to %d)\n", pidx.first, global_idx, from_ls, to_ls, from_segment, to_segment);
          pts_to_change_per_d.at(d_i).push_back(pt);
          pt_seg_from.at(d_i).push_back(from_segment);
          pt_seg_to.at(d_i).push_back(to_segment);
          pts_to_change.insert(pt);
          
          // Update the 'from' segment cache 
          pt_segment_idx.at(d_i).at(pt - 1) = to_segment;
        }
      } else {
        pt_uniq_idx = get_unique_indices(pt_idx_vec.begin(), pt_idx_vec.end());
        std::vector<int> target_pt_ids = std::vector<int>();
        target_pt_ids.reserve(pt_uniq_idx.size());
        for(std::map<int,int>::iterator it = pt_uniq_idx.begin(); it != pt_uniq_idx.end(); ++it) {
          target_pt_ids.push_back(it->first);
        }
        std::map<int, int> pt_target_seg = last_known_target_segment(target_pt_ids, target_index.at(d_i), d_i);
        
        for (auto& pidx: pt_uniq_idx){
          const int pt = pidx.first;
          const int global_idx = start_idx + pidx.second; // global pt index
          const int from_ls = c_from_ls.at(global_idx), to_ls = c_to_ls.at(global_idx);
          
          // Update the level set configuration, retrieve which segment(s) the point is going to/coming from
          index_t new_ls_idx = compute_ls_idx(c_swap_idx.at(global_idx), d_i);
          const int from_segment = pt_segment_idx.at(d_i).at(pt - 1);
          const int to_segment = pt_target_seg.at(pt);
          
          // Record the changes that need to happen
          Rprintf("pt id: %d (global=%d) going from ls %d to %d (segment %d to %d)\n", pidx.first, global_idx, from_ls, to_ls, from_segment, to_segment);
          pts_to_change_per_d.at(d_i).push_back(pt);
          pt_seg_from.at(d_i).push_back(from_segment);
          pt_seg_to.at(d_i).push_back(to_segment);
          pts_to_change.insert(pt);
          
          // Update the 'from' segment cache 
          pt_segment_idx.at(d_i).at(pt - 1) = to_segment;
        }
      }
    }
    // Rcout << "finsihg preprocessing " << std::endl; 
    // Prepare vector of dimension indices
    std::vector<int> dim_idx(d);
    std::iota(dim_idx.begin(), dim_idx.end(), 0);
    
    // Swap the points between the appropriate segments!  
    index_t key_from(d), key_to(d); 
    std::for_each(pts_to_change.begin(), pts_to_change.end(), [&](const int pt){
      std::for_each(dim_idx.begin(), dim_idx.end(), [&](const int d_i){
        std::vector<int> candidate_pts = pts_to_change_per_d.at(d_i);
        ptrdiff_t pos = std::distance(candidate_pts.begin(), std::find(candidate_pts.begin(), candidate_pts.end(), pt));
        if (pos == candidate_pts.size()) { // point doesn't exist in this dimension; use the segment it's at
          // Rprintf("Using point %d's (dim=%d) default current positon %d\n", pt, d_i, pt-1);
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
      std::vector<int>& from_pts = segment_map[key_from];
      std::vector<int>::iterator from_end = std::remove_if(from_pts.begin(), from_pts.end(), [pt](const int x_i){ return(x_i == pt); });
      from_pts.resize(std::distance(from_pts.begin(), from_end));
      segment_map[key_to].push_back(pt);
    });
    
    // Return the level sets that changed in the update. 
    // LogicalVector ls_to_include = LogicalVector(index_set.size(), FALSE);
    // std::for_each(dim_idx.begin(), dim_idx.end(), [&](const int d_i){
    //   
    // });
    // for (int i = current_index.at(d_i))

    
    // If we updated all the way to start, set cover to original disjoint cover
    for (d_i = 0; d_i < d; ++d_i){
      current_index.at(d_i) = target_index.at(d_i);
      const int ls_config = current_index.at(d_i) == -1 ? 0 : swap_idx.at(d_i).at(current_index.at(d_i));
      current_cover.at(d_i) = ls_config;
      index_t new_ls_idx = compute_ls_idx(ls_config, d_i);
      assign_ls_idx(new_ls_idx, ls_idx.at(d_i)); // dynamically update the indices
    }
    
    // return(List::create(_["max_ls_diff"] = wrap(current_cover)));
    
    
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
  

  
  // List get_level_sets(){
  //   
  // }
  
  std::vector<int> get_level_set(const index_t lsmi, const std::vector< index_t >& ls_segment_idx){
    
    // Generate the segment index ranges in the form [start, end) per dimension
    std::vector< std::vector<uint_fast8_t > > segment_ranges(d); 
    std::for_each(dim_idx.begin(), dim_idx.end(), [&](const int d_i){
      const int idx = static_cast<int>(lsmi.at(d_i));
      const int seg_start = ls_segment_idx.at(d_i).at(idx*2);
      const int seg_end = ls_segment_idx.at(d_i).at((idx*2) + 1);
      segment_ranges.at(d_i).resize(std::abs((seg_end - seg_start)));
      std::iota(segment_ranges.at(d_i).begin(), segment_ranges.at(d_i).end(), seg_start);
    });
    
    // Generate the segment endpoint ranges for the given level set
    std::vector< std::vector<int>* > segments_in_ls;
    CartesianProduct(segment_ranges, [&](const index_t current_segment){
      if (segment_map.find(current_segment) != segment_map.end()){
        segments_in_ls.push_back(&segment_map.at(current_segment));
      }
    });
    
    // Merge the results
    return(merge_vectors(segments_in_ls));
  }
  
  // Generates the level sets given the current cover parameterization, i.e. a list of which points (by index) intersect each level set. 
  std::map< index_t, std::vector<int> > gen_level_sets(const std::vector<index_t> query_lsmi, index_t ls_config){
    
    // Generate the segment indices 
    std::vector< index_t > ls_segment_idx = std::vector< index_t >(d);
    std::transform(dim_idx.begin(), dim_idx.end(), ls_segment_idx.begin(), [&](const int d_i){
      return(compute_ls_idx(ls_config.at(d_i), d_i));
    });
    
    // Collect the endpoint indices for each level sets
    const int n_level_sets = query_lsmi.size();
    std::map< index_t, std::vector<int> > res = std::map<index_t, std::vector<int>>(); 
    for (index_t lsmi: query_lsmi){
      // index_t lsmi = kv.first; // level set multi index
      // Rprintf("Current LS: %s\n", index_to_str(lsmi).c_str());
      std::vector<int> ls_pt_membership = get_level_set(lsmi, ls_segment_idx);
      
      // Save the results and move to the next level set 
      // std::string key = index_as_key(lsmi);
      res.emplace(lsmi, ls_pt_membership);
      
      // const std::size_t n_segments = all_segments.at(0).size();
      // index_t current_segment = index_t(d);
      // std::vector< std::vector<int>* > segments;
      // segments.reserve(n_segments);
      // 
      // for (int s_i = 0; s_i < n_segments; ++s_i){
      //   // Extract the current segment
      //   // Rcout << "current segment:";
      //   std::transform(dim_idx.begin(), dim_idx.end(), current_segment.begin(), [&all_segments, s_i](const int d_i){
      //     // Rcout << all_segments.at(d_i).at(s_i) << " ";
      //     return(all_segments.at(d_i).at(s_i));
      //   });
      //   // Rcout << std::endl; 
      // 
      //   if (segment_map.find(current_segment) != segment_map.end()){
      //     segments.push_back(&segment_map.at(current_segment));
      //   }
      // }
      // 
      // // Save the results and move to the next level set 
      // index_t lsmi_key = lsmi;
      // std::transform(lsmi_key.begin(), lsmi_key.end(), lsmi_key.begin(), std::bind2nd(std::plus<uint_fast8_t>(), 1));
      // std::string key = "("+index_to_str(lsmi_key)+")";
      // res[key] = merge_vectors(segments);
      // l_i++;
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
  
  // Generate the segment indices
  std::vector<index_t> get_current_ls_idx(){
    std::vector< index_t > ls_segment_idx = std::vector< index_t >(d);
    std::transform(dim_idx.begin(), dim_idx.end(), ls_segment_idx.begin(), [&](const int d_i){
      return(compute_ls_idx(current_cover.at(d_i), d_i));
    });
    return(ls_segment_idx);
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
  // .method( "get_level_sets", &MultiScale::get_level_sets )
  .method( "get_segment_idx", &MultiScale::get_segment_idx )
  .method( "build_multiscale_configuration", &MultiScale::build_multiscale_configuration )
  .method( "get_segment_map", &MultiScale::get_segment_map )
  .method( "update_multi_cover", &MultiScale::update_multi_cover )
  .method( "test", &MultiScale::test )  
  // .method( "get_ls_that_change", &MultiScale::get_ls_that_change )
  .method( "ls_to_change", &MultiScale::ls_to_change )
  .method( "compute_mapper", &MultiScale::compute_mapper )
  
  ;
}


// Computes Mapper across multiple cover parameterizations. Requires as input:
//  pt_idx := The point indices to change per parameterization
//  from_ls := The level set index the point lies within at the current parameterization
//  to_ls := The level set index the point will intersect at the current parameterization
//  dist_thresh := 1/2th the interval size of each box at the current parameterization
//  clustering_function := the clustering function to apply per parameterization
// [[Rcpp::export]]
void multiscale(){
  // multi_index_t<3> tmp = multi_index_t<3>(3,3,4);
  // std::array<int, 3>::iterator it = tmp.begin(); 
  // for(; it != tmp.end(); ++it){
  // for(auto m: multi_index(1, 2, 3)){
  //   // std::array<int, 3>  m = *it;
  //   Rcout << m[0] << m[1] << m[2] << std::endl; 
  // }
}



/*** R
*/
