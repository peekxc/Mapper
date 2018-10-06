// utility.h 
// Set of inline utility functions.


inline std::size_t index_lower_triangular(std::size_t from, std::size_t to, const std::size_t N){
  if (from < to){ std::swap(from, to); }
  return((N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1));
}
#ifndef INDEX_TO
  #define INDEX_TO(k, n) n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5) // expects 0-based, returns 0-based
#endif
#ifndef INDEX_FROM
  #define INDEX_FROM(k, n, i) k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2 // expects 0-based, returns 0-based
#endif


// Creates a vector with the range [i, j]
template <typename T> 
std::vector<T> seq_ij(const int i, const int j){
  std::size_t sz = std::abs(j - i)+1;
  std::vector<T> rng = std::vector<T>(sz);
  std::iota(rng.begin(), rng.end(), i);
  return(rng);
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

template<typename ForwardIterator>
inline std::map<int, int> get_unique_indices(ForwardIterator first, ForwardIterator last){
  std::map<int, int> pt_to_unique_idx;
  for(std::size_t i = 0; first != last; ++i, ++first){
    auto it = pt_to_unique_idx.find(*first);
    if (it == pt_to_unique_idx.end()) { // value doesn't exist
      pt_to_unique_idx.emplace(*first, i);
    }
  }
  return pt_to_unique_idx;
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
