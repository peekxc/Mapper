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
