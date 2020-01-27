// nerve_utility.cpp
// Contains utility functions related to the nerve construction

#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <algorithm>

using std::vector; 
using std::begin; 
using std::end;

// Moves through two ordered sets to see if they are disjoint
template <typename SetA, typename SetB>
bool sorted_disjoint(const SetA &a, const SetB &b) {
  auto it_a = a.begin();
  auto it_b = b.begin();
  while (it_a != a.end() && it_b != b.end()) {
    switch (*it_a == *it_b ? 0 : *it_a < *it_b ? -1 : 1) {
    case 0:
      return false;
    case -1:
      it_a = std::lower_bound(++it_a, a.end(), *it_b);
      break;
    case 1:
      it_b = std::lower_bound(++it_b, b.end(), *it_a);
      break;
    }
  }
  return true;
}

// Given two random-access iterator ranges, (a1, a2), (b1, b2), return a boolean indicating 
// whether or not they have a non-empty intersection. Does not assumes either is sorted.
template <typename Iter>
bool nonempty(Iter a1, Iter a2, Iter b1, Iter b2){
  using it_cat = typename std::iterator_traits<Iter>::iterator_category;
  static_assert(std::is_same<it_cat, std::random_access_iterator_tag>::value, "Iterator type must be random-access."); 
  using T = typename std::iterator_traits<Iter>::value_type;
    
  // Either empty == they do not have an intersection
  const size_t a_sz = std::distance(a1, a2); 
  const size_t b_sz = std::distance(b1, b2); 
  if (a_sz == 0 || b_sz == 0) { return false; }
  
  // a is much smaller than b => partial sort b, then do binary search on b for each element of a
  if (a_sz * 100 < b_sz) {
    vector<T> b_sort(b_sz);
    std::partial_sort_copy(b1, b2, begin(b_sort), end(b_sort)); // partial-sorted elements of y copied to y_sort
    while (a1 != a2){
      if (std::binary_search(begin(b_sort), end(b_sort), T(*a1))) { return(true); }
      ++a1;
    }
    return(false);
  } else if (b_sz * 100 < a_sz) { // Opposite case
    vector<T> a_sort(a_sz);
    std::partial_sort_copy(a1, a2, begin(a_sort), end(a_sort)); // partial-sorted elements of y copied to y_sort
    while (b1 != b2){
      if (std::binary_search(begin(a_sort), end(a_sort), T(*b1))) { return(true); }
      ++b1;
    }
    return(false);
  }
  
  // Otherwise, sort both, then use lower_bound type approach to potentially skip massive sections.
  vector<T> a_sort(a_sz), b_sort(b_sz);
  std::partial_sort_copy(a1, a2, begin(a_sort), end(a_sort)); // partial-sorted elements of y copied to y_sort
  std::partial_sort_copy(b1, b2, begin(b_sort), end(b_sort)); // partial-sorted elements of y copied to y_sort
  return !sorted_disjoint(a_sort, b_sort);
}

template <typename Iter>
bool nfold_nonempty(vector< std::pair< Iter, Iter > > ranges){
  using it_cat = typename std::iterator_traits<Iter>::iterator_category;
  static_assert(std::is_same<it_cat, std::random_access_iterator_tag>::value, "Iterator type must be random-access.");
  using T = typename std::iterator_traits<Iter>::value_type;

  // Either empty == they do not have an intersection
  const size_t n_rng = ranges.size();
  vector< size_t > rng_sizes = vector< size_t >(n_rng);
  std::transform(begin(ranges), end(ranges), begin(rng_sizes), [](const std::pair<Iter, Iter> rng){
    return std::distance(rng.first, rng.second);
  });
  bool any_empty = std::any_of(begin(rng_sizes), end(rng_sizes), [](const size_t sz){ return sz == 0; });
  if (any_empty){ return(false); };

  // Use multiset to track number of ids
  std::multiset<T> ids; 
  const auto insert_rng = [&ids](Iter a, Iter b){ 
    std::for_each(a, b, [&ids](T elem){ ids.insert(elem); }); 
  };
  for (auto rng: ranges){ insert_rng(rng.first, rng.second); };
  
  // If any ids appeared k times, there's a nonempty intersection between them all 
  bool nonempty_int = std::any_of(ids.begin(), ids.end(), [&ids, &n_rng](const T id){
    return (ids.count(id) == n_rng);
  });
  return nonempty_int;
}

template <typename Iter>
auto intersection(Iter a1, Iter a2, Iter b1, Iter b2) 
  -> vector< typename std::iterator_traits<Iter>::value_type > {
  using T = typename std::iterator_traits<Iter>::value_type;
  std::unordered_set< T > c;
  for (; a1 != a2; ++a1) {
    if (std::find(b1, b2, *a1) != b2) { 
      c.insert(*a1); 
    }
  }
  vector< T > res(begin(c), end(c));
  return(res);
}


// auto iv_intersection(const IntegerVector::const_iterator a1, const IntegerVector::const_iterator a2, 
//                   const IntegerVector::const_iterator b1, const IntegerVector::const_iterator b2) 
//   -> vector< int > {
//   std::unordered_set< int > a = std::unordered_set< int >(a1, a2);
//   std::unordered_set< int > b = std::unordered_set< int >(b1, b2);
//   std::unordered_set< int > c;
//   for (auto i = a.begin(); i != a.end(); i++) {
//     if (b.find(*i) != b.end()) { 
//       c.insert(*i); 
//     }
//   }
//   vector< int > res(begin(c), end(c));
//   return(res);
// }

template <typename Iter>
vector< int > nfold_intersection(vector< std::pair< Iter, Iter > > ranges){
  using it_cat = typename std::iterator_traits<Iter>::iterator_category;
  static_assert(std::is_same<it_cat, std::random_access_iterator_tag>::value, "Iterator type must be random-access.");
  using T = typename std::iterator_traits<Iter>::value_type;
  vector< T > res;
  const size_t n_pairs = ranges.size()-1; 
  for (size_t i = 0; i < n_pairs; ++i){
    const std::pair< Iter, Iter > rng1 = ranges[i];
    const std::pair< Iter, Iter > rng2 = ranges[i+1];
    vector< T > c_int = intersection(rng1.first, rng1.second, rng2.first, rng2.second);
    res.insert(res.end(), begin(c_int), end(c_int));
    if (res.size() == 0){
      return(res);
    }
  }
  return(res);
}

// IntegerVector test_intersection1(const IntegerVector v1, const IntegerVector v2){
//   auto res = intersection(v1.begin(), v1.end(), v2.begin(), v2.end());
//   return(wrap(res));
// }
// IntegerVector test_intersection2(const IntegerVector v1, const IntegerVector v2){
//   IntegerVector res = Rcpp::intersect(v1, v2);
//   return(res);
// }
// LogicalVector test_nonempty(const IntegerVector v1, const IntegerVector v2){
//   const bool is_non_empty = nonempty_intersection(v1.begin(), v1.end(), v2.begin(), v2.end());
//   return(is_non_empty);
// }

/*** R
# x <- sample(x = 1e5, size = 1e4, replace = FALSE)
# y <- sample(x = 1e5, size = 1e4, replace = FALSE)
# 
# microbenchmark::microbenchmark(intersect(x, y))
# microbenchmark::microbenchmark(test_intersection1(x, y))
# microbenchmark::microbenchmark(test_intersection2(x, y))
# microbenchmark::microbenchmark(test_nonempty(x, y))
# 
# test_nonempty(c(1, 2, 3), c(4, 5, 6))
*/
