// utility_rcpp.h
// Utility functions using Rcpp. 

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]   


// Nice simple interface from: https://stackoverflow.com/questions/17973442/index-element-from-list-in-rcpp
template <typename WHAT>
class ListOf : public List {
public:
  template <typename T>
  ListOf( const T& x) : List(x){}
  WHAT operator[](int i){ return as<WHAT>( ( (List*)this)->operator[]( i) ) ; }
  WHAT at(int i){ return as<WHAT>( ( (List*)this)->at(i) ) ; }
};


template <typename T> 
static inline IntegerVector to_ivec(std::vector<T> vec){
  static_assert(std::is_integral<T>::value, "T must be integral type");
  IntegerVector res = IntegerVector(vec.size());
  std::size_t i = 0; 
  std::for_each(vec.begin(), vec.end(), [&i, &res](const T val){
    res.at(i++) = static_cast<int>(val);
  });
  return(res);
}

// Fast partial-sort/binary-search check to see if the intersection between two given vectors has non-zero length
// Loosely based off of bugged version found here: https://stackoverflow.com/questions/21359432/a-c-version-of-the-in-operator-in-r
static inline bool any_is_in(const IntegerVector& x, const IntegerVector& y){
  std::vector<int> y_sort(y.size());
  std::partial_sort_copy (y.begin(), y.end(), y_sort.begin(), y_sort.end()); // partial-sorted elements of y copied to y_sort
  const int nx = x.size();
  for (int i = 0; i < nx; ++i) {
    if (std::binary_search(y_sort.begin(), y_sort.end(), x[i])) {
      return(true); // end the search
    }
  }
  return(false);
}

// rbindlist_int: Takes a list of integer vectors and rbind's them together.
static inline IntegerMatrix rbindlist_int(std::list<IntegerVector>& lst){
  unsigned int n = lst.size();
  if(n == 0) { Rcpp::stop("Invalid sized list."); }
  unsigned int d = lst.front().size();
  Rcpp::IntegerMatrix res = Rcpp::no_init(n, d);
  size_t i = 0;
  for (std::list<IntegerVector>::iterator it = lst.begin(); it != lst.end(); ++it, ++i) {
    if ((*it).size() != d) { Rcpp::stop("Invalid integer vector size."); }
    res(i, _) = *it;
  }
  return res;
}

static inline List resize_list(const List& x, int n){
  const std::size_t lst_sz = x.size();
  List y(n);
  for(int i = 0; i < lst_sz; i++) { y[i] = x[i]; }
  return(y);
}

