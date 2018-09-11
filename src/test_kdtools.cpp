// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
// #include <kdtools.h>// depends(kdtools)
using namespace Rcpp;

template <typename T> 
struct my_object { 
  typedef std::vector<int>::iterator iterator;
  typedef std::vector<int>::const_iterator const_iterator;
  std::vector<int> vec;
  
  my_object(){
    vec = std::vector<int>();
    vec.push_back(1);
    vec.push_back(2);
    // kdtools::kd_sort(std::begin(vec), std::end(vec));
  }

  iterator begin() { return vec.begin(); }
  iterator end() { return vec.end(); }
};

// [[Rcpp::export]]
NumericVector test_it(NumericVector x) {
  my_object<int> v = my_object<int>();
  NumericVector out = NumericVector(5);
  std::vector<int> v2 = { 1, 2, 3, 4, 5 };
  std::copy(std::begin(v), std::end(v), out.begin());
  return(out);
}

/*** R
test_it(42)
*/
