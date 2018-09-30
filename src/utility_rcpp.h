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

List resize_list(const List& x, int n){
  const std::size_t lst_sz = x.size();
  List y(n);
  for(int i = 0; i < lst_sz; i++) { y[i] = x[i]; }
  return(y);
}

