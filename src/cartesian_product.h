#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <cstdlib> // alloca

template <typename T, typename Func> 
inline void CartesianProduct(const std::vector< std::vector<T> >& elems, Func f) {
  
  // Initialize the slots to hold the current iteration value for each depth
  const std::size_t depth = elems.size();
  std::size_t* slots = (std::size_t*) alloca(sizeof(std::size_t) * depth);
  for (std::size_t i = 0; i < depth; i++) { slots[i] = 0; }
  
  // Extract the sizes of each vector in the product
  std::vector<std::size_t> max = std::vector<std::size_t>(depth);
  std::transform(elems.begin(), elems.end(), max.begin(), [](const std::vector<T>& lst){ return(lst.size()); });
  std::vector<T> current_element(depth);
  
  int index = 0, i = 0;
  while (true) {
    
    // Fill the element and apply the lambda 
    i = 0; 
    std::transform(elems.begin(), elems.end(), current_element.begin(), [&i, &slots](const std::vector<T>& v){
      return(v.at(slots[i++]));
    });
    f(current_element);
    
    // Increment
    slots[0]++;
    
    // Carry
    while (slots[index] == max.at(index)) {
      if (index == depth - 1) { return; } // Overflow, we're done
      slots[index++] = 0;
      slots[index]++;
    }
    index = 0;
  }
}

IntegerMatrix make_cartesian_product(const List& vecs);

