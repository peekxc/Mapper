// utility_rcpp.hpp
// Contains pure-template implementations suitable for being inlined.  

// Applies the function Func to all pairwise combinations in the range [first, last)
template<typename Iter, typename Func>
inline void combine_pairwise(Iter first, Iter last, Func func)
{
  for(; first != last; ++first){
    for(Iter next = std::next(first); next != last; ++next){
      func(*first, *next);
    }
  }
}

// Implements a generic n-vector cartesian product for a vector of vectors using an iterative design pattern. 
// Individual items are put into a fixed-sized vector and given to a passed in Lambda function. Limited to products 
// of vectors having the same type. Original design based on: 
// https://stackoverflow.com/questions/18732974/c-dynamic-number-of-nested-for-loops-without-recursion/30808351
template <typename T, typename Func> 
inline void CartesianProduct(const std::vector< std::vector<T> >& elems, Func&& f) {
  
  // Initialize the slots to hold the current iteration value for each depth
  const std::size_t depth = elems.size();
  std::size_t* slots = (std::size_t*) alloca(sizeof(std::size_t) * depth);
  for (std::size_t i = 0; i < depth; i++) { slots[i] = 0; }
  
  // Extract the sizes of each vector in the product
  std::vector<std::size_t> max = std::vector<std::size_t>(depth);
  std::transform(elems.begin(), elems.end(), max.begin(), [](const std::vector<T>& lst){ return(lst.size()); });
  std::vector<T> current_element(depth);
  
  std::size_t index = 0, i = 0;
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