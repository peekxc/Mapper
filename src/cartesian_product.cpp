
#include "cartesian_product.h"

// Implements a generic n-vector cartesian product for a vector of vectors. Individual items are put into a fixed-sized vector and given to
// a passed in Lambda function. Limited to products of vectors having the same type. Original design based loosely on: 
// https://stackoverflow.com/questions/18732974/c-dynamic-number-of-nested-for-loops-without-recursion/30808351
// template <typename T, typename Func> 
// void CartesianProduct(const std::vector< std::vector<T> >& elems, Func f) {
//   
//   // Initialize the slots to hold the current iteration value for each depth
//   const std::size_t depth = elems.size();
//   std::size_t* slots = (std::size_t*) alloca(sizeof(std::size_t) * depth);
//   for (std::size_t i = 0; i < depth; i++) { slots[i] = 0; }
//   
//   // Extract the sizes of each vector in the product
//   std::vector<std::size_t> max = std::vector<std::size_t>(depth);
//   std::transform(elems.begin(), elems.end(), max.begin(), [](const std::vector<T>& lst){ return(lst.size()); });
//   std::vector<T> current_element(depth);
//   
//   int index = 0, i = 0;
//   while (true) {
//     
//     // Fill the element and apply the lambda 
//     i = 0; 
//     std::transform(elems.begin(), elems.end(), current_element.begin(), [&i, &slots](const std::vector<T>& v){
//       return(v.at(slots[i++]));
//     });
//     f(current_element);
//     
//     // Increment
//     slots[0]++;
//     
//     // Carry
//     while (slots[index] == max.at(index)) {
//       if (index == depth - 1) { return; } // Overflow, we're done
//       slots[index++] = 0;
//       slots[index]++;
//     }
//     index = 0;
//   }
// }

// [[Rcpp::export]]
IntegerMatrix make_cartesian_product(const List& vecs){
  std::vector< std::vector<int> > vv(vecs.size());
  std::size_t n_rows = 1;
  for (int i = 0; i < vecs.size(); ++i){
    IntegerVector v = vecs.at(i);
    vv.at(i) = as< std::vector<int> >(v);
    n_rows = n_rows * v.size();
  }

  // Copy to result
  int i = 0;
  IntegerMatrix res = IntegerMatrix(n_rows, vecs.size());
  CartesianProduct(vv, [&i, &res](std::vector<int> idx){
    res(i++, _) = as<IntegerVector>(wrap(idx));
  });
  return(res);
}

// void test_cart_product(const List& test) {
//   const int d = test.size();
//   std::vector< std::vector<int> > elements(d);
//   for (int d_i = 0; d_i < d; ++d_i){
//     IntegerVector vec = test.at(d_i);
//     std::vector<int> dim_values = as< std::vector<int> >(vec);
//     elements.at(d_i) = dim_values;
//   }
//   
//   CartesianProduct(elements, [](const std::vector<int> tuple){
//     Rcout << tuple.at(0) << tuple.at(1) << std::endl;
//   });
//   
// }


/*** R
# timesTwo(42)
*/
