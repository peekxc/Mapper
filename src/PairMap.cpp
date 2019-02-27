#include <Rcpp.h>
using namespace Rcpp;

#include <cstdint>
namespace IntSize {
  template <class > struct next_size;
  template <class T> using next_size_t = typename next_size<T>::type;
  template <class T> struct tag { using type = T; };
  
  template <> struct next_size< std::uint8_t >  : tag< std::uint16_t > { };
  template <> struct next_size< std::uint16_t > : tag< std::uint32_t > { };
  template <> struct next_size< std::uint32_t > : tag< std::uint64_t > { };
}

// PairMap - Stores key-value pairs contiguous where the 'keys' are 
// represented by a pair of integral-types. Uses a pairing function to 
// 'hash' the key into a single integer. Keys are stored in a sorted, contiguous manner, 
// which them act as index maps into the values array. 
// This data structure emphasizes efficient access time and iteration at the expense of
// (potentially) inefficient insertion and deletion. 
template < typename I, typename V >
struct PairMap {
  static_assert(std::is_integral<I>::value, "Integral key type required.");
  static_assert(std::is_unsigned<I>::value, "Integral key type must be unsigned.");
  using K = IntSize::next_size_t<I>; // Guarenteed to be at least twice as large as I
  std::vector< K > keys;
  std::vector< V > values; 
  

  void insert(I a, I b, V val){
    
  }
  
  
private: 
  K pair(const I x, const I y) const {  return(x >= y ? x * x + x + y : x + y * y); }
  std::pair<I, I> unpair (const K z){
    const K z_floor = std::floor(std::sqrt(z)); 
    const K z_upper = std::pow(z_floor, 2);
    return(
      z - z_upper < z_floor ? 
      std::make_pair(z - z_upper, z_floor) : 
      std::make_pair(z_floor, z - z_upper - z_floor)
    );
  }
};

// TODO 
NumericVector test_pair_map(NumericVector x) {
  return x * 2;
}


/*** R

*/
