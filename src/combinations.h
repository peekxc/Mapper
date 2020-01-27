#include <type_traits>
#include <vector>
#include <functional>
#include <iterator>


template<typename T>
using it_diff_t = typename std::iterator_traits<T>::difference_type;

using std::begin;
using std::end; 
using std::vector; 
using std::size_t;

// Rotates two discontinuous ranges to put *first2 where *first1 is.
//     If last1 == first2 this would be equivalent to rotate(first1, first2, last2),
//     but instead the rotate "jumps" over the discontinuity [last1, first2) -
//     which need not be a valid range.
//     In order to make it faster, the length of [first1, last1) is passed in as d1,
//     and d2 must be the length of [first2, last2).
//  In a perfect world the d1 > d2 case would have used swap_ranges and
//     reverse_iterator, but reverse_iterator is too inefficient.
template <class It>
void rotate_discontinuous(
	It first1, It last1, it_diff_t< It > d1,
  It first2, It last2, it_diff_t< It > d2)
{
  using std::swap;
  if (d1 <= d2){ std::rotate(first2, std::swap_ranges(first1, last1, first2), last2); }
  else {
		It i1 = last1;
		while (first2 != last2)
			swap(*--i1, *--last2);
		std::rotate(first1, i1, last1);
  }
}

// Call f() for each combination of the elements [first1, last1) + [first2, last2)
//    swapped/rotated into the range [first1, last1).  As long as f() returns
//    false, continue for every combination and then return [first1, last1) and
//    [first2, last2) to their original state.  If f() returns true, return
//    immediately.
//  Does the absolute mininum amount of swapping to accomplish its task.
//  If f() always returns false it will be called (d1+d2)!/(d1!*d2!) times.
template < typename It, typename Lambda >
bool combine_discontinuous(
	It first1, It last1, it_diff_t< It > d1,  
 	It first2, It last2, it_diff_t< It > d2,
	Lambda&& f, it_diff_t< It > d = 0)
{
	using D = it_diff_t< It >;
	using std::swap;
	if (d1 == 0 || d2 == 0){ return f(); }
	if (d1 == 1) {
		for (It i2 = first2; i2 != last2; ++i2) {
			if (f()){ return true; }
			swap(*first1, *i2);
		}
	}
	else {
		It f1p = std::next(first1), i2 = first2;
		for (D d22 = d2; i2 != last2; ++i2, --d22){
			if (combine_discontinuous(f1p, last1, d1-1, i2, last2, d22, f, d+1))
				return true;
			swap(*first1, *i2);
		}
	}
	if (f()){ return true; }
	if (d != 0){ rotate_discontinuous(first1, last1, d1, std::next(first2), last2, d2-1); }
	else { rotate_discontinuous(first1, last1, d1, first2, last2, d2); }
	return false;
}


template < typename Lambda, typename It > 
struct bound_range { 
	Lambda f_;
	It first_, last_;
	bound_range(Lambda& f, It first, It last) : f_(f), first_(first), last_(last) {}
	bool operator()(){ return f_(first_, last_); } 
	bool operator()(It, It) { return f_(first_, last_); }
};

template <class It, class Function>
Function for_each_combination(It first, It mid, It last, Function&& f) {
	bound_range<Function&, It> wfunc(f, first, mid);
	combine_discontinuous(first, mid, std::distance(first, mid),
												mid, last, std::distance(mid, last),
												wfunc);
	return std::move(f);
}

template < class I, class Function >
void for_each_combination(I n, I k, Function&& f) {
  static_assert(std::is_integral<I>::value, "Must be integral type.");
  using It = typename vector< I >::iterator;
  vector< I > seq_n(n);
  std::iota(begin(seq_n), end(seq_n), 0);
  for_each_combination(begin(seq_n), begin(seq_n)+k, end(seq_n), [&f](It a, It b){
    return f(a, b);
  });
  return;
}


template < class I, class Function >
void for_each_combination_idx(I n, I k, Function&& f) {
  static_assert(std::is_integral<I>::value, "Must be integral type.");
  using It = typename vector< I >::iterator;
  vector< I > seq_n(n);
  std::iota(begin(seq_n), end(seq_n), 0);
  for_each_combination(begin(seq_n), begin(seq_n)+k, end(seq_n), [k, &f](It a, It b){
    vector< I > cc(k);
    std::transform(a, b, begin(cc), [](const I num){ return num; });
    f(cc);
    return false; 
  });
  return;
}

