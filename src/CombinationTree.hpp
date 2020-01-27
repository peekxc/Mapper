// Combination Tree implementation
// Supports various storage formats + ability to prune entire branches of the tree
// via passed predicate that works on partial combinations. 
// Most code originally from: https://github.com/mraggi/discreture
template < class IntType, class Predicate, class ContainerType = std::vector<IntType>>
class CombinationTree {
public:
    static_assert(std::is_integral<IntType>::value,
                  "Template parameter IntType must be integral");
    static_assert(std::is_signed<IntType>::value,
                  "Template parameter IntType must be signed");
    using value_type = ContainerType;
    using combination = value_type;
    using difference_type = std::ptrdiff_t;
    using size_type = difference_type;
    struct iterator;
    using const_iterator = iterator;

public:
	  IntType n_, k_;
    iterator begin_, end_;
    Predicate pred_;
    
    CombinationTree(IntType n, IntType k, Predicate p)
        : n_(n), k_(k), begin_(n, k, p), end_(p, true), pred_(p) {
    }

    struct iterator  {
        IntType n_{0}, k_{0};
        combination data_{};
        bool at_end_{true};
        Predicate pred_;

        friend class CombinationTree;
    		
        iterator(Predicate p, bool last) : data_(), at_end_(last), pred_(p){}

        iterator(IntType n, IntType k, Predicate p)
            : n_(n), k_(k), data_(), at_end_(false), pred_(p)
        {
            data_.reserve(k_);

            while (DFSUtil(data_, pred_, n_, k_)) {
              if (data_.size() == static_cast<size_t>(k_)) {
                return;
              }
            }
            at_end_ = true;
        }

        inline bool is_at_end(IntType n) const {
        	return at_end_;
        }
        
        iterator& operator++(){
          while (DFSUtil(data_, pred_, n_, k_))  {
            if (data_.size() == static_cast< size_t >(k_))
              return(*this);
          }
          at_end_ = true;
        	return(*this);
        }
	 			const combination& operator*() const { return data_; }
       	
    		bool operator!=(const iterator& it) const {
    			return(!(*this == it));
    		}
    		bool operator==(const iterator& it) const {
            if (at_end_ != it.at_end_)
                return false;
            if (at_end_)
                return true;
            return data_ == it.data_;
        }

    }; // end class iterator

    const iterator& begin() const { return begin_; }
    const iterator& end() const { return end_; }

private:

    static bool augment(combination& comb,
                        Predicate pred,
                        IntType n_,
                        IntType k_,
                        IntType start = 0)
    {
        if (comb.empty()) {
          if (start < n_ - k_ + 1) {
            comb.push_back(start);
            return true;
          }
          return false;
        }

        auto last = comb.back();
        auto guysleft = k_ - comb.size();
        start = std::max(static_cast<IntType>(last + 1), start);

        for (size_t i = start; i < n_ - guysleft + 1; ++i) {
            comb.push_back(i);
            if (pred(comb))
                return true;
            comb.pop_back();
        }

        return false;
    }

    static bool DFSUtil(combination& comb, Predicate pred, IntType n_, IntType k_){
      if (comb.size() < k_) {
          if (augment(comb, pred, n_, k_))
              return true;
      }

      auto last = comb.back();

      // If it can't be augmented, be it because size is already k or else, we
      // have to start backtracking
      while (!comb.empty()) {
          last = comb.back();
          comb.pop_back();
          if (augment(comb, pred, n_, k_, last + 1))
              return true;
      }
      return false;
  	}
}; // end class CombinationTree

template < class PartialPredicate, class idx_t = int32_t >
auto PrunedCombinationTree(const idx_t N, const idx_t K, PartialPredicate p) -> CombinationTree< idx_t, PartialPredicate > {
	return CombinationTree< idx_t, PartialPredicate >(N, K, p);
}
