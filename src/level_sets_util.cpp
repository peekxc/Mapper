#include <Rcpp.h>
using namespace Rcpp;

// Given 
// Creates a list of elements, one per block, where each element is itself a list containing: 
// 1) lsfi := the level set flat index of the level set to update 
// 2) pt_idx := the points to add to the current level set, in the current block 
// Note that this is particularly tricky to handle due to non-orthogonal level sets  
// [[Rcpp::export]]
NumericVector createUpdateBlocks(const List& G, const List& overlap_blocks, const int n_lvl_sets, const int n_blocks, const int n) {
  const int d = G.size(); 
  //List res1 = List(n_blocks);
  std::vector< std::vector< IntegerVector > > res = std::vector< std::vector< IntegerVector > >(n_blocks);
  // List res = List(n_lvl_sets); // list to store the resulting point indices to update  
  std::vector< List > G_ = as< std::vector< List > >(G);
  for (int i = 0; i < n; ++i){
    for (int d_i = 0; d_i < d; ++d_i){
      const List overlap_info = G[d_i];
      const IntegerMatrix& block_ids = overlap_blocks[d_i];
      const IntegerMatrix& target_lsfi = overlap_info["target_lsfi"];
      const int n = block_ids.nrow(); 
    }
    
    // IntegerMatrix::ConstRow lvl_sets_to_intersect = target_lsfi(i, _);
    // IntegerMatrix::ConstRow::const_iterator pt_lsfi = target_lsfi(i, _).begin();
    // IntegerMatrix::ConstRow::const_iterator pt_block_id = block_ids(i, _).begin();
    // for (; pt_lsfi != target_lsfi(i, _).end(); ++pt_lsfi, ++pt_block_id){
    //   const int pt_target_lsfi = *pt_lsfi;
    //   res[*pt_block_id][(*pt_lsfi)].push_back(i);
    // }
  }
  return(NumericVector::create());
}


/*** R
# timesTwo(42)
*/
