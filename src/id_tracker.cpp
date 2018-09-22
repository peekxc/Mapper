#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]   

// Very simple class that allows tracking the allocation/deallocation of sequentiall contiguous indices. 
// That is, upon removing ids, newly requested ids are guarenteed to fall along the first 
struct ID_Generator {
  std::vector<bool> id_vec; 
  ID_Generator() : id_vec(std::vector<bool>()){ }
  IntegerVector get_new_ids(const std::size_t n_ids){
    if (n_ids == 0){ return(IntegerVector::create()); }
    if (id_vec.size() == 0){ 
      id_vec.resize(n_ids, false);
      return(Rcpp::seq(1, n_ids)); 
    }
    if (n_free() < n_ids){  id_vec.resize(id_vec.size() + n_ids, true); }
    IntegerVector res = IntegerVector(); 
    std::vector<std::size_t> idx_removed = std::vector<std::size_t>();
    std::size_t n_to_get = n_ids; 
    
    // Loop through, finding the first indices (1-based) that are true
    for(std::size_t cc = 0; n_to_get > 0 && cc < id_vec.size(); ++cc){
      if (id_vec.at(cc)){ 
        idx_removed.push_back(cc);
        res.push_back(cc+1);
        n_to_get--;
      }
    }
    std::for_each(idx_removed.begin(), idx_removed.end(), [&](const int i){ id_vec.at(i) = false; });
    return(res);
  }
  void remove_ids(const IntegerVector& ids_to_remove){
    std::for_each(ids_to_remove.begin(), ids_to_remove.end(), [&](const int i){
      id_vec.at(i-1) = true;
    });
  }
  std::size_t n_free(){
    return(std::count(id_vec.begin(), id_vec.end(), true));
  }
  void print_available(){
    std::for_each(id_vec.begin(), id_vec.end(), [](const bool id_status){
      Rcout << static_cast<int>(id_status);
    });
    Rcout << std::endl; 
  }
  
};

RCPP_MODULE(id_tracker_module) {
  Rcpp::class_<ID_Generator>("ID_Generator")
  .constructor()
  .method( "get_new_ids", &ID_Generator::get_new_ids )
  .method( "remove_ids", &ID_Generator::remove_ids )
  .method( "print_available", &ID_Generator::print_available )
  ;  
}


/*** R
id_gen <- ID_Generator$new()
all(id_gen$get_new_ids(3) == c(1, 2, 3))
all(id_gen$get_new_ids(3) == c(4, 5, 6))
id_gen$print_available()
id_gen$remove_ids(1:3)
id_gen$print_available()
*/
