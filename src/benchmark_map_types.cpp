#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <unordered_map>
#include <map>
#include <random>

typedef std::vector< uint_fast8_t > index_t;

std::string index_to_str2(index_t idx){
  std::string res_key{ idx.begin(), idx.end() };
  return(res_key);
}

void bench_map(const IntegerMatrix& keys, const IntegerVector& resolution) {
  
  std::map<index_t, std::size_t> ls_map; 
  std::unordered_map<std::string, std::size_t> ls_map2;

  std::size_t n_keys = keys.nrow(); 
  index_t c_key = index_t(keys.ncol());
  for (std::size_t i = 0; i < n_keys; ++i){
    IntegerVector key_vec = keys.row(i);
    std::transform(key_vec.begin(), key_vec.end(), c_key.begin(), [](const int idx){ return(static_cast<uint_fast8_t>(idx)); });
    ls_map.emplace(c_key, i);
    ls_map2.emplace(index_to_str2(c_key), i);
  }

  
  // Run the benchmark
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  auto gen_key = [&resolution, &c_key, &distribution, &generator](){
    std::transform(resolution.begin(), resolution.end(), c_key.begin(), [&distribution, &generator](const int k_i){
      return(static_cast<uint_fast8_t>(distribution(generator) * (k_i)) + 1);
    });
  };



  const int n = 500000;
  const int d = resolution.length();
  // std::vector<std::size_t> res = std::vector<std::size_t>(n);

  auto start = std::chrono::system_clock::now();
  for (std::size_t i = 0; i < n; ++i){
    gen_key();
    // for (int d_i = 0; d_i < d; ++d_i){ Rcout << static_cast<int>(c_key[d_i]) << " "; }
    // Rcout << std::endl; 
    std::size_t testing = ls_map.at(c_key);
  }
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  Rcout << elapsed.count() << "ms" << std::endl; 
  
  std::string umap_key; 
  start = std::chrono::system_clock::now();
  for (std::size_t i = 0; i < n; ++i){
    gen_key();
    // for (int d_i = 0; d_i < d; ++d_i){ Rcout << static_cast<int>(c_key[d_i]) << " "; }
    // Rcout << std::endl; 
    // umap_key = { c_key.begin(), c_key.end() };
    umap_key.assign(c_key.begin(), c_key.end());
    std::size_t testing = ls_map2.at(umap_key);
  }
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  Rcout << elapsed.count() << "ms" << std::endl; 
  
  // TODO: test just raw computational version using accumulate w/ reverse iterators
  // start = std::chrono::system_clock::now();
  // for (std::size_t i = 0; i < n; ++i){
  //   gen_key();
  //   
  //   std::size_t testing = c_key
  // }
  // end = std::chrono::system_clock::now();
  // elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  // Rcout << elapsed.count() << "ms" << std::endl; 
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
resolution <- rep(7L, 3L)
index_set <- as.matrix(do.call(expand.grid, lapply(resolution, seq)))
bench_map(index_set, resolution)
*/
