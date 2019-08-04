//----------------------------------------------------------------------
//                        Disjoint-set data structure 
// File:                        union_find.cpp
//----------------------------------------------------------------------
// Copyright (c) 2018 Matt Piekenbrock. All Rights Reserved.
//
// Class definition based off of data-structure described here:  
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure

#include "UnionFind.h"

UnionFind::UnionFind(const size_t _size) : size(_size), parent(_size), rank(_size)
{
  if (_size <= 0){ stop("size must be positive."); }
  for (size_t i = 0; i < size; ++i)
  { parent[i] = i, rank[i] = 0; }
}
  
// Destructor not needed w/o dynamic allocation
UnionFind::~UnionFind() { }
  
SEXP UnionFind::as_XPtr(){
  Rcpp::XPtr< UnionFind> p(this, false); // don't register the finalizer to allow passing to other methods
  return(p);
}

// Union x and y
void UnionFind::Union(const size_t x, const size_t y) {
  if (x >= size){ stop("x out of range"); }
  if (y >= size){ stop("y out of range"); }
  const size_t xRoot = Find(x);
  const size_t yRoot = Find(y);
  if (xRoot == yRoot)
   return; 
  else if (rank[xRoot] > rank[yRoot])
    parent[yRoot] = xRoot; 
  else if (rank[xRoot] < rank[yRoot]) 
    parent[xRoot] = yRoot; 
  else if (rank[xRoot] == rank[yRoot])
  {
    parent[yRoot] = parent[xRoot];
    rank[xRoot] = rank[xRoot] + 1;
  }
}

// Union all indices in the the integer vector
// void UnionFind::UnionAll(const IntegerVector idx){
//   if (idx.size() <= 1){ return; }
//   const int n_pairs = idx.size()-1;
//   for (int i = 0; i < n_pairs; ++i){ Union(idx[i], idx[i+1]); }
// }
void UnionFind::UnionAll(const vector< size_t >& idx){
  if (idx.size() <= 1){ return; }
  const size_t n_pairs = idx.size()-1;
  for (size_t i = 0; i < n_pairs; ++i){ Union(idx[i], idx[i+1]); }
}

// IntegerVector UnionFind::FindAll(const IntegerVector idx){
//   if (idx.size() == 0){ return IntegerVector::create(); }
//   const int n = idx.size();
//   IntegerVector cc = Rcpp::no_init(n);
//   std::transform(idx.begin(), idx.end(), cc.begin(), [=](const int i){
//     return(Find(i));
//   });
//   return(cc);
// }

vector< size_t > UnionFind::FindAll(const vector< size_t >& idx){
  using idx_v = vector< size_t >;
  if (idx.size() == 0){ return idx_v(); }
  const size_t n = idx.size();
  idx_v cc = idx_v(n);
  std::transform(idx.begin(), idx.end(), cc.begin(), [=](const size_t i){
    return(Find(i));
  });
  return(cc);
}

// Find operation
const size_t UnionFind::Find(const size_t x) {
  if (x < 0 || x >= size){ return(-1); }
  if (parent[x] == x){ return x; } else {
    parent[x] = Find(parent[x]);
    return parent[x];
  }
}

// Add more sets to the UF
void UnionFind::AddSets(const size_t n_sets){
  parent.resize(size + n_sets);
  std::iota(parent.begin() + size, parent.end(), size); // parent initialized incrementally
  rank.resize(size + n_sets, 0); // rank all 0 
  size += n_sets; 
}


// Return new integer vector representing the connected components
IntegerVector UnionFind::getCC(){
  IntegerVector cc = Rcpp::no_init(size);
  for (size_t i = 0; i < size; ++i){ cc[i] = Find(i); }
  return(cc);
}

// Simple method to print the CCs on one line
void UnionFind::printCC(){
  for (size_t j = 0; j < size; ++j){ Rcout << Find(j) << " "; }
  Rcout << std::endl;
}

// Export as an Rcpp module
RCPP_MODULE(union_find_module) {
  Rcpp::class_<UnionFind>("UnionFind")
  .constructor< std::size_t >()
  .field_readonly( "size", &UnionFind::size)
  .field_readonly( "parent", &UnionFind::parent)
  .field_readonly( "rank", &UnionFind::rank)
  .method("as_XPtr", &UnionFind::as_XPtr)
  .method("print", &UnionFind::printCC )
  .method("connected_components", &UnionFind::getCC)
  .method("find", &UnionFind::Find)
  .method("find_all", &UnionFind::FindAll)
  .method("union", &UnionFind::Union)
  .method("union_all", &UnionFind::UnionAll)
  .method("add_sets", &UnionFind::AddSets)
  ;
}