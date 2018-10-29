//----------------------------------------------------------------------
//                        Disjoint-set data structure 
// File:                        union_find.cpp
//----------------------------------------------------------------------
// Copyright (c) 2018 Matt Piekenbrock. All Rights Reserved.
//
// Class definition based off of data-structure described here:  
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure

#include "UnionFind.h"

UnionFind::UnionFind(const int _size) : size(_size), parent(_size), rank(_size)
{
  if (_size <= 0){ stop("size must be positive."); }
  for (int i = 0; i < size; ++i)
  { parent[i] = i, rank[i] = 0; }
}
  
// Destructor not needed w/o dynamic allocation
UnionFind::~UnionFind() { }
  
SEXP UnionFind::as_XPtr(){
  Rcpp::XPtr< UnionFind> p(this, false); // don't register the finalizer to allow passing to other methods
  return(p);
}

  
void UnionFind::Union(const int x, const int y) {
  if (x < 0 || x >= size){ stop("x out of range"); }
  if (y < 0 || y >= size){ stop("y out of range"); }
  const int xRoot = Find(x);
  const int yRoot = Find(y);
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

void UnionFind::UnionAll(const IntegerVector idx){
  if (idx.size() <= 1){ return; }
  // IntegerVector::const_iterator it = idx.begin(); 
  const int n = idx.size();
  for (int i = 0; i < n-1; ++i){
    Union(idx[i], idx[i+1]);
  }
}

IntegerVector UnionFind::FindAll(const IntegerVector idx){
  if (idx.size() == 0){ return IntegerVector::create(); }
  const int n = idx.size();
  IntegerVector cc = Rcpp::no_init(n);
  for (int i = 0; i < n; ++i){
    cc[i] = Find(idx[i]);
  }
  return(cc);
}

const int UnionFind::Find(const int x) {
  if (x < 0 || x >= size){ return(-1); }
  if (parent[x] == x)
    return x; 
  else
  {
    parent[x] = Find(parent[x]);
    return parent[x];
  }
}

// Return new integer vector representing the connected components
IntegerVector UnionFind::getCC(){
  IntegerVector cc = Rcpp::no_init(size);
  for (unsigned int i = 0; i < size; ++i){ cc[i] = Find(i); }
  return(cc);
}

// Simple method to print the CCs on one line
void UnionFind::printCC(){
  for (int j = 0; j < size; ++j){ Rcout << Find(j) << " "; }
  Rcout << std::endl;
}

// Export as an Rcpp module
RCPP_MODULE(union_find_module) {
  Rcpp::class_<UnionFind>("UnionFind")
  .constructor<const int>()
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
  ;
}