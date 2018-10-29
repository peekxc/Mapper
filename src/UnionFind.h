//----------------------------------------------------------------------
//                        Disjoint-set data structure 
// File:                        union_find.h
//----------------------------------------------------------------------
// Copyright (c) 2018 Matt Piekenbrock. All Rights Reserved.
//
// Class definition based off of data-structure described here:  
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure

#include <Rcpp.h>
using namespace Rcpp;

struct UnionFind {
  const int size; 
  Rcpp::IntegerVector parent;
  Rcpp::IntegerVector rank;
  
  UnionFind(const int _size);
  ~UnionFind();
  SEXP as_XPtr();
  void Union(const int x, const int y); 
  void UnionAll(const IntegerVector idx);
  const int Find(const int x); 
  IntegerVector FindAll(const IntegerVector idx);
  IntegerVector getCC();
  void printCC();
}; // class UnionFind