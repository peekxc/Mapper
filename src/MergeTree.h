#include <Rcpp.h>
using namespace Rcpp;

#include "SimplexTree.h"

using std::size_t;
using std::map;
using std::pair;
using std::vector;

struct super_node {
  size_t id; 
  double height;
};
struct super_edge {
  vector< size_t > ids; 
  vector< double > heights; 
};

// Merge tree declaration
// class MergeTree {
// public: 
//   SimplexTree relations; // the relations composing the merge tree 
//   vector< super_node > super_nodes; // super node heights 
//   map< pair< size_t, size_t >, super_edge > super_edges; // super edges
//   
//   // Ctor + Dtor
//   MergeTree();
//   ~MergeTree();
// 
//   // Constructs the merge tree.
//   SEXP as_XPtr();
//   List construct(SEXP smesh, const NumericVector& h);
//   List export_to_list();
// };
