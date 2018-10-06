// Simple but limited implementation of the Simplex tree data structure using Rcpp + STL
// Original Reference: Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: 
// An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#define sptr std::shared_ptr // Shared pointer

typedef unsigned int uint;
#include <map>
#include <unordered_map>
#include <queue>
#include <vector>

// Node structure stored by the simplex tree. Contains the following fields:
//  label := unsigned integer representing the id of simplex it represents
//  parent := node pointer to its parent in the trie
//  children := connected simplexes whose labels > the current simplex's label
struct node {
  uint label;
  std::shared_ptr<node> parent;
  std::map< uint, std::shared_ptr<node> > children;
  node(uint id, std::shared_ptr<node> c_parent) : label(id), parent(c_parent){ 
    children = std::map< uint, std::shared_ptr<node> >(); 
  }
};
typedef std::shared_ptr<node> node_ptr; 


// Simplex tree data structure. 
struct SimplexTree {
  node_ptr root; // empty face; initialized to id = 0, parent = nullptr
  std::unordered_map< std::string, std::vector<node_ptr> > level_map; // maps strings of the form "<id>-<depth>" to a vector of node pointers
  std::vector<uint> n_simplexes;
    
  // Constructor + Destructor 
  SimplexTree();
  ~SimplexTree();
  
  // Utilities
  SEXP as_XPtr();
  void record_new_simplexes(const uint k, const uint n);// record keeping
  
  // User-facing API 
  int find_vertex(const int v_id);
  IntegerVector vertex_available(uint n_vertices);
  IntegerVector adjacent_vertices(const int v);
  void add_vertices(const uint v_i);
  void remove_vertices(IntegerVector vertex_ids);
  void remove_edge(IntegerVector labels);
  void insert_simplex(std::vector<uint> labels);
  bool find_simplex(const IntegerVector& simplex);
  void expansion(const uint k);
  void print_tree();
  
  // Export utilities
  IntegerMatrix as_adjacency_matrix(); // Exports the 1-skeleton as an adjacency matrix 
  List as_adjacency_list(bool one_based = false); // Exports the 1-skeleton as an adjacency matrix 
  IntegerMatrix as_edge_list(bool one_based = false); // Exports the 1-skeleton as an edgelist 
  
  // Recursive helper functions
  void add_child(node_ptr c_parent, uint child_label, uint depth);
  void remove_child(node_ptr c_parent, uint child_label, uint depth);
  void add_children(node_ptr c_parent, const std::vector<uint>& new_children, uint depth);
  void insert(uint* labels, const size_t i, const size_t n_keys, node_ptr c_node, const uint depth);
  node_ptr find (uint label);
  node_ptr find (std::vector<int> simplex);
  // void as_list_helper(List& res, node_ptr c_node, uint depth, uint* simplex, const size_t n_keys);
  
  // Utility 
  uint get_dfs_height(node_ptr cnode, uint c_height);
  void print_level(node_ptr cnode, uint level);
  std::vector<uint> getLabels(const std::map<uint, node_ptr>& level, const uint offset = 0);
  uint intersection_size(std::vector<uint> v1, std::vector<uint> v2);
  std::vector<uint> intersection(std::vector<uint> v1, std::vector<uint> v2);
  void expand(std::map<uint, node_ptr>& v, const uint k, uint depth, uint* simplex);
};
