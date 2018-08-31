// Implementation of the Simplex tree data structure using Rcpp + STL
// Original Reference: Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: 
// An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

typedef unsigned int uint;
#include <map>
#include <unordered_map>
#include <queue>
#include <vector>

// Node structure. Contains the following fields:
//  label := unsigned integer representing the id of simplex it represents
//  parent := node pointer to its parent in the trie
//  children := connected simplexes whose labels > the current simplex's label
struct node {
  uint label;
  node* parent;
  std::map<uint, node*> children;
  node(uint id, node* c_parent) : label(id), parent(c_parent){ children = std::map<uint, node*>(); }
  void createChild(uint id){ children.insert(std::make_pair(id, new node(id, this))); }
};


// Simplex tree data structure. 
struct SimplexTree {
  node* root; // empty face; initialized to id = 0, parent = nullptr
  std::unordered_map< std::string, std::vector<node*> > level_map; // maps strings of the form "<id>-<depth>" to a vector of node pointers
  
  // Constructor
  SimplexTree();
  
  // API functions + associated internal recursive overloads
  void add_vertices(const int v_i);
  void add_children(node* c_parent, const std::vector<int>& new_children);
  void attach_child(node* c_parent, node* child);

  void insert_simplex(std::vector<uint> labels);
  void insert_simplex_int(uint* labels, const size_t n_keys);
  void insert(uint* labels, const size_t i, const size_t n_keys, node* c_node, const int depth);

  bool find_simplex(const IntegerVector& simplex);
  node* find (int label);
  node* find (std::vector<int> simplex);
  
  
  // Utility 
  int get_dfs_height(node* cnode, int c_height);
  void print_level(node* cnode, int level);
  void print_tree();
  std::vector<int> getLabels(const std::map<uint, node*>& level, const int offset = 0);
  int intersection_size(std::vector<int> v1, std::vector<int> v2);
  std::vector<int> intersection(std::vector<int> v1, std::vector<int> v2);
  
  // Expansion function
  void expansion(const int k);
  void expand(std::map<uint, node*>& v, const int k, int depth, uint* simplex);
  
  // Export utilities
  IntegerMatrix as_adjacency_matrix(); // Exports the 1-skeleton as an adjacency matrix 
  List as_adjacency_list(); // Exports the 1-skeleton as an adjacency matrix 
  IntegerMatrix as_edge_list(); // Exports the 1-skeleton as an edgelist 
};
