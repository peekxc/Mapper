// Implementation of the Simplex tree data structure using Rcpp + STL
// Original Reference: Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes."
//                     Algorithmica 70.3 (2014): 406-427.

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

uint* extractKeys(uint key){
  const uint n_keys = std::ceill(std::log10l(key+1));
  uint* key_array = new uint[n_keys];
  for (int i = n_keys - 1; i >= 0; i--) { key_array[i] = key % 10; key/= 10; }
  return(key_array);
}

// void updateCousins(node* c_node, std::vector< node* >& eq_cousins){
//   if (c_node == nullptr || c_node->parent == nullptr){ return; }
//   node* node_parent = c_node->parent;
//   std::map< uint, node* > parent_siblings = node_parent->parent->children; // the grandparents children are
//   std::map< uint, node* >::iterator simplex_it;
//   eq_cousins.clear(); // note that node pointers are *not* deleted (nor should they be); they are just removed from the container
//   for (std::map< uint, node* >::iterator parent_it = parent_siblings.begin(); parent_it != parent_siblings.end(); ++parent_it){
//     // If the current parent doesn't equal the nodes parent, it's a sibling of the node's parent
//     node* next_parent = parent_it->second;
//     if (next_parent != node_parent){
//       std::map< uint, node* > potential_cousins = next_parent->children;
//       // The current node's parents siblings children are the nodes cousins if they share the same label
//       if ((simplex_it = potential_cousins.find(c_node->label)) != potential_cousins.end()){
//         eq_cousins.push_back(simplex_it->second);
//       }
//     }
//   }
// }

// a node B is a cousin of node A if:
// 1) node B and node A have the same label, and
// 2) node A and B's parents are siblings
struct simplex_tree {
  node* root; // empty face; initialized to id = 0, parent = nullptr
  std::unordered_map< std::string, std::vector<node*> > level_map; // maps strings of the form "<id>-<depth>" to a vector of node pointers

  simplex_tree(const int n_nodes) {  // Constructs the simplex tree
    root = new node(0, nullptr);
    for (size_t i = 0; i < n_nodes; ++i){ root->children.emplace(i+1, new node(i+1, root)); }
    level_map = std::unordered_map< std::string, std::vector<node*> >();
  }


  void attachChild(node* c_parent, node* child){
    c_parent->children[child->label] = child; // .insert(std::pair<uint, node*>(child->label, child));
  }

  void insert_simplex(uint* labels, const size_t n_keys){
    insert(labels, 0, n_keys, root, 0); // start the recursion from the root
  }

  void insert(uint* labels, const size_t i, const size_t n_keys, node* c_node, const int depth){
    if (i >= n_keys || labels == nullptr || c_node == nullptr){ return; } // base case + safety checks
    // Create a set of (i)-simplexes as children of the current node, if they don't already exist
    for (int j = i; j < n_keys; ++j){
      size_t j_exists = c_node->children.count(labels[j]);
      if (!bool(j_exists)){
        // Rcout << "Creating new node " << labels[j] << ", parent: " << c_node->label << std::endl;
        node* new_node = new node(labels[j], c_node);
        attachChild(c_node, new_node);
        if (depth > 0){ // keep track of nodes which share ids at the same depth
          std::string key = std::to_string(labels[j]) + "-" + std::to_string(depth);
          level_map[key].push_back(new_node);
        }
      }
    }
    // Recurse on the subtrees of the current node
    for (int j = i; j < n_keys; ++j){
      // Rcout << "Recursing with child node: " << labels[j] <<  " of parent: " << c_node->label << " and grandparent: " << (c_node->parent == nullptr ? 0 : c_node->parent->label)  << std::endl;
      node* child_node = c_node->children.at(labels[j]);
      insert(labels, j + 1, n_keys, child_node, depth + 1);
    }
  }

  // Given an integer label, searches the tree to see if the simplex exists. If so, the simplex
  // (node) is returned, else a nullptr is returned.
  node* find (uint label){
    node* c_node = root;
    uint* int_keys = extractKeys(label);
    const size_t n_keys = std::ceill(std::log10l(label+1));
    std::map<uint,node*>::iterator node_it;
    for (size_t i = 0; i < n_keys; ++i){
      if ((node_it = c_node->children.find(int_keys[i])) != c_node->children.end()){
        c_node = node_it->second;
      } else { return nullptr; }
    }
    return(c_node);
  }

  // Basic breadth-first printing
  void print_tree(){
    std::queue<node*, std::list<node*> > queue = std::queue<node*, std::list<node*> >();
    std::unordered_map<node*, bool> visited = std::unordered_map<node*, bool>();
    queue.push(root);
    while (true) {
      int nodeCount = queue.size(); // size of the current level to print
      if (nodeCount == 0){ break; }

      // Dequeue all nodes of current level and Enqueue all nodes of next level
      while (nodeCount > 0) {
        node* c_node = queue.front();
        Rcout << c_node->label << " ";
        queue.pop();
        // Enqueue all children
        for (std::map<uint, node*>::iterator it = c_node->children.begin(); it != c_node->children.end(); ++it){
          node* child_node = it->second;
          if (!visited.count(child_node)){ queue.push(child_node); }
        }
        nodeCount--;
      }
      Rcout << std::endl;
    }
  }

  IntegerVector getChildrenLabels(node* c_node){
    IntegerVector labels = IntegerVector();
    for (auto const& child_node: c_node->children){ labels.push_back(child_node.first); }
    return(labels);
  }
  List getTreeLevel(node* c_node){
    if (c_node->children.empty()){
      return R_NilValue;
    } else{
      List res = List::create(getChildrenLabels(c_node));
      // for (auto const& child_node: c_node->children){
      //   res.push_back(getChildrenLabels(child_node.second));
      // }
      return(res);
    }
  }
  void printChildren(node* c_node, bool addresses = false){
    Rcout << "Printing children of: " << c_node->label << " --> ";
    std::map<uint, node*> children = c_node->children;
    for (std::map<uint, node*>::iterator it = children.begin(); it != children.end(); ++it){
      Rcout << it->first;
      if (addresses){  Rcout << "(" << it->second << ")"; }
      Rcout << ", ";
    }
    Rcout << std::endl;
  }
  void debug(){

    // Prints the level mapping, and the id + parents of each node
    for(auto kv : level_map) {
      Rcout << kv.first << ": ";
      for (std::vector<node*>::iterator node_it = kv.second.begin(); node_it != kv.second.end(); ++node_it){
        Rcout << "(id: " << (*node_it)->label << ", parent: " << (*node_it)->parent->label << "), ";
      }
      Rcout << std::endl;
    }
  }

  // Performs an expansion of order k, thus reconstructing a k-skeleton from the 1-skeleton alone.
  // void expansion(int k){
  //   root->children.begin()
  //   for ()
  // }
  //
  // bool intersect(int dim){
  //   if (dim == 0){ return false; }
  //   if (dim == 1){ return true; }
  //   else ()
  // }

  // Convert simplex tree to a list
  // List childrenToList(node* c_node){
  //   List res = getTreeLevel(c_node);
  //   if (res.size() == 0){
  //     return(res);
  //   } else {
  //     List new_depth = List::create();
  //     for (auto const& child_node: c_node->children){
  //       List child_res = toList(child_node.second);
  //       if (child_res.size() > 0){ new_depth.push_back(child_res); }
  //     }
  //     return (new_depth);
  //   }
  //   return(res);
  // }
};


// [[Rcpp::export]]
List construct_simplex_tree(const IntegerMatrix& el, const int n_nodes) {

  simplex_tree* mapper_tree = new simplex_tree(0);

  // const size_t n = el.nrow();
  // for (size_t i = 0; i < n; ++i){
  //   Rcpp::IntegerMatrix::ConstRow edge = el.row(i);
  //   uint simplex_edge[2] = { (uint) edge[0], (uint) edge[1] };
  //   mapper_tree->insert_simplex(simplex_edge, 2);
  // }

  uint simplex_edge[2] = { (uint) 1, (uint) 2 };
  mapper_tree->insert_simplex(simplex_edge, 2);

  simplex_edge[0] = 3;
  simplex_edge[1] = 4;
  mapper_tree->insert_simplex(simplex_edge, 2);

  uint simplex_tri[3] = { (uint) 1, (uint) 2 , (uint) 3 };
  mapper_tree->insert_simplex(simplex_tri, 3);

  simplex_tri[0] = 7;
  simplex_tri[1] = 8;
  simplex_tri[2] = 9;
  mapper_tree->insert_simplex(simplex_tri, 3);
  // (uint* labels, const size_t i, const size_t n_keys, node* c_node){


  // mapper_tree->insert(simplex_edge, 2);
  //
  // uint simplex_edge2[3] = { (uint) 1, (uint) 2, (uint) 3 };
  // mapper_tree->insert(simplex_edge2, 3);
  //
  // uint simplex_edge3[3] = { (uint) 2, (uint) 3, (uint) 4 };
  // mapper_tree->insert(simplex_edge3, 3);

  mapper_tree->print_tree();
// mapper_tree->toList(mapper_tree->root)
  // mapper_tree->debug();
  return(List::create());
}


/*** R
  # g <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
  m_ref <- Mapper::mapper(X = circ, filter_values = f_x, number_intervals = 5, overlap = 0.20, return_reference = TRUE)
  ls_pairs <- t(combn(length(m_ref$cover$level_sets), 2))
  node_lsfi <- sapply(m_ref$G$nodes, function(node) attr(node, "level_set"))
  ls_node_map <- lapply(seq(length(m_ref$cover$level_sets)), function(lvl_set_idx) {
    node_indices <- which(node_lsfi == lvl_set_idx)
    if (length(node_indices) == 0){ return(NULL) } else { return(node_indices) }
  })
  el <- Mapper:::edgeList_int(ls_pairs = ls_pairs, nodes = m_ref$G$nodes, ls_node_map = ls_node_map)
  sorted_el <- t(apply(el, 1, sort))
  wut <- Mapper:::construct_simplex_tree(sorted_el, length(m_ref$G$nodes))
*/
