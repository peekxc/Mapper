// Implementation of the Simplex tree data structure using Rcpp + STL
// Original Reference: Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: 
// An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.

#include "SimplexTree.h"

// void updateCousins(node_ptr c_node, std::vector< node_ptr >& eq_cousins){
//   if (c_node == nullptr || c_node->parent == nullptr){ return; }
//   node_ptr node_parent = c_node->parent;
//   std::map< uint, node_ptr > parent_siblings = node_parent->parent->children; // the grandparents children are
//   std::map< uint, node_ptr >::iterator simplex_it;
//   eq_cousins.clear(); // note that node pointers are *not* deleted (nor should they be); they are just removed from the container
//   for (std::map< uint, node_ptr >::iterator parent_it = parent_siblings.begin(); parent_it != parent_siblings.end(); ++parent_it){
//     // If the current parent doesn't equal the nodes parent, it's a sibling of the node's parent
//     node_ptr next_parent = parent_it->second;
//     if (next_parent != node_parent){
//       std::map< uint, node_ptr > potential_cousins = next_parent->children;
//       // The current node's parents siblings children are the nodes cousins if they share the same label
//       if ((simplex_it = potential_cousins.find(c_node->label)) != potential_cousins.end()){
//         eq_cousins.push_back(simplex_it->second);
//       }
//     }
//   }
// }

// Constructs the simplex tree
SimplexTree::SimplexTree() {  
  root = node_ptr(new node(0, nullptr));
  // level_map = std::unordered_map< std::string, std::vector<node_ptr> >();
  n_simplexes = std::vector<uint>(); 
}

SimplexTree::~SimplexTree() { 
  std::vector<uint>().swap(n_simplexes);
  // delete_tree(root);
}

void SimplexTree::remove_edge(IntegerVector edge){
  if (edge.size() != 2){ stop("Invalid query. 'remove_edge' takes as input a vector of length 2."); }
  if (edge[0] > edge[1]){ int tmp = edge[0]; edge[0] = edge[1]; edge[1] = tmp; }
  bool edge_exists = find_simplex(edge);
  if (!edge_exists){ return; }
  else {
    std::shared_ptr<node> v_ptr = root->children.at(edge[0]);
    remove_child(v_ptr, edge[1], 1);
  }
}

SEXP SimplexTree::as_XPtr(){
  Rcpp::XPtr< SimplexTree> p(this, false);
  return(p);
}

void SimplexTree::record_new_simplexes(const uint k, const uint n){
  if (n_simplexes.size() < k+1){
    n_simplexes.resize(k+1);
  }
  n_simplexes.at(k) += n;
}
  
// Emplace v_i new 0-simplices. Ids are automatically assigned. 
void SimplexTree::add_vertices(const uint v_i){
  const size_t n = root->children.size();
  for (size_t i = 0; i < v_i; ++i){ 
    size_t v_id = n+i+1;
    root->children.emplace(v_id, node_ptr(new node(v_id, root))); 
  }
  record_new_simplexes(0, v_i);
}
  

// Remove a child node directly to the parent if it exists
void SimplexTree::remove_child(node_ptr c_parent, uint child_label, uint depth){
  std::map<uint, node_ptr>& parents_children = c_parent->children;
  std::map<uint, node_ptr>::iterator key_lb = parents_children.find(child_label); 
  if (key_lb == parents_children.end()) { return; } 
  else { parents_children.erase(key_lb); }
  record_new_simplexes(depth, -1);
}
  
// Attachs a child node directly to the parent if it doesn't exist
void SimplexTree::add_child(node_ptr c_parent, uint child_label, uint depth){
  std::map<uint, node_ptr>& parents_children = c_parent->children;
  std::map<uint, node_ptr>::iterator key_lb = parents_children.find(child_label); 
  if (key_lb == parents_children.end()) { 
    node_ptr new_child = node_ptr(new node(child_label, c_parent));
    parents_children.insert(key_lb, std::map<uint, node_ptr>::value_type(child_label, new_child));
  }
  record_new_simplexes(depth, 1);
}

// Creates a new set of child nodes for the given parent. createChild check for redundancy with insert. 
// Returns all of the children of the node upon completion. 
void SimplexTree::add_children(node_ptr c_parent, const std::vector<uint>& new_children, uint depth){
  std::for_each(new_children.begin(), new_children.end(), [&](const uint& id){
    add_child(c_parent, id, depth);
  });
}
  
// Rcpp wrapper needs STL for API access
inline void SimplexTree::insert_simplex(std::vector<uint> labels){
  std::sort(labels.begin(), labels.end()); // Demand that labels are sorted on insertion! 
  insert((uint*) &labels[0], 0, labels.size(), root, 0); // start the recursion from the root
}


void SimplexTree::insert(uint* labels, const size_t i, const size_t n_keys, node_ptr c_node, const uint depth){
  if (i >= n_keys || labels == nullptr || c_node == nullptr){ return; } // base case + safety checks
  // Create a set of (i)-simplexes as children of the current node, if they don't already exist
  for (int j = i; j < n_keys; ++j){
    size_t j_exists = c_node->children.count(labels[j]);
    if (!bool(j_exists)){
      // Rcout << "Creating new node " << labels[j] << ", parent: " << c_node->label << std::endl;
      node_ptr new_node = node_ptr(new node(labels[j], c_node));
      add_child(c_node, labels[j], depth);
      if (depth > 0){ // keep track of nodes which share ids at the same depth
        std::string key = std::to_string(labels[j]) + "-" + std::to_string(depth);
        // level_map[key].push_back(new_node);
      }
    }
  }
  // Recurse on the subtrees of the current node
  for (int j = i; j < n_keys; ++j){
    // Rcout << "Recursing with child node: " << labels[j] <<  " of parent: " << c_node->label << " and grandparent: " << (c_node->parent == nullptr ? 0 : c_node->parent->label)  << std::endl;
    node_ptr child_node = c_node->children.at(labels[j]);
    insert(labels, j + 1, n_keys, child_node, depth + 1);
  }
}

// Rcpp wrapper to the find function
bool SimplexTree::find_simplex(const IntegerVector& simplex){
  if (simplex.size() == 0){ return false; }
  if (simplex.size() == 1){ return(simplex.at(0) >= 1 && simplex.at(0) <= root->children.size()); } 
  else {
    std::vector<int> simplex_query = as< std::vector<int> >(simplex);
    node_ptr res = find(simplex_query);
    return(bool(res));
  }
}
  
// Overload to get the top node 
node_ptr SimplexTree::find (uint label){
  std::map<uint, node_ptr>::iterator it = root->children.find(label);
  if (it != root->children.end()){
    return(it->second); 
  } else { return(nullptr); }
}
  
// Given an integer label, searches the tree to see if the simplex exists. If so, the simplex
// (node) is returned, else a nullptr is returned.
node_ptr SimplexTree::find (std::vector<int> simplex){
  node_ptr c_node = root;
  std::map<uint,node_ptr>::iterator node_it;
  for (size_t i = 0; i < simplex.size(); ++i){
    if ((node_it = c_node->children.find(simplex[i])) != c_node->children.end()){
      c_node = node_it->second;
    } else { return nullptr; }
  }
  return(c_node);
}
  

// utility to get the maximumheight /longest path any a given node.
uint SimplexTree::get_dfs_height(node_ptr cnode, uint c_height){
  std::map<uint, node_ptr>::iterator it; 
  uint max_height = 0; 
  for (it = cnode->children.begin(); it != cnode->children.end(); ++it){ 
    int tmp = get_dfs_height(it->second, c_height + 1); 
    if (tmp > max_height){ max_height = tmp; }
  }
  return((uint) std::max(c_height, max_height));
}
  
void SimplexTree::print_level(node_ptr cnode, uint level){
  if (cnode == nullptr || cnode == NULL) return;
  if (level == 0) Rprintf(" %d", cnode->label);
  else if (level > 0)
  {
    std::map<uint, node_ptr>::iterator it;
    for (it = cnode->children.begin(); it != cnode->children.end(); ++it){ 
      print_level(it->second, level-1);
    }
  }
}
  
// Basic breadth-first printing. Each level is prefixed with '.' <level> number of times, followed by the 
// the ids of the nodes at that breadth-level enclosed within parenthesis, i.e. ..( 2 3 4 ) 
void SimplexTree::print_tree(){
  std::map<uint, node_ptr>::iterator it;
  for (it = root->children.begin(); it != root->children.end(); ++it){ 
    uint h = get_dfs_height(it->second, 0); 
    Rcout << it->first << " (h = " << h << "): ";
    for (int i = 1; i <= h; ++i){ 
      for (int j = 1; j <= i; ++j){  Rcout << "."; }
      Rcout << "("; print_level(it->second, i); Rcout << " )";
    }
    Rcout << std::endl;
  }
}
  
std::vector<uint> SimplexTree::getLabels(const std::map<uint, node_ptr>& level, const uint offset){
  // return from iterator begin to end labels as vector
  std::vector<uint> labels;
  labels.reserve(level.size() - offset);
  std::map<uint, node_ptr>::const_iterator it = level.begin(); 
  std::advance(it, offset);
  std::transform (it, level.end(), back_inserter(labels), 
    [&](std::pair<uint, node_ptr> const& pair) { return pair.first; });
  return(labels);
}
  
uint SimplexTree::intersection_size(std::vector<uint> v1, std::vector<uint> v2){
  std::unordered_set<uint> s(v1.begin(), v1.end());
  uint res = count_if(v2.begin(), v2.end(), [&](uint k) {return s.find(k) != s.end(); });
  return(res);
}

// v1 and v2 should already by sorted!!
std::vector<uint> SimplexTree::intersection(std::vector<uint> v1, std::vector<uint> v2){
  std::vector<uint> v3;
  std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));
  return(v3);
}

// Performs an expansion of order k, thus reconstructing a k-skeleton from the 1-skeleton alone.
void SimplexTree::expansion(const uint k){
  std::for_each(root->children.begin(), root->children.end(), [&](const std::pair<uint, node_ptr>& c_node){
    uint simplex[k];
    expand(c_node.second->children, k, 0, simplex);
  });
}
  
// Expand operation compares a given 'head' nodes children to its siblings. 
// If they have a non-empty intersection, then the intersection is added as a child to the head node. 
void SimplexTree::expand(std::map<uint, node_ptr>& v, const uint k, uint depth, uint* simplex){
  
  if (v.size() <= 1){ return; } // Current level only has one node; intersection will be empty
  
  // For each child node
  std::map<uint, node_ptr>::const_iterator v_it = v.begin(); 
  for (v_it = v.begin(); v_it != v.end(); ++v_it){
    
    // Get the 'head' nodes 
    node_ptr rel_head = v_it->second; // *relative* head
    node_ptr root_head = find(v_it->first); // *root* head
    
    // If the root/0-simplex of the head doesn't have children, we're done
    if (root_head->children.size() == 0){ return; }
    
    // Get the (1-offset) siblings of the relative head, and the labels of the children of root head
    std::vector<uint> siblings = getLabels(v, 1);
    std::vector<uint> children = getLabels(root_head->children, 0);
    std::vector<uint> sc_int = intersection(children, siblings);
    
    IntegerVector sib1 = wrap(siblings);
    IntegerVector children1 = wrap(children);
    IntegerVector sc_int1 = wrap(sc_int);
    Rcout << sib1 << std::endl; 
    Rcout << sc_int1 << std::endl; 
    if (sc_int.size() > 0){
      Rcout << "Adding children " << sc_int1 << " to node " << v_it->first << " (son of " << v_it->second->parent->label << ")" << std::endl; 
      add_children(rel_head, sc_int, depth + 2);
      expand(rel_head->children, k, depth + 1, simplex);
    }
  }
}
  
// Exports the 1-skeleton as an adjacency matrix 
IntegerMatrix SimplexTree::as_adjacency_matrix(){
  const size_t n = root->children.size();
  IntegerMatrix res = IntegerMatrix(n, n);
  std::for_each(root->children.begin(), root->children.end(), [&](std::pair<uint, node_ptr> v){
    const std::map<uint, node_ptr> vc = v.second->children;
    std::for_each(vc.begin(), vc.end(), [&](std::pair<uint, node_ptr> child){
      res.at(v.first-1, child.first-1) = 1; 
      res.at(child.first-1, v.first-1) = 1; 
    });
  });
  return(res);
}
  
// Exports the 1-skeleton as an adjacency matrix 
List SimplexTree::as_adjacency_list(){
  const size_t n = root->children.size();
  std::vector< std::vector<uint> > res(n); // output
  std::for_each(root->children.begin(), root->children.end(), [&](std::pair<uint, node_ptr> v){
    const std::map<uint, node_ptr> vc = v.second->children;
    std::vector<uint> adjacencies = std::vector<uint>();
    std::for_each(vc.begin(), vc.end(), [&](std::pair<uint, node_ptr> child){
      res.at(v.first-1).push_back(child.first);
      res.at(child.first-1).push_back(v.first);
    });
  });
  return(wrap(res));
}
  
// Exports the 1-skeleton as an edgelist 
IntegerMatrix SimplexTree::as_edge_list(){
  const size_t n = root->children.size();
  uint n_edges = 0; 
  std::for_each(root->children.begin(), root->children.end(), [&](std::pair<uint, node_ptr> v){
    n_edges += v.second->children.size();
  });
  
  int i = 0; 
  IntegerMatrix res = no_init_matrix(n_edges, 2);
  std::for_each(root->children.begin(), root->children.end(), [&](std::pair<uint, node_ptr> v){
    const std::map<uint, node_ptr> vc = v.second->children;
    std::for_each(vc.begin(), vc.end(), [&](std::pair<uint, node_ptr> child){
      res.row(i++) = IntegerVector::create(v.first, child.first); 
    });
  });
  return(res);
}
  
  // IntegerMatrix export_k_simplexes(const int k, Integer){
  //   
  // }
  
  // Exports the k-skeleton as a list
  // List as_list(){
  //   
  // }
  // 
  // void SimplexTree::as_list_helper(List& res, node_ptr c_node, uint depth, uint* simplex, const size_t n_keys){
  //   
  //   // Base case
  //   if (c_node->children.size() == 0){
  //     IntegerMatrix& simplices = res.at(depth);
  //     
  //   }
  // }


// List construct_simplex_tree(const IntegerMatrix& el, const int n_nodes) {

  // uint simplex_edge[2] = { (uint) 1, (uint) 2 };
  // mapper_tree->insert_simplex(simplex_edge, 2);
  // 
  // simplex_edge[0] = 3;
  // simplex_edge[1] = 4;
  // mapper_tree->insert_simplex(simplex_edge, 2);
  // 
  // uint simplex_tri[3] = { (uint) 1, (uint) 2 , (uint) 3 };
  // mapper_tree->insert_simplex(simplex_tri, 3);
  // 
  // simplex_tri[0] = 7;
  // simplex_tri[1] = 8;
  // simplex_tri[2] = 9;
  // mapper_tree->insert_simplex(simplex_tri, 3);
  
  // (uint* labels, const size_t i, const size_t n_keys, node_ptr c_node){


  // mapper_tree->insert(simplex_edge, 2);
  //
  // uint simplex_edge2[3] = { (uint) 1, (uint) 2, (uint) 3 };
  // mapper_tree->insert(simplex_edge2, 3);
  //
  // uint simplex_edge3[3] = { (uint) 2, (uint) 3, (uint) 4 };
  // mapper_tree->insert(simplex_edge3, 3);

//   mapper_tree->print_tree();
// // mapper_tree->toList(mapper_tree->root)
//   // mapper_tree->debug();
//   return(List::create());
// }

// Exposed Rcpp Module 
RCPP_MODULE(simplex_tree_module) {
  Rcpp::class_<SimplexTree>("SimplexTree")
  .constructor()
  .field_readonly("n_simplexes", &SimplexTree::n_simplexes)
  .method( "as_XPtr", &SimplexTree::as_XPtr)
  .method( "add_vertices", &SimplexTree::add_vertices)
  .method( "insert_simplex", &SimplexTree::insert_simplex)
  .method( "find_simplex", &SimplexTree::find_simplex)
  .method( "remove_edge", &SimplexTree::remove_edge)
  .method( "print_tree", &SimplexTree::print_tree )
  .method( "expansion", &SimplexTree::expansion )
  .method( "as_adjacency_matrix", &SimplexTree::as_adjacency_matrix )
  .method( "as_adjacency_list", &SimplexTree::as_adjacency_list)
  .method( "as_edge_list", &SimplexTree::as_edge_list)
  ;
}

/*** R
library("Mapper")
n_vertices <- 5L
stree <- Mapper::simplex_tree()
stree$add_vertices(n_vertices)
rm(stree)
gc()

stree <- Mapper::simplex_tree()
stree_ptr <- stree$as_XPtr()
stree$insert_simplex(as.integer(c(1, 2)))
stree$insert_simplex(as.integer(c(1, 3)))
stree$insert_simplex(as.integer(c(2, 3)))
stree$insert_simplex(as.integer(c(2, 4)))
stree$insert_simplex(as.integer(c(2, 5)))
stree$insert_simplex(as.integer(c(3, 4)))
stree$insert_simplex(as.integer(c(3, 5)))
stree$insert_simplex(as.integer(c(4, 5)))

stree$print_tree()

stree$expansion(2)

stree$insert_simplex(as.integer(c(1, 2, 3)))
stree$insert_simplex(as.integer(c(2, 3, 4)))
stree$insert_simplex(as.integer(c(2, 4, 5)))
stree$insert_simplex(as.integer(c(3, 4, 5)))
stree$print_tree()

stree$insert_simplex(as.integer(c(1, 2, 3, 4, 5)))
stree$print_tree()

# ## Generate noisy points around the perimeter of a circle 
# n <- 150
# t <- 2*pi*runif(n)
# r <- runif(n, min = 2, max = 2.1)
# circ <- cbind(r*cos(t), r*sin(t))
# 
# ## Define filter values equal to the distance from each point to the left-most point in the circle 
# left_pt <- circ[which.min(circ[, 1]),]
# f_x <- sqrt(colSums(apply(circ, 1, function(pt) (pt - left_pt)^2)))
# 
# # g <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
# m_ref <- Mapper::mapper(X = circ, filter_values = f_x, number_intervals = 5, overlap = 0.20, return_reference = TRUE)
# ls_pairs <- t(combn(length(m_ref$cover$level_sets), 2))
# node_lsfi <- sapply(m_ref$G$nodes, function(node) attr(node, "level_set"))
# ls_node_map <- lapply(seq(length(m_ref$cover$level_sets)), function(lvl_set_idx) {
#   node_indices <- which(node_lsfi == lvl_set_idx)
#   if (length(node_indices) == 0){ return(NULL) } else { return(node_indices) }
# })
# el <- Mapper:::edgeList_int(ls_pairs = ls_pairs, nodes = m_ref$G$nodes, ls_node_map = ls_node_map)
# sorted_el <- t(apply(el, 1, sort))
# wut <- Mapper:::construct_simplex_tree(sorted_el, length(m_ref$G$nodes))
*/
