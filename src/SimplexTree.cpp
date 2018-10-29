// Simple but limited implementation of the Simplex tree data structure using Rcpp + STL
// Original Reference: Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: 
// An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.

#include "SimplexTree.h"

// Constructs the simplex tree
SimplexTree::SimplexTree() {  
  root = node_ptr(new node(0, nullptr));
  level_map = std::unordered_map< std::string, std::vector<node_ptr> >();
  n_simplexes = std::vector<uint>(); 
}

SimplexTree::~SimplexTree() { 
  std::vector<uint>().swap(n_simplexes);
  // delete_tree(root);
}

SEXP SimplexTree::as_XPtr(){
  Rcpp::XPtr< SimplexTree> p(this, false);
  return(p);
}

// Search the level map (cousins) to quickly get the adjacency relations
IntegerVector SimplexTree::adjacent_vertices(const int v){

  // Resulting vector to return
  IntegerVector res = IntegerVector(); 
  
  // First extract the vertices which labels > v
  std::string key = std::to_string(v) + "-1";
  std::unordered_map< std::string, std::vector<node_ptr> >::iterator it = level_map.find(key);
  if (it != level_map.end()){
    // Get the cousins parents 
    std::vector<node_ptr> cousins = (*it).second;
    std::for_each(cousins.begin(), cousins.end(), [&res](const node_ptr& c_node ){
      res.push_back((*c_node).parent->label);
    });
  }
  
  // Then get the vertices with labels < v
  std::map< uint, node_ptr > top_vertices = root->children;
  if (top_vertices.find(v) != top_vertices.end()){
    node_ptr vertex = top_vertices.at(v);
    std::map< uint, node_ptr > v_children = vertex->children;
    for (auto& kv: v_children){
      res.push_back(kv.first);
    }
  }
  
  // Return 
  return(unique(res));
}

void SimplexTree::record_new_simplexes(const uint k, const uint n){
  if (n_simplexes.size() < k+1){
    n_simplexes.resize(k+1);
  }
  n_simplexes.at(k) += n;
}

// Finds the 0-based index of the vertex in the top nodes, or -1 otherwise.
int SimplexTree::find_vertex(const int v_id){
  std::map< uint, node_ptr > top_vertices = root->children;
  auto it = top_vertices.find(v_id); 
  if (it != top_vertices.end()){ 
    return(std::distance(top_vertices.begin(), it));
  }
  return(-1);
}
  
// Returns a 0-based integer vector of the first ids available to allocate new vertices
IntegerVector SimplexTree::vertex_available(uint n_vertices){
  std::map< uint, node_ptr > top_vertices = root->children;
  uint max = top_vertices.size() + n_vertices;
  IntegerVector new_idx = IntegerVector();
  for (uint cc = 0; cc < max; ++cc){
    if (top_vertices.find(cc) == top_vertices.end()){
      new_idx.push_back(cc);
      if (new_idx.size() == n_vertices){ break; }
    }
  }
  return(new_idx);
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

// Removes a vertex from the simplex tree, including all edges connected to it.
void SimplexTree::remove_vertices(IntegerVector vertex_ids){
  if (vertex_ids.size() == 0){ return; }
  std::map< uint, std::shared_ptr<node> > top_vertices = root->children;
  IntegerVector edge_to_remove = IntegerVector::create(0, 0);
  for (IntegerVector::const_iterator vid = vertex_ids.begin(); vid != vertex_ids.end(); ++vid){
    
    // First remove any of its cofaces
    remove_vertex_cofaces(*vid);
    
    // Then remove the edges with labels > the query vertex
    for (auto& top_vertex: top_vertices) {
      if (top_vertex.first >= *vid){ continue; }
      edge_to_remove[0] = top_vertex.first;
      edge_to_remove[1] = *vid;
      remove_edge(edge_to_remove);
    }
    
    // Remove the vertex itself and its edges with label < the query vertex  
    if (root->children.find(*vid) != root->children.end()){
      node_ptr c_node = root->children.at(*vid);
      int n_children = c_node->children.size(); 
      c_node->children.clear();
      record_new_simplexes(1, -n_children);
      remove_child(root, *vid, 0);
    }
    
    // Remove any reference in the level map
    std::string key = std::to_string(*vid) + "-1";
    level_map.erase(key);
  }
}

// Removes all cofaces containing vertex v
void SimplexTree::remove_vertex_cofaces(const int v){
  using simplices = std::map< uint, std::shared_ptr<node> >;
  simplices top_nodes = root->children;
  simplices::iterator it = top_nodes.find(v);
  if (it != top_nodes.end()){
    
    // First, find the query vertex's connected edges with label > than the vertex 
    simplices v_children = (*it).second->children;
    for (auto& child: v_children){
      
      // Then get those vertices cousins
      std::string key = std::to_string(child.first)+"-1";
      if (level_map.find(key) != level_map.end()){
        std::vector< node_ptr >& adj_vertices = level_map[key];
        
        // If the query vertex is the parent of any such vertices, remove it 
        adj_vertices.erase(
          std::remove_if(adj_vertices.begin(), adj_vertices.end(), [&v](node_ptr vi){
            return(vi->parent->label == v);
          }), 
          adj_vertices.end()
        );
      }
    }
  }
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
  

// Remove a child node directly to the parent if it exists
void SimplexTree::remove_child(node_ptr c_parent, uint child_label, uint depth){
  //Rcout << c_parent << std::endl;
  if (c_parent == nullptr){ return; }
  std::map<uint, node_ptr>& parents_children = c_parent->children;
  std::map<uint, node_ptr>::iterator key_lb = parents_children.find(child_label);
  if (key_lb == parents_children.end()) { return; }
  else {
    // TODO: Remove the for_each loop and call .clear() instead, crashing because of iterator invalidation
    // node_ptr& c_node = key_lb->second;
    // if (c_node->children.size() > 0){
    //   std::for_each(c_node->children.begin(), c_node->children.end(), [&](std::map<uint, node_ptr>::value_type& pair)   {
    //     remove_child(c_node, pair.first, depth+1);
    //   });
    // }
    // uint key = key_lb->first; 
    // Finally, remove the intended node
    parents_children.erase(key_lb);
  }
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
        level_map[key].push_back(new_node);
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
  if (simplex.size() == 1){ 
    return(root->children.find(simplex.at(0)) != root->children.end()); 
  } 
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

// Prints all the cofaces at a given depth
void SimplexTree::print_cofaces(int depth){
  // root->children(); 
  for (auto& kv: level_map){
    std::string key = kv.first;
    std::vector<node_ptr> val = kv.second;
    Rcout << key << ": ";
    std::for_each(val.begin(), val.end(), [](const node_ptr v){
      Rcout << v->label << "-^" << v->parent->label << ", ";
    });
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

// Experimental k-expansion algorithm.
// Performs an expansion of order k, thus reconstructing a k-skeleton from the 1-skeleton alone.
void SimplexTree::expansion(const uint k){
  std::for_each(root->children.begin(), root->children.end(), [&](const std::pair<uint, node_ptr>& c_node){
    //uint simplex[k];
    //expand(c_node.second->children, k, 0, simplex);
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
      // Rcout << "Adding children " << sc_int1 << " to node " << v_it->first << " (son of " << v_it->second->parent->label << ")" << std::endl; 
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
      const int i = find_vertex(v.first), j = find_vertex(child.first);
      res.at(i, j) = 1; 
      res.at(j, i) = 1; 
    });
  });
  return(res);
}
  
// Exports the 1-skeleton as an adjacency matrix 
List SimplexTree::as_adjacency_list(){
  const size_t n = root->children.size();
  std::unordered_map< std::string, std::vector<uint> > res(n); // output
  std::for_each(root->children.begin(), root->children.end(), [&](std::pair<uint, node_ptr> v){
    const std::map<uint, node_ptr> vc = v.second->children;
    std::vector<uint> adjacencies = std::vector<uint>();
    std::for_each(vc.begin(), vc.end(), [&](std::pair<uint, node_ptr> child){
      res[std::to_string(v.first)].push_back(child.first);
      res[std::to_string(child.first)].push_back(v.first);
    });
  });
  return(wrap(res));
}
  
// Exports the 1-skeleton as an edgelist 
IntegerMatrix SimplexTree::as_edge_list(){
  // const size_t n = root->children.size();
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




// Exposed Rcpp Module 
RCPP_MODULE(simplex_tree_module) {
  Rcpp::class_<SimplexTree>("SimplexTree")
  .constructor()
  .field_readonly("n_simplexes", &SimplexTree::n_simplexes)
  .method( "as_XPtr", &SimplexTree::as_XPtr)
  .method( "add_vertices", &SimplexTree::add_vertices)
  .method( "remove_vertices", &SimplexTree::remove_vertices)
  .method( "vertex_available", &SimplexTree::vertex_available)
  .method( "adjacent_vertices", &SimplexTree::adjacent_vertices)
  .method( "insert_simplex", &SimplexTree::insert_simplex)
  .method( "find_simplex", &SimplexTree::find_simplex)
  .method( "remove_edge", &SimplexTree::remove_edge)
  .method( "print_tree", &SimplexTree::print_tree )
  .method( "print_cofaces", &SimplexTree::print_cofaces )
  // .method( "expansion", &SimplexTree::expansion )
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
# stree_ptr <- stree$as_XPtr()
stree$insert_simplex(as.integer(c(1, 2)))
stree$insert_simplex(as.integer(c(1, 3)))
stree$insert_simplex(as.integer(c(2, 3)))
stree$insert_simplex(as.integer(c(2, 4)))
stree$insert_simplex(as.integer(c(2, 5)))
stree$insert_simplex(as.integer(c(3, 4)))
stree$insert_simplex(as.integer(c(3, 5)))
stree$insert_simplex(as.integer(c(4, 5)))

stree$print_tree()

stree$remove_vertices(4)

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
