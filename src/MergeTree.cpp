#include "MergeTree.h"
#include "UnionFind.h"

template <typename T>
std::vector<size_t> order(const std::vector<T> &v) {
  // Initialize original index locations
  std::vector< std::size_t > idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  // Sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
  return idx;
}

template <typename T>
std::pair< std::vector<size_t>, std::vector<T> > ordered_sort(const std::vector<T> &v) {
  std::vector< std::size_t > idx(v.size());  // Initialize original index locations
  iota(idx.begin(), idx.end(), 0);
  std::vector<T> res = std::vector<T>(v.size()); 
  std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
  std::transform(idx.begin(), idx.end(), res.begin(), [&v](const size_t i){ 
    return(v[i]);
  }); // Sort indexes based on comparing values in v
  return std::make_pair(idx, v);
}

// map between point in the simplex tree and it's corresponding height + disjoint set index
struct id_map {
  size_t z_idx; //simplex tree id + 0-based index 
  double height; // point height
  bool operator < (const id_map& o) const {
    return (height < o.height);
  }
};

List export_to_list(const vector<super_node>& sn, const map< pair<size_t, size_t>, super_edge >& se){
  // Collect nodes
  const size_t n = sn.size();
  IntegerVector node_ids = Rcpp::no_init(n);
  NumericVector node_heights = Rcpp::no_init(n);
  for (size_t i = 0; i < n; ++i){
    node_ids[i] = sn[i].id;
    node_heights[i] = sn[i].height;
  }
  
  // Collect edges 
  List edges = List();
  for (auto& edge: se){
    std::ostringstream ss;
    ss << edge.first.first << "," << edge.first.second;
    std::string edge_str = ss.str();
    edges[edge_str] = List::create(_["ids"] = wrap(edge.second.ids), _["heights"] = wrap(edge.second.heights));
  }
  
  // Return  
  List res = List::create(
    _["nodes"] = List::create(_["ids"] = node_ids, _["heights"] = node_heights), 
    _["edges"] = edges
  );
  return(res);
};

// template <typename T1, typename T2>
// struct Bimap{
//   map<T1, T2> x_to_y;
//   map<T2, T1> y_to_x;
//   
// };

// Construct the merge tree attempt #2
// [[Rcpp::export]]
List construct_merge_tree2(SEXP smesh, const NumericVector& h, SEXP merge_tree_res){
  using idx_v = vector< size_t >;
  Rcpp::XPtr<SimplexTree> stree_ptr(smesh); // Collect the simplex tree (mesh of vertices)

  // Sort and order the height function values
  const size_t n = h.size();
  if (stree_ptr->n_simplexes[0] != n){ stop("The number of height values must be equal to the number of 0-simplices in the mesh."); }
  idx_v vids = as< idx_v >(stree_ptr->get_vertices());

  // Generate id -> point maps
  size_t i = 0;
  map< size_t, id_map > pts = map< size_t, id_map >( );
  generate_n(std::inserter(pts, pts.begin()), n, [&h, &i, &vids](){
    id_map pt{ i, h[i] };
    return(std::make_pair(vids[i++], pt));
  });

  // Order vertices by height
  idx_v h_order = order(as< vector<double> >(h));

  // TO get adjacent vertices
  auto adj_v = [&stree_ptr](const size_t v_i){ return(stree_ptr->adjacent_vertices(v_i)); };

  UnionFind ds = UnionFind(n); // CCs
  IntegerVector lowest_vertex = stree_ptr->get_vertices();
  Rcpp::XPtr<SimplexTree> mt_ptr(merge_tree_res);
  idx_v edge = idx_v(2);

  // Rprintf("retrieving v_i from i:%d, idx: %d\n", i, h_order.at(i));
  for (int i = n-1; i >= 0; --i){
    const size_t v_i = vids.at(h_order.at(i));
    const double vi_h = pts[v_i].height;
    Rprintf("v_i: %d, height: %lf\n", v_i, vi_h);
  
    // Union all lower adjacent vertices
    IntegerVector v_adj = adj_v(v_i);
    IntegerVector::const_iterator v_it;
    for (v_it = v_adj.begin(); v_it != v_adj.end(); ++v_it){
      const size_t v_j = size_t(*v_it);
      const size_t vi_z = pts[v_i].z_idx, vj_z = pts[v_j].z_idx;
      if ( pts[v_j].height < vi_h || ds.Find(vi_z) == ds.Find(vj_z)){
        continue;
      }
      Rprintf("vj: %d, height: %lf, cc: %d\n", v_j, pts[v_j].height, ds.Find(vj_z));
      // Else, make edge
      Rprintf("Unioned %d,%d, adding edge %d, %d\n", vi_z, vj_z, v_i, lowest_vertex.at(ds.Find(vj_z)));
      edge[0] = v_i; edge[1] = lowest_vertex.at(ds.Find(vj_z));
      mt_ptr->insert_simplex(edge);
      ds.Union(vi_z, vj_z);

      Rprintf("inserted simplex, updating: lv[%d]=%d (lv.size == %d)\n", ds.Find(vj_z), v_i, lowest_vertex.size());
      size_t vj_cc = ds.Find(vj_z);
      Rprintf("vj_cc:%d\n",vj_cc);
      Rprintf("lv[vj_cc]:%d\n",lowest_vertex.at(vj_cc));
      Rcout << "Assigning" << std::endl;
      lowest_vertex.at(vj_cc) = v_i;
    }
  }
  return(List::create());
}

// Construct the merge tree

// List construct_merge_tree(SEXP smesh, const NumericVector& h, SEXP merge_tree_res) {
//   using idx_v = vector< size_t >;
//   Rcpp::XPtr<SimplexTree> stree_ptr(smesh); // Collect the simplex tree (mesh of vertices)
//   
//   // Sort and order the height function values 
//   const size_t n = h.size(); 
//   if (stree_ptr->n_simplexes[0] != n){ stop("The number of height values must be equal to the number of 0-simplices in the mesh."); }
//   idx_v vids = as< idx_v >(stree_ptr->get_vertices()); 
//   
//   // Generate point maps (sorted by height)
//   size_t i = 0; 
//   map< size_t, id_map > pts = map< size_t, id_map >( );
//   generate_n(std::inserter(pts, pts.begin()), n, [&h, &i, &vids](){
//     id_map pt{ i, h[i] };
//     return(std::make_pair(vids[i++], pt)); 
//   });
//   
//   // Order vertices by height
//   idx_v h_order = order(as< vector<double> >(h));
//   
//   // To extract adjacenct vertices
//   auto adj_v = [&stree_ptr](const size_t v_i){
//     return(stree_ptr->adjacent_vertices(v_i));
//   };
//   // To get the vertices that are of lower height than a given vertex v_i
//   auto lower_v = [&adj_v, &pts](const size_t v_i){
//     idx_v adj = as< idx_v >(adj_v(v_i));
//     adj.erase(std::remove_if(adj.begin(), adj.end(), [&v_i, &pts](const size_t v_j){
//       return( pts[v_i].height < pts[v_j].height );
//     }), adj.end());
//     return(adj);
//   };
//   // Convert IDS to zero-based indices used by disjoint set 
//   auto id_to_zidx = [&pts](const idx_v& ids){
//     idx_v res = idx_v(ids.size());
//     std::transform(ids.begin(), ids.end(), res.begin(), [&pts](const size_t vi){
//       return(pts[vi].z_idx);
//     });
//     return(res);
//   };
//   // To detect whether a point is local min
//   auto is_local_min = [&adj_v, &pts](const int v_i){
//     const IntegerVector j = adj_v(v_i);
//     return(std::all_of(j.begin(), j.end(), [&v_i, &pts](const int v_j){
//       return( pts[v_i].height < pts[v_j].height );
//     }));
//   };
//   
//   // Detect whether a point is a critical point 
//   UnionFind ds = UnionFind(n);
//   // IntegerVector top = ds.getCC(); // Maintain map between CCs and highest vertex id
//   auto is_critical = [&ds, &pts, &lower_v, &id_to_zidx](const int v_i){
//     idx_v j_v = lower_v(v_i); // First, extract lower-height points
//     IntegerVector j_v_int = wrap(j_v);
//     Rcout << "v_i: " << v_i << " adjacencies: " << j_v_int << std::endl;
//     Rcout << "v_i height: " << pts[v_i].height << ", v_j height: ";
//     for (size_t j = 0; j < j_v.size(); ++j){ // short-circuit if more than 2 points lower
//       size_t v_j = static_cast< size_t >(j_v[j]);
//       Rcout << pts[v_j].height << ", ";
//     }
//     Rcout << std::endl;
// 
// 
//     IntegerVector j_comp = ds.FindAll(id_to_zidx(j_v)); // then get their components
//     return(bool(as< IntegerVector >(unique(j_comp)).size() >= 2)); // 2 or more cc's == critical point
//     
//     // std::remove_if(j_v)
//     // const size_t i_cc = ds.Find(pts[v_i].z_idx);
//     // const double i_h = pts[v_i].height;
//     // size_t n_lower = 0, n_j = j_v.size(); 
//     // for (size_t j = 0; n_lower < 2 && j < n_j; ++j){ // short-circuit if more than 2 points lower
//     //   size_t v_j = static_cast< size_t >(j_v[j]);
//     //   if (pts[v_j].height < i_h && i_cc != ds.Find(pts[v_j].z_idx)){
//     //     ++n_lower;
//     //   }
//     // }
//     // return(n_lower >= 2);
//   };
//   
//   // Maintain lower map between extrema points and regular point. This is use to 
//   // generate the super edges later. 
//   map< size_t, idx_v > top; 
//   vector< size_t > crit_map = vector< size_t >(n);
//   std::iota(crit_map.begin(), crit_map.end(), 0);
//   
//   // Prepare result objects 
//   Rcpp::XPtr<SimplexTree> relations(merge_tree_res); // Result
//   auto super_nodes = vector< super_node >(); // super nodes  
//   auto super_edges = map< pair< size_t, size_t >, super_edge >(); // super edges
//   
//   // Iterate via Carr's algorithm 
//   std::vector< size_t > node = { 0 };
//   for (idx_v::iterator v_it = h_order.begin(); v_it != h_order.end(); ++v_it){
//     size_t v_i = vids[*v_it]; 
//     const id_map& vi_info = pts[v_i];
//     
//     // If local min, create new super node
//     if (is_local_min(v_i)){
//       Rprintf("v_%d is a local min\n", v_i);
//       node[0] = v_i;
//       relations->insert_simplex(node);
//       super_node v = { v_i, vi_info.height };
//       super_nodes.push_back(v);
//       top[v_i] = idx_v();
//     } 
//     // Else if critical point, such as a saddle, create a super node + super edge(s)
//     else if (is_critical(v_i)){
//       Rprintf("v_%d is a critical point\n", v_i);
//       node[0] = v_i;
//       relations->insert_simplex(node);
//       super_node v = { v_i, vi_info.height };
//       super_nodes.push_back(v);
//       top[v_i] = idx_v();
//       
//       const idx_v lv = lower_v(v_i);
//       const size_t vi_cc = ds.Find(vi_info.z_idx);
//       std::for_each(lv.begin(), lv.end(), [&crit_map, &ds, &v_i, &vi_cc, &pts]( const size_t v_j ){
//         const size_t vj_cc = ds.Find(pts[v_j].z_idx);
//         if (vi_cc != vj_cc){
//           ds.Union(pts[v_i].z_idx, pts[v_j].z_idx);
//           size_t vi_cc = ds.Find(pts[v_i].z_idx);
//           crit_map.at(vi_cc) = v_i;
//           // super_edge e { vector< size_t >(), vector< double >() };
//           // super_edges.insert( std::make_pair(std::make_pair(vi_cc, vj_cc), e) );
//         }
//       });
//     } // Otherwise, v_i is a regular point, push into appropriate super edge
//     else { 
//       const idx_v lv = lower_v(v_i);
//       //const size_t vi_top = top[ds.Find(vi_info.z_idx)], vj_top = top[ds.Find(pts[lv.at(0)].z_idx)];
//       
//       // Union the two points 
//       ds.Union(vi_info.z_idx, pts[lv.at(0)].z_idx);
//       size_t crit_id = crit_map[ds.Find(pts[lv.at(0)].z_idx)];
//       top[crit_id].push_back(v_i);
//       // crit_map[]
//       // Add point to super edge
//       // super_edge& e = super_edges[std::make_pair(vi_cc, vj_cc)];
//       // e.ids.push_back(v_i);
//       // e.heights.push_back(vi_info.height);
//     }
//   }
//   
//   ds.printCC();
//   
//   // Export 
//   List res = export_to_list(super_nodes, super_edges), top_res = List();
//   for (auto& kv: top){ top_res[std::to_string(kv.first)] = wrap(kv.second); }
//   res["top"] = top_res;
//   return(res);
// }
//   
// Exposed Rcpp Module 
// RCPP_MODULE(merge_tree_module) {
//   Rcpp::class_<MergeTree>("MergeTree")
//   .constructor()
//   .method( "as_XPtr", &MergeTree::as_XPtr)
//   .method( "construct", &MergeTree::construct)
//   .method( "export_to_list", &MergeTree::export_to_list)
//   ;
// }
  
/*** R
set.seed(1234)
x <- c(rnorm(10, mean = -1, sd = 0.5), 
       rnorm(15, sd = 0.5), 
       rnorm(10, mean = 1, sd = 0.5), 
       rnorm(5, mean = 2.5, sd = 0.5))

## Get density estimate
bw <- ks::hpi(x)/2.5
d_x <- ks::kde(x, h = bw)
pt_est <- ks::kde(x, h = bw, eval.points = x)

## Get mesh of data points connecting neighbors + density estimate as height function 
h <- pt_est$estimate
n <- length(h)
x_idx <- order(x)-1L ## gets order of data (w/ domain == support of the density)
stree <- simplex_tree()
invisible(mapply(function(i, j){ stree$insert_simplex(c(i, j)) }, x_idx[1:(n-1)], x_idx[2:n]))
vids <- stree$get_vertices()

## Plot the density; you can see the critical points 
plot(d_x)
rug(x)
points(cbind(x, h), col = "red", cex = 0.5)
text(cbind(x, h), col = "red", labels = vids, pos = 1, cex = 0.45)

## Extract the merge tree
mtree <- Mapper:::MergeTree$new()
mtree$construct(stree$as_XPtr(), h)
test <- mtree$export_to_list()

## Plot critical points
points(cbind(x[match(test$nodes$ids, vids)], test$nodes$heights), col = "green", cex = 0.55)

## To igraph 
vids <- jt$get_vertices()
el <- jt$as_edge_list()
g <- igraph::graph_from_edgelist(t(apply(el, 1, function(e){ match(e, vids) })), directed = FALSE)
igraph::vertex_attr(g, "label") <- vids


## Uses low-poly version of the stanford bunny
bunny <- readobj::read.obj("~/Downloads/Bunny-LowPoly.obj")

layout(matrix(c(1,2), nrow = 1))


## Plot critical points 
crit <- jt$get_vertices()+1L
points(cbind(x, pt_est$estimate)[crit,], col = "green", cex = 0.5)

## This computes the critical points 
# cenv <- new.env(parent = emptyenv())
# cenv$critical_pts <- c()
# uf <- Mapper::union_find(length(h))
# for (idx in (order(h)-1L)){
#   adj <- stree$adjacent_vertices(idx)
#   is_lower <- h[idx+1L] < h[adj+1L]
#   is_local_min <- all(is_lower)
#   if (is_local_min){
#     cenv$critical_pts <- c(cenv$critical_pts, idx)
#     next
#   }
#   ## If two of more neighbors of v_i are in different CCs and their heights are lower than h_i, 
#   ## then v_i is a saddle point, under the assumption of linearity between simplexes. 
#   is_saddle <- sum((uf$find(idx) != uf$find_all(adj)) & (!is_lower)) >= 2
#   uf$union_all(c(idx, adj))
#   if (is_saddle){
#     cenv$critical_pts <- c(cenv$critical_pts, idx)
#     next
#   }
# }

## Plot 
plot(d_x)
rug(x)
points(cbind(x, pt_est$estimate), col = "red", cex = 0.5)
text(cbind(x, pt_est$estimate), col = "red", labels = seq(length(x)), pos = 3, cex = 0.45)

crit_idx <- unname(id_to_idx[cenv$critical_pts])+1L
points(cbind(x, pt_est$estimate)[crit_idx,], col = "green", cex=0.55)

# wut <- Mapper:::join_tree(pt_est$estimate, stree$as_XPtr())
# points(, col = "red", cex = 0.5)

*/
//   
// }
// std::for_each(pts.begin(), pts.end(), [](const id_map& v){
//   
//   if (is_local_min(v.)){
//     
//   }
// });
// 
// 
// 
// // std::vector<double> h_ = as< std::vector<double> >(h);
// // std::vector<size_t> idx = order(h_); // indices start from 0, ordered by height
// 
// // Map between 0-based index and 
// // std::map<int, int> idx_to_id = std::map<int, int>(), id_to_idx = std::map<int, int>(); 
// // const IntegerVector vids = stree_ptr->get_vertices();
// // for (size_t i = 0; i < n, ++i){ idx_to_id[i] = vids[i]; }
// // for (size_t i = 0; i < n, ++i){ id_to_idx[vids[i]] = idx_h.first[i]; }
// 
// // Disjoint set 
// UnionFind uf = UnionFind(n); // track connected components
// // IntegerVector cc = uf.getCC();
// 
// 
// 
// 
// // Get the join tree reference
// // The join tree will use the ids in the simplicial mesh to label the super nodes.
// // and super edges. A 0-based index map between these labels and a disjoint-set is required
// // to track components efficiently.  
// SimplexTree& join_tree = *jt_res_ptr;
// std::vector<uint> super_node = std::vector<uint>(1), super_edge = std::vector<uint>(2);
// 
// 
// // std::vector<int> top = as< std::vector<int> >(uf.getCC());
// 
// // Ascend 
// for (size_t i = 0; i < n; ++i){
//   size_t v_i = idx[i];
//   
//   // If v_i is a local minimum, create a new super node in the graph
//   if (is_local_min(v_i)){
//     Rprintf("Local min: %d\n", v_i);
//     super_node[0] = v_i;
//     join_tree.insert_simplex(super_node);
//   }
//   // Else if v_i is a critical point (but not a minima), create a new super node
//   // in the tree, and connect it to the prior two components
//   else if (is_critical(v_i)){
//     Rprintf("Critical merge: %d\n", v_i);
//     // Create merge node 
//     super_node[0] = v_i;
//     join_tree.insert_simplex(super_node);
//     
//     // Connect to other components
//     const IntegerVector j_v = stree_ptr->adjacent_vertices(v_i);
//     std::for_each(j_v.begin(), j_v.end(), [&uf, &v_i, &join_tree, &super_edge, &top, &h](const int v_j){
//       if (uf.Find(v_j) != uf.Find(v_i) && h[v_j] < h[v_i] ){
//         super_edge[0] = v_i, super_edge[1] = top[uf.Find(v_j)];
//         uf.Union(v_i, v_j);
//         Rprintf("Branch: %d --> %d\n", super_edge[0], super_edge[1]);
//         join_tree.insert_simplex(super_edge);
//         
//         // Update top
//         top[uf.Find(v_i)] = v_i;
//       }
//     });
//   } else {
//     IntegerVector j_v = stree_ptr->adjacent_vertices(v_i);
//     std::for_each(j_v.begin(), j_v.end(), [&uf, &v_i, &join_tree, &h, &top](const int v_j){
//       if (h[v_j] < h[v_i]){ 
//         int c_min = top[uf.Find(v_j)];
//         uf.Union(v_i, v_j); 
//         top[uf.Find(v_i)] = top[uf.Find(v_j)] = c_min; // Always mark top with *lowest* vertex
//       }
//     });
//   }
// }


// for (int i = 0; i < n; ++i){
//   const int v_i = idx_to_id[idx_h.first[i]]; // 0-based index into h
//   const double h_i = h[id_to_idx[v_i]];
//   IntegerVector adj_j = stree_ptr->adjacent_vertices(v_i);
//   std::for_each(adj_j.begin(), adj_j.end(), [&uf, &h_i](const int v_j){
//     const double h_j = h[id_to_idx[v_j]];
//     if (h_j < h_i){
//       const int i_comp = uf.Find(id_to_idx[v_i]); 
//       const int j_comp = uf.Find(id_to_idx[v_j]); 
//       if (i_comp != j_comp){
//         uf.Union(id_to_idx[v_i], id_to_idx[v_j]);
//         
//         // New merge node. 
//         IntegerVector cc = uf.getCC();
//         
//         // Create new merge node with id == v_i
//         cc_heighest_v[uf.Find(id_to_idx[v_j])];
//       }
//     }
//       
//     
//       uf.Union(v_i, j);
// 
//       cc_heighest_v[v_i], uf.Find(idx);
//     }
//   });
//   IntegerVector cc = uf.getCC();
// }

// BAD IDEA
// Map between vertex ids and height function values
// std::map<double, int> id_to_h = std::map<double, int>();
// const IntegerVector vids = stree_ptr->get_vertices();
// if (vids.size() != h.size()){ stop("h length != number of 0-simplices in mesh."); }
// for (size_t i = 0; i < vids.size(), ++i){ id_to_h[h[i]] = vids[i]; }


// Constructor
// MergeTree::MergeTree(){
//   relations = SimplexTree();
//   super_nodes = vector< super_node >(); // super nodes  
//   super_edges = map< pair< size_t, size_t >, super_edge >(); // super edges
// }
// MergeTree::~MergeTree(){ }

// Export as XPtr
// SEXP MergeTree::as_XPtr(){
//   Rcpp::XPtr< MergeTree> p(this, false);
//   return(p);
// }