# ## Preprocessing
# X_pt <- as.matrix(readr::read_delim(file = "~/GHmatch/GHMatch/nonrigid3d/horse5.vert", delim = " ", col_names = FALSE))
# Y_pt <- as.matrix(readr::read_delim(file = "~/GHmatch/GHMatch/nonrigid3d/horse10.vert", delim = " ", col_names = FALSE))
# 
# ## Compute the euclidean distances
# X_dist <- parallelDist::parallelDist(X_pt)
# Y_dist <- parallelDist::parallelDist(Y_pt)
# 
# ## Make the eps-ball graph 
# eps_graph <- function(x){
#  mst <- dbscan:::prims(x, n = attr(x, "Size"))
#  eps <- max(mst[,3]) + max(mst[,3])*0.01 
#  x[x > eps] <- Inf # Make paths farther than max-eps unreachable
#  return(x)
# }
# X_eps <- eps_graph(X_dist)
# Y_eps <- eps_graph(Y_dist)
# 
# ## Make the graph(s)
# X_graph <- igraph::graph_from_adjacency_matrix(as.matrix(X_eps), mode = "undirected", weighted = TRUE)
# Y_graph <- igraph::graph_from_adjacency_matrix(as.matrix(Y_eps), mode = "undirected", weighted = TRUE)
# 
# ## Choose number of prototypes 
# n <- 10L
# X_proto <- Mapper::landmarks(X_pt, n = n)
# Y_proto <- Mapper::landmarks(Y_pt, n = n)
# 
# ## Compute the shortest paths 
# D1 <- as.dist(igraph::distances(X_graph, v = X_proto, to = X_proto, algorithm = "dijkstra"))
# D2 <- as.dist(igraph::distances(Y_graph, v = Y_proto, to = Y_proto, algorithm = "dijkstra"))
# 
# ## Compute the GH distance
# gamma <- Mapper::gromov_hausdorff(D1, D2)
# 
# 
# index_lt <- function(from, to, N){
#  from <- max(c(from, to)); to <- min(c(from, to))
#  return((N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1))
# }
# sigma <- 4
# mu <- 8
# iter <- 15L
gh_match <- function(D1, D2, iter, sigma, mu){
 stopifnot("dist" %in% class(D1) && "dist" %in% class(D2))
 stopifnot(attr(D1, "Size") == attr(D2, "Size"))
 n <- attr(D1, "Size")
 
 ## gromov-hausdorff distance
 gamma <- Mapper:::gromov_hausdorff(D1, D2)
 
 ## Constraints 
 A <- matrix(0L, nrow=2*n, ncol=n^2)
 for (i in 1:n){
  A[i,((i-1)*n+1):(i*n)] <- 1L
  for (j in 1:n){ A[i+n, (j-1)*n+i] <- 1L }
 }
 b <- matrix(1L, nrow = 2*n, ncol = 1) 
 m <- nrow(A)  # number of constraints
 
 ## lagrangian
 lagrangian <- function(C, A, b, Y, lambda, sigma){
  f <- function(Y){ 
   sum(diag(C %*% (Y %*% t(Y))))-t(lambda) %*% (A %*% Y-b) + 0.5*as.vector(sigma)*norm(A %*% Y-b, "2")^2 
  }
  g <- function(Y){ 2.0 * C %*% Y - t(A) %*% lambda + as.vector(sigma) * t(A) %*% (A %*% Y-b) }
  list(f=f,g=g)
 }
 
 # initial guess for Y
 Y <- matrix(1/2, nrow = n^2, ncol = 1);
 
 # initial multiplier
 lambda <- matrix(0L, nrow = m, ncol = 1);
 N <- n^2
 
 # store the iterations
 Z <- matrix(0L, nrow=N, ncol=iter)
 Z[, 1L] <- Y
 
 # feasibility
 feas <- matrix(0L, nrow=iter)
 feas[[1]] <- norm(A %*% Y-b, type = "2")
 obj <- matrix(0L, nrow=iter)
 obj[[1]]  <- sum(diag(gamma %*% Y %*% t(Y)))
 
 ## optimizations options
 opt_control <- nloptr::nl.opts(list(
  algorithm = "NLOPT_LD_MMA",
  xtol_rel = 1.0e-4,
  print_level = 2,
  check_derivatives = TRUE,
  check_derivatives_print = "all"
 ))
 
 # Run the optimization 
 lb <- matrix(0, N, 1); ub <- matrix(1, N, 1)
 for (i in 2L:iter){
  obj_f <- lagrangian(gamma, A, b, Y, lambda, sigma)
  new_Y <- nloptr::auglag(
   x0 = Y, fn = obj_f$f, gr = obj_f$g,
   lower = lb, upper = ub, 
   localsolver = "LBFGS", localtol = 1e-4, nl.info = TRUE, 
   opts = list(),
  )
  
  # update feasibility
  feas[i,1] <- norm(A %*% Y - b, "2")
  
  # update multipliers
  lambda <- lambda - as.vector(sigma) * (A %*% Y - b)
  sigma <- mu * sigma
  Z[,i] <- new_Y$par
  obj[i,1] <- sum(diag(gamma %*% Y %*% t(Y)))
 }
 
}