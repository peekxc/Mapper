#' gromov hausdorff 
#' @description This function implements the discrete, extended version of gromov hausdorff distance proposed by 
#' Memoli in section 7 of (1) for comparing two measure-metric spaces.
#' @param d_X either a dist object, or (n x n) numeric metric distance matrix over \code{X}. 
#' @param d_Y either a dist object, or (m x m) numeric metric distance matrix over \code{Y}. 
#' @param mu_X (n-length) numeric vector giving the probability measure at each point of \code{X}. 
#' @param mu_Y (m-length) numeric vector giving the probability measure at each point of \code{Y}.  
#' @param flb_only If only crude approximation is needed, the lower bound can be returned. See details.
#' @param flb_solver Which LOP ROI plugin solver to use to solve the lower bound approximation.
#' @param roi_opt Outer options to pass to ROI_solve
#' @param control Inner options to pass to ROI_solves 'control' parameter. 
#' @references 1. Mémoli, Facundo. "On the use of Gromov-Hausdorff distances for shape comparison." (2007).
#' @details This function provides an implementation of the computational objective \emph{(Pp)} listed in Section 7 of (1), which 
#' calculates an Lp approximition of Gromov-Hausdorff (i.e. wasserstein relaxation) distance between two measure-metric spaces. The first lower bound (FLB)
#' is first computed via an LOP, and then the solution is used as the initialization point for the corresponding quadratic optimization 
#' routine. This code utilizes the \emph{R Optimization Infrastructure}(\pkg{ROI}), and relies on having available at least one plugin 
#' available for to solve the initial lower bound (a linear program) and the final optimization (cast as a generic nonlinear optimization). 
#' Both solutions are returned to the user, unless \code{flb_only} is set to true, in which case only the FLB is returned. 
#' 
#' \strong{NOTE:} The primary quadratic optimization problem (QOP) is non-convex, and requires solving for (n^2 x m^2) variables.
#' Thus, the function may be very computationally expensive to run. The FLB is an LOP of the same size. 
#' 
#' @export 
gromov_hausdorff <- function(d_X, d_Y, mu_X, mu_Y, p = 1,
                             flb_solver="glpk", flb_only=FALSE, 
                             options = list(localsolver="COBYLA", localtol = 1e-7), 
                             control = list(maxeval=300), 
                             return_optimizer = FALSE){
  roi_installed <- requireNamespace("ROI", quietly = TRUE)
  if (!roi_installed){ stop("This function requires the R Optimization Infrastructure (ROI) package to be installed.") }
  stopifnot((sum(mu_X) - 1.0) < sqrt(.Machine$double.eps))
  stopifnot((sum(mu_Y) - 1.0) < sqrt(.Machine$double.eps))
  if (is.matrix(d_X)) { stopifnot(dim(d_X)[[1]] == dim(d_X)[[2]]) }
  else { stopifnot(is(d_X, "dist")); d_X <- as.matrix(d_X) }
  if (is.matrix(d_Y)) { stopifnot(dim(d_Y)[[1]] == dim(d_Y)[[2]]) }
  else { stopifnot(is(d_Y, "dist")); d_Y <- as.matrix(d_Y) }
  
  ## Useful variables
  { n_x <- nrow(d_X); n_y <- nrow(d_Y) }
  idx <- matrix(seq(n_x * n_y), nrow = n_x, ncol = n_y)
  
  ## Make constraints
  A <- gh_make_A(idx-1L)
  constraints <- ROI::L_constraint(L = A, dir = rep("==", n_x + n_y), rhs = c(mu_X, mu_Y))
  
  ## Make objective for FLB
  s <- function(d_X, p, lambda){ apply(d_X, 1, function(d_i){ sum((d_i^p) * lambda)^(1/p) }) }
  sp_X <- s(d_X, p = 1, lambda = mu_X)
  sp_Y <- s(d_Y, p = 1, lambda = mu_Y)
  flb_objective <- ROI::L_objective(as.vector(outer(1L:n_x, 1L:n_y, FUN = Vectorize(function(i,j){ abs(sp_X[i] - sp_Y[j]) }))))
  
  ## Bounds on mu
  mu_bnds <- ROI::V_bound(li = seq(n_x*n_y), ui = seq(n_x*n_y), lb = rep(0, n_x*n_y), ub = rep(1,n_x*n_y))
  
  ## The lower-bound -- as a LOP
  lop <- ROI::OP(objective = flb_objective, constraints = constraints, bounds = mu_bnds)
  
  ## Need to have a solver available
  available_solvers <- ROI::ROI_applicable_solvers(lop)
  if (is.null(available_solvers)){
    stop("No applicable solvers found with ROI! The FLB is a linear optimization.")
  } 
  if (!missing(flb_solver)){ stopifnot(flb_solver %in% available_solvers) }
  
  ## Solve the first lower bound
  flb <- ROI::ROI_solve(lop, solver = flb_solver)
  if (flb_only) { return(wrap_solution(flb, idx, return_optimizer)) }
  
  ## Make the Q matrix of distances
  Q <- gh_make_Q(d_X, d_Y)
  
  ## Setup objective function(s), gradients, etc. 
  N <- n_x * n_y
  f <- function(mu){
    as.vector({ structure(mu, dim=c(1L, N)) %*% Q %*% structure(mu, dim=c(N, 1L)) })
  }
  grad_f <- function(mu){ Q %*% matrix(mu, ncol = 1) }
  heq <- function(mu){ (A %*% mu) - c(mu_X, mu_Y) } # rowSums(relist(mu, idx)) - c(mu_X, mu_Y)
  # heqjac <- function(x) { nloptr::nl.jacobian(x, heq) }
  
  ## Minor fix to the FLB bounds to make optimization more stable. This usually isn't necessary. 
  flb$solution[flb$solution < 0] <- 0.0
  flb$solution[flb$solution > 1] <- 1.0
  
  ## Set solver options
  { lb <- rep(0, n_x*n_y); ub <- rep(1, n_x*n_y) }
  default_opts <- list(x0=flb$solution, fn = f, lower = lb, upper = ub, heq = heq)
  default_opts$control <- modifyList(list(maxeval=300, stopval=flb$objval), control)
  options <- modifyList(default_opts, val = options)
  
  ## Attach gradient function if COBYLA isn't used
  if (options$localsolver %in% c("LBFGS", "MMA", "SLSQP")){
    options$gr <- grad_f
  }

  ## Use nloptr augmented lagrangian approach
  al_res <- do.call(nloptr::auglag, options)

  ## Return results
  return(list(lop_res = wrap_roi_solution(flb, idx, return_optimizer), 
              qop_res = wrap_nloptr_solution(al_res, idx, return_optimizer)))
}

wrap_nloptr_solution <- function(res, idx, return_opt){
  tmp <- list(
    gh = res$value,
    mu = res$par, 
    correspondences = list(
      xy = apply(relist(res$par, idx), 1, which.max),
      yx = apply(relist(res$par, idx), 2, which.max)
    )
  )
  if (return_opt){ tmp$optim_status <- res[c("global_solver", "local_solver", "convergence", "message")] }
  return(tmp)
}

wrap_roi_solution <- function(res, idx, return_opt){
  tmp <- list(
    gh = res$objval, 
    mu = res$solution,
    correspondences = list(
      xy = apply(relist(res$solution, idx), 1, which.max),
      yx = apply(relist(res$solution, idx), 2, which.max)
    )
  )
  if (return_opt){ tmp$optim_status <- res[c("status", "message")] }
  return(tmp)
}
  
## Gromov Hausdorff objective module (can be loaded anywhere)
Rcpp::loadModule("gh_module", TRUE)

## Extra code 
# A <- outer(X = seq(n_x + n_y), Y = seq(n_x * n_y), FUN = Vectorize(function(i, j){
#     if (i <= n_x){ ifelse(j %in% idx[i,], 1L, 0L) }
#     else { ifelse(j %in% idx[,(i-n_x)], 1L, 0L) }
# }))
## Make the control
# control <- list(
#   algorithm = "NLOPT_LN_AUGLAG_EQ",
#   maxeval = 10L,
#   local_opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-4, stopval = 0, heq = heq, heqjac = heqjac),
#   print_level = 0,
#   check_derivatives = FALSE
# )
# The indices to use to compute the distances, filling Q row-wise
# rhs_idx <- cbind(as.vector(row(idx)), as.vector(col(idx)))
# master_idx <- do.call(rbind, lapply(seq(nrow(rhs_idx)), function(i){
#   cbind( t(replicate(n = nrow(rhs_idx), rhs_idx[i,])), rhs_idx)
# }))
## The distances to use
# Q <- matrix(apply(master_idx, 1, function(ijkl){
#   { i <- ijkl[1]; j <- ijkl[2]; k <- ijkl[3]; l <- ijkl[4] }
#   abs(d_X[i,k] - d_Y[j, l])
# }), nrow = n_x*n_y, ncol = n_x*n_y, byrow = TRUE)
## Solve the QOP -- using general non-linear optimizer, since there is no non-convex solvers available
## -- Unfortunately, the interface isn't quite as flexible as one would hope
# q_obj <- ROI::Q_objective(Q = Q)
# qop <- ROI::OP(objective = q_obj, constraints = constraints, bounds = mu_bnds)
# qop_res <- ROI::ROI_solve(qop, solver = "nloptr", x0 = flb$solution, control = control)
