

invisible(sapply(c("ROI", "ROI.plugin.glpk", "nloptr"), function(lib) { library(lib, character.only = TRUE) }))
rotate <- function(X, theta){
  theta <- theta*pi/180 # radians
  r <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  t(r %*% t(X))
}
reflect <- function(X, r_axis=0){
  if (r_axis == 0){ return(cbind(X[,1], -X[,2])) }
  if (r_axis == 1){ return(cbind(-X[,1], X[,2])) }
}

add_noise <- FALSE
n <- 45L
X <- replicate(2, rnorm(n))
eps <- min(dist(X))
if (add_noise){
  Y <- rotate(X, 45) + runif(n, min = 1.5*eps, max = eps*2.5)
} else {
  Y <- rotate(X, 45)
}

{ d_X <- as.matrix(dist(X)); d_Y <- as.matrix(dist(Y)) }
mu_X <- rep(1/nrow(X), nrow(X))
mu_Y <- rep(1/nrow(Y), nrow(Y))
res <- Mapper::gromov_hausdorff(d_X = d_X, d_Y = d_Y, mu_X = mu_X, mu_Y = mu_Y)

## Plot the points + translated Y
diam_X <- max(apply(d_X, 2, max))
translated_Y <- cbind(Y[,1] + diam_X, Y[,2])
plot(rbind(X, translated_Y), pch = 20, col = c(rep("black", nrow(X)), rep("red", nrow(Y))), asp = 1)


text(X, labels = 1L:n, col = "black", pos = 1)
text(translated_Y, labels = seq(nrow(Y)), col = "red", pos = 1)

idx <- matrix(seq(nrow(X) * nrow(Y)), nrow = nrow(X), ncol = nrow(Y))
smpl <- sample(1:nrow(X), size = 15L)
correspondence <- res$qop_res$correspondences$xy
for ( i in smpl){
  segments(X[i,1], X[i,2], translated_Y[correspondence[i],1], translated_Y[correspondence[i],2], col = "blue")
}




## Mix up Y 
rand_idx <- sample(1:nrow(Y))
Y <- Y[rand_idx,]
d_Y <- as.matrix(dist(Y))
rand_idx
apply(relist(flb$solution, idx), 1, which.max)
apply(relist(soln$solution, idx), 1, which.max)


## Plot both, and the matching 
plot(rbind(X, Y), asp = 1, col = rep("black", nrow(X)), rep("red", nrow(Y)))
idx <- matrix(seq(nrow(X) * nrow(Y)), nrow = nrow(X), ncol = nrow(Y))
correspondence <- apply(relist(res$solution, idx), 1, which.max)
for (i in 1:nrow(X)){
  segments(X[i,1], X[i,2], Y[correspondence[i],1], Y[correspondence[i],2], col = "blue")
}

make_counter <- function(start, end, reset=TRUE){
  cc <<- start
  return(function(){
    if (cc >= end){ if(reset) cc <<- start else return(end) }
    else { cc <<- cc+1L; return(cc-1L) }
  })
}


centroid <- apply(X, 2, min) + apply(X, 2, function(x) diff(range(x)))/2
diameter <- max(dist(X))
eps <- diameter * 0.05
fps <- 30




test_f <- function(...){
  args <- pryr::dots(...)
  eval(args[[1]])
}
test_f({ print("hello") }, {print(x)})


## Testing rotation animations
library(ggvis)

ui <- shinyUI(pageWithSidebar(headerPanel = headerPanel("GGvis"),
  sidebarPanel(
    radioButtons("rotate", label = "Rotate", choices = c("1", "2"))
  ),
  mainPanel(
    # uiOutput("ggvis_ui"),
    ggvisOutput("ggvis"), 
    tags$script('
    $(document).on("keydown", function (e) {
       Shiny.onInputChange("keydown", e.which);
    });
    ') 
  )
))
server <- shinyServer(function(input, output, session) {
  counter <- make_counter(0, 90, reset = FALSE)
  rv <- reactiveValues(keydown=0L)
  X_anim <- reactive({
    if (rv$keydown == 0){
      structure(data.frame(X), names = c("V1", "V2"))
    } else if (rv$keydown == 1){
      invalidateLater(2000/90, NULL)
      step <- counter()
      structure(data.frame(rotate(X, step)), names = c("V1", "V2"))
    } else if (rv$keydown == 2){
      
    }
  })
  X_anim %>% ggvis(~V1, ~V2) %>% layer_points() %>% 
    scale_numeric("x", domain = c(centroid[1]-diameter/2-eps, centroid[1]+diameter/2+eps)) %>% 
    scale_numeric("y", domain = c(centroid[2]-diameter/2-eps, centroid[2]+diameter/2+eps)) %>% 
    add_axis("x", title = "") %>% add_axis("y", title = "") %>%
    set_options(duration = 0) %>% 
    bind_shiny(plot_id = "ggvis")
  observeEvent(input$keydown, {
    rv$keydown <- rv$keydown+1L
    # print("hello")
  })
})
shiny::shinyApp(ui, server)

animation::saveGIF({
  X_anim %>% ggvis(~V1, ~V2) %>% layer_points() %>% 
    scale_numeric("x", domain = c(centroid[1]-diameter/2-eps, centroid[1]+diameter/2+eps)) %>% 
    scale_numeric("y", domain = c(centroid[2]-diameter/2-eps, centroid[2]+diameter/2+eps)) %>% 
    add_axis("x", title = "") %>% add_axis("y", title = "") %>%
    set_options(duration = 0)
}, img.name = "rotate_90.gif")


# 
# 
# # wut <- lpSolveAPI::get.objective(lp_program)
#   # lpSolveAPI::get.variables(lp_program)
#   # outer(X = seq(n_x + n_y), Y = seq(n_x * n_y), FUN = Vectorize(function(i, j){
#   #   lpSolveAPI::get.mat(lp_program, i, j)
#   # }))
#   
#   x <- replicate(2, rnorm(10))
#   theta <- 1
# 
#   
#   reactive()
#   
#   x <- data.frame(x1=rnorm(10), x2=rnorm(10))
#   #heta <- input_slider(min = 0, max = 90, value = 0, animate = TRUE)
#   r <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
#   theta_slider <- input_slider(min = 0, max = 90, value = 0, animate = TRUE)
#   x %>% ggvis( props(.props = list(theta = input_slider(min = 0, max = 90, value = 0))) ) %>% 
#     layer_points(props(x = ~x1 + theta,  y = ~x2))
#   f <- function(theta) {
#    
#     
#     x %>% ggvis(x = ~x1, y = ~x2) %>% layer_points()
#   }
# 
#   ggvis::prop()
# 
#   x %>% ggvis::ggvis() %>% ggvis::layer_points()
# 
# }
# 
# 
# 
allPermMatrices <- function(n){
  x <- array(0L, dim = c(n, n, factorial(n)))
  permutations <- rbind(seq(n), permute::allPerms(n))
  for (i in 1:dim(x)[3]){ x[,,i] <- diag(n)[permutations[i,],] }
  return(x)
}

## Test case
{ x <- replicate(2, rnorm(5)); y <- replicate(2, rnorm(5)) }
{ dist_x <- as.matrix(dist(x)); dist_y <- as.matrix(dist(y)) }
R <- allPermMatrices(nrow(x))
max_gamma <- function(dx, dy, r){
  res <- 0
  for (i in 1L:nrow(dx)){
    for (j in 1L:nrow(dy)){
      for (k in 1L:nrow(dx)){
        for (l in 1L:nrow(dy)){
          tmp <- r[i,j]*r[k,l]*abs(dx[i, k] - dy[j, l])
          if (tmp > res){ res <- tmp }
        }
      }
    }
  }
  return(res)
}
r_dist <- sapply(seq(dim(R)[3]), function(i){ max_gamma(dist_x, dist_y, R[,,i]) })
gh_dist <- (1/2)*r_dist[which.min(r_dist)]
max_gamma(dist_x, dist_y, R[,,2])

idx_p <- function(dx, dy, r){
  indices <- list()
  for (i in 1L:nrow(dx)){
    for (j in 1L:nrow(dy)){
      for (k in 1L:nrow(dx)){
        for (l in 1L:nrow(dy)){
          if (r[i,j]*r[k,l] == 1){
            indices <- append(indices, list(data.frame(x1=i,x2=k,y1=j,y2=l)))
          }
        }
      }
    }
  }
  return(indices)
}

t_y <- cbind(y[,1]+3L, y[,2])

animation::saveGIF({
  for (i in seq(dim(R)[[3]])){
    plot(rbind(x, t_y), col = c(rep("blue", nrow(x)), rep("red", nrow(y))), pch=20, cex=1.75, asp=1, 
         xlab = "", ylab="", main = sprintf("R_i = %d",i))
    idx <- as.matrix(do.call(rbind, idx_p(dist_x, dist_y, R[,,i])))
    apply(idx, 1, function(ii){
      x_seg <- x[c(ii[1],ii[2]),]
      y_seg <- t_y[c(ii[3],ii[4]),]
      segments(x0=x_seg[1,1], x1=x_seg[2,1], y0=x_seg[1,2], y1=x_seg[2,2], col=adjustcolor("blue", alpha.f = 0.50))
      segments(x0=y_seg[1,1], x1=y_seg[2,1], y0=y_seg[1,2], y1=y_seg[2,2], col=adjustcolor("red", alpha.f = 0.50))
      segments(x0=x_seg[1,1], x1=y_seg[1,1], y0=x_seg[1,2], y1=y_seg[1,2], col=adjustcolor("orange", alpha.f = 0.50))
      segments(x0=x_seg[2,1], x1=y_seg[2,1], y0=x_seg[2,2], y1=y_seg[2,2], col=adjustcolor("orange", alpha.f = 0.50))
    })
  }
}, interval=0.20)

animation::saveGIF({
  for (i in seq(dim(R)[[3]])){
    plot(grid.arrange(tableGrob(R[,,i])))
  }
}, interval = 0.20)






# 
# 
# gh_dist_approx <- (1/2) * res$opt_result$objective
# 
# ## Objective function
# # make_obj_f <- function(X_dist, Y_dist){
# #   cc <- 1L
# #   function(mu){
# #     # print(cc)
# #     c_sum <- 0
# #     for (i in 1L:n_x){
# #       for (i_p in 1L:n_x){
# #         for (j in 1L:n_y){
# #           for (j_p in 1L:n_y){
# #             { c_idx <- idx[i,j]; c_idx2 <- idx[i_p,j_p] }
# #             c_sum <- c_sum + mu[c_idx]*mu[c_idx2]*abs(d_X[i,j]-d_Y[i_p,j_p])
# #           }
# #         }
# #       }
# #     }
# #     cc <<- cc + 1L
# #     return(c_sum)
# #   }
# # }


rotate <- function(X, theta){
  theta <- theta*pi/180 # radians
  r <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  t(r %*% t(X))
}

## Example comparing a rigid isomorphism
n <- 10L
mixup_Y <- FALSE
Y_idx <- if (mixup_Y)  sample(seq(n)) else seq(n) 
{ X <- replicate(2, rnorm(n)); Y <- rotate(X, 90)[Y_idx,] }
{ d_X <- as.matrix(dist(X)); d_Y <- as.matrix(dist(Y)) }
{ mu_X <- rep(1/nrow(X), nrow(X)); mu_Y <- rep(1/nrow(Y), nrow(Y)) }

## Plot the points + translated Y
diam_X <- max(apply(d_X, 2, max))
translated_Y <- cbind(Y[,1] + diam_X, Y[,2])
plot(rbind(X, translated_Y), pch = 20, col = c(rep("black", nrow(X)), rep("red", nrow(Y))))
text(X, labels = 1L:n, col = "black", pos = 1)
text(translated_Y, labels = Y_idx, col = "red", pos = 1)

## GH distance should be close to 0
soln <- Mapper::gromov_hausdorff(d_X, d_Y, mu_X, mu_Y)

## Plot the mapping between X and Y
x_to_y <- soln$qop_res$correspondences$xy
for (i in 1:nrow(X)){
  segments(X[i,1], X[i,2], translated_Y[x_to_y[i],1], translated_Y[x_to_y[i],2], col = "blue")
}

## Examples with n != m
eps <- min(dist(X))
{ X <- replicate(2, rnorm(n)); Y <- rbind(rotate(X, 90) + runif(n, max = 0.5*eps), rotate(X, 90) + runif(n, max = 0.25*eps))  }
{ d_X <- as.matrix(dist(X)); d_Y <- as.matrix(dist(Y)) }
{ mu_X <- rep(1/nrow(X), nrow(X)); mu_Y <- rep(1/nrow(Y), nrow(Y)) }

## Plot the points + translated Y
diam_X <- max(apply(d_X, 2, max))
translated_Y <- cbind(Y[,1] + diam_X, Y[,2])
plot(rbind(X, translated_Y), pch = 20, col = c(rep("black", nrow(X)), rep("red", nrow(Y))))
text(X, labels = 1L:n, col = "black", pos = 1)
text(translated_Y, labels = 1L:n, col = "red", pos = 1)

## Plot the mapping between X and Y
soln <- Mapper::gromov_hausdorff(d_X, d_Y, mu_X, mu_Y)
idx <- matrix(seq(nrow(X) * nrow(Y)), nrow = nrow(X), ncol = nrow(Y))
x_to_y <- apply(relist(soln$solution, idx), 1, which.max)
for (i in 1:nrow(X)){
  segments(X[i,1], X[i,2], translated_Y[x_to_y[i],1], translated_Y[x_to_y[i],2], col = "blue")
}

## Plot the mapping between Y and X
y_to_x <- apply(relist(soln$solution, idx), 2, which.max)
for (i in 1:nrow(Y)){
  segments(translated_Y[i,1], translated_Y[i,2], X[y_to_x[i],1], X[y_to_x[i],2], col = "blue")
}


## ----- Testing different optimizers -----

## Select all derivative-free optimizers to use with the augmented-lagrangian
opt <- nloptr::nloptr.get.default.options()
optimizers <- grep(x = strsplit(opt$possible_values[[1]], ", ")[[1]], pattern = ".*_[L|G][N|D]_.*", value = TRUE)
optimizers <- optimizers[!optimizers %in% c("NLOPT_LN_AUGLAG", "NLOPT_LN_AUGLAG_EQ")]

## Generate noisy data that *should* still have a perfect bijection, where min_f is known
n <- 10L
X <- replicate(2, rnorm(n))
min_f <- min(dist(X))/2
noise <- runif(n, max = min_f)
rand_idx <- sample(seq(n))
Y <- rotate(X[rand_idx,], 90) + noise
{ d_X <- as.matrix(dist(X)); d_Y <- as.matrix(dist(Y)) }
{ mu_X <- rep(1/nrow(X), nrow(X)); mu_Y <- rep(1/nrow(Y), nrow(Y)) }

res <- lapply(optimizers, function(opt_f){
  control <- list(
    algorithm = "NLOPT_LN_AUGLAG_EQ",
    maxeval = 10L,
    local_opts = list(algorithm = opt_f, xtol_abs = 1e-4, stopval = 0),
    print_level = 3, 
    check_derivatives = TRUE
  )
  Mapper::gromov_hausdorff(d_X, d_Y, mu_X, mu_Y, control = control, return_optimizer = TRUE)
})
sapply(res, function(x) x$gh)

# mu_idx <- matrix(seq(n_x*n_y), nrow = n_x, ncol = n_y)



library('nloptr')
f <- function(x) {
  return( 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2 )
}
f.gradient <- function(x) {
  return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
             200 * (x[2] - x[1] * x[1])) )
}
x <- ROI::OP( objective = ROI::F_objective(f, n = 2L, G = f.gradient),
              bounds = ROI::V_bound(ld = -3, ud= 3, nobj = 2L) )
nlp <- ROI::ROI_solve(x, solver = "nloptr", control = list(x0 = c(-1.2, 1), algorithm="NLOPT_LD_MMA"))

make_f <- function(Q, sep=FALSE){
  if (sep){
    return(list(
      "objective" = function(mu){ as.vector(matrix(mu, nrow = 1) %*% Q %*% matrix(mu, ncol = 1)) }, 
      "gradient" = function(mu){ Q %*% matrix(mu, ncol = 1) }
    ))
  }
  return(function(mu){
    list("objective" = as.vector(matrix(mu, nrow = 1) %*% Q %*% matrix(mu, ncol = 1)),
         "gradient" = Q %*% matrix(mu, ncol = 1))
  })
}

qf <- make_f(Q, sep = TRUE)
alt_obj <- ROI::F_objective(f, n = n_x * n_y, G = grad_f)
q_qut <- ROI::Q_objective(Q)

f_obj <- ROI::F_objective(qf$objective, n = n_x * n_y, G = qf$gradient)

x <- ROI::OP( objective = f_obj,
              bounds = mu_bnds, constraints = constraints )
nlp <- ROI::ROI_solve(x, solver = "nloptr", 
                      control = list(x0 = flb$solution, algorithm="NLOPT_LD_MMA", 
                                     eval_g_eq = function(mu){ }))


test <- nloptr::auglag(x0 = flb$solution, fn = qf$objective, gr = qf$gradient,
                       lower = rep(0L, n_x * n_y), upper = rep(1L, n_x + n_y),
                       heq = c(mu_X, mu_Y), localsolver = "MMA")


alabama::auglag(flb$solution, fn = qf$objective, gr = qf$gradient, 
                heq = function(mu){
                  return(mu - 1.0)
                }
)

