library("Mapper")
library("testthat")
testthat::context("Testing mapper")
# skip_on_cran()

# 
# m <- MapperRef$new(runif(1000))
# m$use_cover(filter_values = cbind(runif(1000), runif(1000)), number_intervals=10, percent_overlap=60)
# m$construct_k_skeleton(k=2L)
# 
# ## Load chemical diabetes data set
# library(locfit)
# data("chemdiab", package = "locfit")
# 
# ## L-Infinity eccentricity 
# X <- tourr::rescale(as.matrix(chemdiab[, -6]))
# f_ecc <- matrix(apply(as.matrix(dist(X)), 1, max))
# 
# ## Make the mapper 
# chem_mapper <- function(f_x){
#   stopifnot(is.matrix(f_x))
#   H <- apply(X, 2, stats::bw.nrd0)
#   kde_res <- ks::kde(x = f_x, H = diag(H, nrow = ncol(f_x), ncol = ncol(f_x)), eval.points = f_x, binned = FALSE)
#   m <- Mapper::MapperRef$new(X = X)
#   m$use_cover(matrix(kde_res$estimate), typename = "fixed rectangular", number_intervals = 8L, percent_overlap = 35)
#   m$use_distance_measure(measure = "euclidean")
#   m$use_clustering_algorithm(cl = "single", num_bins = 10L)
#   m$compute_k_skeleton(k = 1L)
#   return(m)
# }
# 
# current_mapper <- chem_mapper(f_ecc)
# 
# ## Enable multiscale
# m$enable_multiscale()
# 
# wut <- m$update_mapper(percent_overlap = 30, stats = TRUE)
# 
# tourr::draw_tour_axes
# 
# ## Start a scene
# options(rgl.useNULL = TRUE)
# open3d(useNULL = TRUE)
# rgl_scene <- scene3d()
# rgl::rgl.close()
# 
# vc <- colnames(chemdiab[, -6])
# ui <- shinyUI(pageWithSidebar(
#   headerPanel = headerPanel("\"Visualization is insight; not pictures\""),
#   sidebarPanel = sidebarPanel(
#     registerSceneChange(),
#     sliderInput("alpha_index", label = "Alpha", min = 0, max = 1, value = 0, step = 1/25, 
#                 animate = animationOptions(interval=40*2, loop=FALSE)), 
#     selectInput("dim1_1", label = "From Dim 1", choices = vc, selected = vc[[1]], width = "30%"), 
#     selectInput("dim1_2", label = "From Dim 2", choices = vc, selected = vc[[2]], width = "30%"),
#     selectInput("dim1_3", label = "From Dim 3", choices = vc, selected = vc[[3]], width = "30%"),
#     selectInput("dim2_1", label = "To Dim 1", choices = vc, selected = vc[[1]], width = "30%"), 
#     selectInput("dim2_2", label = "To Dim 2", choices = vc, selected = vc[[2]], width = "30%"),
#     selectInput("dim2_3", label = "To Dim 3", choices = vc, selected = vc[[3]], width = "30%"),
#     actionButton("build_mapper", label = "Build Mapper")
#   ), 
#   mainPanel = mainPanel(
#     rgl::rglwidgetOutput("rgl_out", width = "100%", height = "50%"),
#     grapher::grapherOutput("grapher", width = "100%", height = "50%"),
#     style = "height: 100vh !important;"
#   )
# ))
# 
# server <- shinyServer(function(input, output, session) {
#   ## Initialize the scene
#   plot3d(rgl_scene)
#   
#   ## Set callback to appropriately cleanup the RGL session on exit
#   dev <- rgl.cur(); save <- options(rgl.inShiny = TRUE)
#   on.exit(options(save))
#   session$onSessionEnded(function() { rgl.set(dev); rgl.close() })
#   
#   ## Reactives
#   from_dims <- reactive({ c(input$dim1_1, input$dim1_2, input$dim1_3) })
#   to_dims <- reactive({ c(input$dim2_1, input$dim2_2, input$dim2_3) })
#   cpath <- reactive({
#     req(from_dims())
#     req(from_dims())
#     tourr::geodesic_path(as.matrix(X[, from_dims()]), as.matrix(X[, to_dims()]))
#   })
#   
#   ## Reactive values to stroe the current points
#   rgl_react <- reactiveValues()
#   rgl_react$data_pts <- NULL # rgl::points3d(tourr::rescale(path$interpolate(0)))
#   output$rgl_out <- rgl::renderRglwidget({ 
#     rgl_widget <<- rgl::rglwidget()
#     return(rgl_widget)
#   })
#   rgl::plot3d(NULL, xlim = c(0,1), ylim=c(0, 1), zlim=c(0, 1))
#   
#   ## Render mapper
#   output$grapher <- grapher::renderGrapher({
#     current_mapper$as_grapher()
#   })
#   
#   ## Observer to modify the projection along a geodesic
#   observeEvent(input$alpha_index, {
#     # req(input$alpha_index)
#     if (!is.null(rgl_react$data_pts)){ rgl::delFromSubscene3d(rgl_react$data_pts) }
#     new_pts <- rgl::points3d(tourr::rescale(cpath()$interpolate(input$alpha_index)), color = palette()[chemdiab$cc])
#     session$sendCustomMessage(
#       "sceneChange",
#       sceneChange("rgl_out", delete = rgl_react$data_pts, add = new_pts)
#     )
#     rgl_react$data_pts <- new_pts
#   })
#   
#   observeEvent(input$build_mapper, {
#     req(input$build_mapper)
#     # f_x <- tourr::rescale(cpath()$interpolate(input$alpha_index))
#     # current_mapper <<- chem_mapper(f_x)
#     # output$grapher <- grapher::renderGrapher({ current_mapper$as_grapher() })
#   })
# })
# 
# shiny::shinyApp(ui, server)
# 
# 
# ggvis::ggvis()
# plot(path$interpolate(0), col = chemdiab$cc)
# 
# normalize <- function(x){ (x - min(x))/(diff(range(x))) }
# it <- 1:100
# for (i in normalize(it)){
#   rgl::points3d(path$interpolate(i))
# }
# rgl::delFromSubscene3d()
# rgl::plot3d(path$interpolate(0))
# rgl::clear3d(type = "shapes")
# 
# rgl::rgl.dev.list()
# rgl::points3d()
# 
# x <- do.call(rbind, list(replicate(2, rnorm(25, mean = 0, sd = 1/3)), 
#                          replicate(2, rnorm(25, mean = 1.5, sd = 1/3)), 
#                          replicate(2, rnorm(25, mean = 3, sd = 1/3))))
# plot(x, col = c(rep_len(1, 25), rep_len(2, 25), rep_len(3, 25)))
# hcl <- hclust(dist(x), method = "single")
# diam <- max(dist(x))
# Mapper::cutoff_first_bin(hcl, diam = diam, num_bins = 10L)
# 
