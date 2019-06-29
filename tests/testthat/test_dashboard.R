# ## Test dashboard 
# library("Mapper")
# data("noisy_circle")
# left_pt <- noisy_circle[which.min(noisy_circle[, 1]),]
# f_x <- matrix(apply(noisy_circle, 1, function(pt) (pt - left_pt)[1]))
# m <- mapper(X = noisy_circle, filter_values = f_x, 
#             cover_params = list(typename="restrained rectangular",number_intervals=10L, percent_overlap=50),
#             measure = "euclidean", 
#             cluster_params = list(cl="single", num_bins=10L), return_reference = TRUE)
# dashboard(m, as.data.frame(noisy_circle))