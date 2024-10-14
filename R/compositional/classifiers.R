# 
# # Sim study -------------------------------------------------------------------------
# 
# library(tidyverse)
# library(magrittr)
# library(pbapply)
# library(caret)
# library(Compositional)
# library(CompositionalML)
# library(scales)
# library(ggh4x)
# 
# a_vals <- seq(-1, 1, by = 0.1)
# 
# params <- expand_grid(
#   model = c("log"),
#   param0 = 2, 
#   param1 = c(1.1, 1.5, seq(from = 2, to = 6, by = 1)), 
#   d = 5,
#   k_frac = c(0.05, 0.15),
#   classifier = c("knn")) %>%
#   mutate(param_index = row_number())
# 
# res <- pblapply(seq_len(nrow(params)), function(i) {
#   
#   # get parameter values
#   d <- params$d[i]
#   model <- params$model[i] 
#   param0 <- params$param0[i] 
#   param1 <- params$param1[i] 
#   k_frac <- params$k_frac[i]
#   
#   nreps <- 25
#   k_vals <- 1:25
#   
#   # perform simulations and run tests
#   out <- replicate(n = nreps, expr = {
#     
#     data <- sim_X_classify(n = 10000, d, param0, param1, model)
#     data_partition <- extreme_data_partition(data = data, train_frac = 0.75, k_frac = k_frac)
#     
#     xnew <- select(data_partition$test_ext_data, starts_with("theta"))
#     x <- select(data_partition$train_ext_data, starts_with("theta"))
#     ina <- pull(data_partition$train_ext_data, class)
#     inanew <- pull(data_partition$test_ext_data, class)
#     
#     knn_tune <- alfaknn.tune(x = x, ina = ina, a = a_vals, k = k_vals)
#     
#     optimal_alfa <- knn_tune$best_a
#     optimal_k <- knn_tune$best_k
#     aitchison_k <-  k_vals[which.max(knn_tune$ela["alpha=0", ])]
#     euclidean_k <-  k_vals[which.max(knn_tune$ela["alpha=1", ])]
#     
#     optimal_ce <- alfa.knn(xnew = xnew, x = x, ina = ina, a = optimal_alfa, k = optimal_k) %>%
#       ModelMetrics::ce(inanew, as.factor(.))
#     
#     aitchison_ce <- alfa.knn(xnew = xnew, x = x, ina = ina, a = 0, k = aitchison_k) %>%
#       ModelMetrics::ce(inanew, as.factor(.))
#     
#     euclidean_ce <- alfa.knn(xnew = xnew, x = x, ina = ina, a = 1, k = euclidean_k) %>%
#       ModelMetrics::ce(inanew, as.factor(.))
#     
#     return(list("optimal_alfa" = optimal_alfa,
#                 "optimal_k" = optimal_k,
#                 "aitchison_k" = aitchison_k,
#                 "euclidean_k" = euclidean_k,
#                 "optimal_ce" = optimal_ce,
#                 "aitchison_ce" = aitchison_ce,
#                 "euclidean_ce" = euclidean_ce))
#   })
#   
# }) %>%
#   lapply(t) %>%
#   lapply(function(M) as.data.frame(M) %>% mutate(rep = seq_len(nrow(M)))) %>%
#   bind_rows(.id = "param_index") %>%
#   unnest(c(optimal_alfa, optimal_k, aitchison_k, euclidean_k, optimal_ce, aitchison_ce, euclidean_ce)) %>%
#   mutate(param_index = as.integer(param_index)) %>%
#   full_join(params, ., by = "param_index") %>%
#   pivot_longer(cols = optimal_alfa:euclidean_ce,
#                names_to = c("method", ".value"),
#                names_pattern = "(.+)_(.+)")
# 
# 
# 
# data <- sim_X_classify(n = 10000, d = 5, param0 = 2, param1 = 1.1, model = "log")
# data_partition <- extreme_data_partition(data = data, train_frac = 0.75, k_frac = 0.05)
# xnew <- select(data_partition$test_ext_data, starts_with("theta"))
# x <- select(data_partition$train_ext_data, starts_with("theta"))
# ina <- pull(data_partition$train_ext_data, class)
# inanew <- pull(data_partition$test_ext_data, class)
# knn_tune_log <- alfaknn.tune(x = x, ina = ina, a = a_vals, k = 1:25, graph = FALSE)
# 
# 
# 
# res %>%
#   mutate(param1 = as.factor(param1)) %>%
# ggplot(aes(x = method, y = ce, fill = method)) +
#   geom_boxplot() +
#   facet_grid(as.factor(k_frac) ~ param1, labeller = label_parsed) +
#   scale_fill_manual(values = c("red", "blue", "darkgrey"),
#                     labels = c(expression(bold(A)*"itchison," ~ alpha == 0),
#                                expression(bold(E)*"uclidean," ~ alpha == 1),
#                                expression(bold(O)*"ptimal" ~ alpha))) +
#   scale_x_discrete(labels = c("A", "E", "O")) +
#   scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)), breaks = breaks_pretty(n = 4)) +
#   geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
#   xlab("Classifier") +
#   ylab("Classification error") +
#   labs(fill = expression(k * "-NN" * (alpha) ~ "classifier")) +
#   theme_light()
# 
# res %>%
#   filter(method == "optimal", k_frac == 0.15) %>%
#   mutate(param1 = as.factor(param1)) %>%
#   ggplot(aes(x = param1, y = alfa)) +
#   geom_boxplot(fill = "lightgrey") +
#   xlab(expression(vartheta[1])) +
#   ylab(expression("Optimal" ~ alpha)) +
#   scale_y_continuous(limits = c(-1, 1), expand = c(0, 0.01)) +
#   geom_hline(yintercept = 1, linetype = "dashed", colour = "blue") +
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
#   theme_light()
# 
# p3 <- ggplot(filter(res, distance == "alfa"), aes(x = k, fill = as.factor(d))) +
#   geom_histogram(position = "identity", alpha = 0.8) + 
#   facet_grid(~ as.factor(param1))
# 
# p1
# p2 / p3
# 
# 
# # NFL Combine -----------------------------------------------------------------------
# 
# set.seed(123)
# 
# data <- load_nfl_combine() %>% select(FortyYard:class)
# data_partition <- extreme_data_partition(data = data, train_frac = 0.75, k_frac = 0.1)
# 
# xnew <- select(data_partition$test_ext_data, starts_with("theta"))
# x <- select(data_partition$train_ext_data, starts_with("theta"))
# ina <- pull(data_partition$train_ext_data, class)
# inanew <- pull(data_partition$test_ext_data, class)
# 
# k_vals_nfl <- 1:20
# 
# knn_tune_nfl <- alfaknn.tune(x = x, ina = ina, a = a_vals, k = k_vals_nfl, graph = FALSE)
# 
# optimal_alfa <- knn_tune_nfl$best_a
# optimal_k <- knn_tune_nfl$best_k
# aitchison_k <-  k_vals_nfl[which.max(knn_tune_nfl$ela["alpha=0", ])]
# euclidean_k <-  k_vals_nfl[which.max(knn_tune_nfl$ela["alpha=1", ])]
# 
# optimal_ce <- alfa.knn(xnew = xnew, x = x, ina = ina, a = optimal_alfa, k = optimal_k) %>%
#   ModelMetrics::ce(inanew, as.factor(.))
# 
# aitchison_ce <- alfa.knn(xnew = xnew, x = x, ina = ina, a = 0, k = aitchison_k) %>%
#   ModelMetrics::ce(inanew, as.factor(.))
# 
# euclidean_ce <- alfa.knn(xnew = xnew, x = x, ina = ina, a = 1, k = euclidean_k) %>%
#   ModelMetrics::ce(inanew, as.factor(.))
# 
# list("optimal_alfa" = optimal_alfa,
#      "optimal_k" = optimal_k,
#      "aitchison_k" = aitchison_k,
#      "euclidean_k" = euclidean_k,
#      "optimal_ce" = optimal_ce,
#      "aitchison_ce" = aitchison_ce,
#      "euclidean_ce" = euclidean_ce)
# 
# 
# # Ecoli -----------------------------------------------------------------------
# 
# set.seed(1)
# 
# data <- load_ecoli() %>% select(-sequence)
# data_partition <- extreme_data_partition(data = data, train_frac = 0.7, k_frac = 0.20)
# 
# xnew <- select(data_partition$test_ext_data, starts_with("theta"))
# x <- select(data_partition$train_ext_data, starts_with("theta"))
# ina <- pull(data_partition$train_ext_data, class)
# inanew <- pull(data_partition$test_ext_data, class)
# 
# k_vals_ecoli <- 1:15
# 
# knn_tune_ecoli <- alfaknn.tune(x = x, ina = ina, a = a_vals, k = k_vals_ecoli, graph = FALSE)
# 
# optimal_alfa <- knn_tune_ecoli$best_a
# optimal_k <- knn_tune_ecoli$best_k
# aitchison_k <-  k_vals_ecoli[which.max(knn_tune_ecoli$ela["alpha=0", ])]
# euclidean_k <-  k_vals_ecoli[which.max(knn_tune_ecoli$ela["alpha=1", ])]
# 
# optimal_ce <- alfa.knn(xnew = xnew, x = x, ina = ina, a = optimal_alfa, k = optimal_k) %>%
#   ModelMetrics::ce(inanew, as.factor(.))
# 
# aitchison_ce <- alfa.knn(xnew = xnew, x = x, ina = ina, a = 0, k = aitchison_k) %>%
#   ModelMetrics::ce(inanew, as.factor(.))
# 
# euclidean_ce <- alfa.knn(xnew = xnew, x = x, ina = ina, a = 1, k = euclidean_k) %>%
#   ModelMetrics::ce(inanew, as.factor(.))
# 
# list("optimal_alfa" = optimal_alfa,
#      "optimal_k" = optimal_k,
#      "aitchison_k" = aitchison_k,
#      "euclidean_k" = euclidean_k,
#      "optimal_ce" = optimal_ce,
#      "aitchison_ce" = aitchison_ce,
#      "euclidean_ce" = euclidean_ce)
# 
# 
# # Plots -----------------------------------------------------------------------------
# 
# lattice.options(
#   axis.padding = list(factor = 0),
#   bottom.padding = list(factor = 0),
#   top.padding = list(factor = 0),
#   key.sub.padding = list(factor = 0),
#   key.axis.padding = list(factor = 0),
#   main.key.padding = list(factor = 0),
#   layout.heights = list(bottom.padding = list(x=0), top.padding = list(x=0)),
#   layout.widths = list(left.padding = list(x = 0), right.padding = list(x = 0))
# )
# 
# p1 <- knn_tune_log$ela %>%
#   set_colnames(1:25) %>%
#   set_rownames(a_vals) %>% 
#   levelplot(xlab = expression(alpha), 
#             ylab = expression(k),
#             col.regions = rev(heat.colors(100)),
#             main = "Simulated",
#             region.type = "contour",
#             scales = list(x = list(labels = c(-1, -0.5 ,0 , 0.5, 1), at = which(a_vals %in% c(-1, -0.5 ,0 , 0.5, 1))),
#                           y = list(labels = seq(from = 1, to = max(1:25), by = 2), at = which(1:25 %in% seq(from = 1, to = max(1:25), by = 2)))),
#             colorkey = list(space = "right")) %>%
#   update(aspect = 0.7)
# 
# p2 <- knn_tune_nfl$ela %>%
#   set_colnames(k_vals_nfl) %>%
#   set_rownames(a_vals) %>% 
#   levelplot(xlab = expression(alpha), 
#             ylab = expression(k),
#             col.regions = rev(heat.colors(100)),
#             main = "NFL Combine",
#             region.type = "contour",
#             scales = list(x = list(labels = c(-1, -0.5 ,0 , 0.5, 1), at = which(a_vals %in% c(-1, -0.5 ,0 , 0.5, 1))),
#                           y = list(labels = seq(from = 1, to = max(k_vals_nfl), by = 2), at = which(k_vals_nfl %in% seq(from = 1, to = max(1:25), by = 2)))),
#                           colorkey = list(space = "right")) %>%
#               update(aspect = 0.7)
# 
# p3 <- knn_tune_ecoli$ela %>%
#   set_colnames(k_vals_ecoli) %>%
#   set_rownames(a_vals) %>% 
#   levelplot(xlab = expression(alpha), 
#             ylab = expression(k),
#             col.regions = rev(heat.colors(100)),
#             main = "Ecoli",
#             region.type = "contour",
#             scales = list(x = list(labels = c(-1, -0.5 ,0 , 0.5, 1), at = which(a_vals %in% c(-1, -0.5 ,0 , 0.5, 1))),
#                           y = list(labels = seq(from = 1, to = max(k_vals_ecoli), by = 2), at = which(k_vals_ecoli %in% seq(from = 1, to = max(k_vals_ecoli), by = 2)))),
#             colorkey = list(space = "right")) %>%
#   update(aspect = 0.7)
# 
# 
# p1 <- load_nfl_combine() %>%
#   select(c(FortyYard, ThreeCone, Shuttle, class)) %>%
#   mutate(R = rowSums(select(., -class))) %>%
#   slice_max(R, prop = 0.025) %>%
#   ggtern(aes(FortyYard, ThreeCone, Shuttle, colour = class, size = log2(R))) +
#   geom_mask() +
#   geom_point(alpha = 0.5) +
#   xlab("") +
#   ylab("") +
#   Llab("FortyYard") +
#   Tlab("ThreeCone") +
#   Rlab("Shuttle") +
#   labs(colour = "Position", size = "Radius") +
#   scale_colour_manual(labels = c("On The Line", "Off The Line"), values = c("blue", "red")) +
#   theme_hidegrid() +
#   theme_hidelabels() +
#   theme(legend.position = "right")
# 
# p2 <- load_nfl_combine() %>%
#   select(c(Bench, Vertical, BroadJump, class)) %>%
#   mutate(R = rowSums(select(., -class))) %>%
#   slice_max(R, prop = 0.025) %>%
#   ggtern(aes(Bench, Vertical, BroadJump, colour = class, size = log2(R))) +
#   geom_mask() +
#   geom_point(alpha = 0.5) +
#   xlab("") +
#   ylab("") +
#   Llab("Bench") +
#   Tlab("Vertical") +
#   Rlab("BroadJump") +
#   labs(colour = "Position", size = "Radius") +
#   scale_colour_manual(labels = c("On The Line", "Off The Line"), values = c("blue", "red")) +
#   theme_hidegrid() +
#   theme_hidelabels() +
#   theme(legend.position = "right")
# 
# ggarrange(print(p1), print(p2), ncol = 2, common.legend = TRUE, legend = "right")
# 
# 
# 
# data <- sim_X_classify(n = 15000, d = 10, param0 = 500, param1 = 500, model = "maxlin")
# data_partition <- extreme_data_partition(data = data, train_frac = 0.75, k_frac = 0.1)
# 
# xnew <- select(data_partition$test_ext_data, starts_with("theta"))
# x <- select(data_partition$train_ext_data, starts_with("theta"))
# ina <- pull(data_partition$train_ext_data, class)
# inanew <- pull(data_partition$test_ext_data, class)
# 
# k_vals <- 2:15
# knn_tune <- alfaknn.tune(x = x, ina = ina, a = a_vals, k = k_vals, graph = TRUE)
# alfa.knn(xnew = xnew, x = x, ina = ina, a = 0.9, k = 7) %>%
#   ModelMetrics::ce(inanew, as.factor(.))
