rm(list = ls())

sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/shrinkage-tpdm", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

library(tidyverse)
library(magrittr)
library(mev)
library(pbapply)
library(dplyover)
library(scales)
library(ggpubr)


# Run method ------------------------------------------------------------------------

set.seed(1)

k_fun <- function(n) floor(4 * sqrt(n))
kappa <- 0.75

# n \in seq(from = 10000, to = 100000, by = 10000)
# B = 5
# nrep = 5
# => ~ 30 mins
# => 100 rep in ~ 10 hours

params <- expand_grid(rho = c(0.15),
                      grid_size = 6,
                      n = seq(from = 10000, to = 200000, by = 10000),
                      B = 5,
                      rep = seq_len(100))

res <- pblapply(seq_len(nrow(params)), function(ind) {
  
  grid_size <- params$grid_size[ind]
  n <- params$n[ind]
  B <- params$B[ind]
  rho <- params$rho[ind]
  rep <- params$rep[ind]
  
  k <- k_fun(n)
  
  d <- grid_size^2 
  W <- create_rook_matrix(grid_size)
  W_eigenvalues <- eigen(W)$values
  rho_max <- min(1 / min(abs(W_eigenvalues)), 1 / max(abs(W_eigenvalues)))
  coord <- expand.grid(seq(from = 0.5, to = grid_size - 0.5, by = 1),
                       seq(from = 0.5, to = grid_size - 0.5, by = 1))
  dist_matrix <- dist(coord, diag = TRUE, upper = TRUE) %>% as.matrix()
  
  # compute true TPDM
  A <- solve(diag(d) - rho * W)
  D <- apply(A, 1, function(x) sqrt(sum(x^2))) %>% diag()
  A_tilde <- solve(D) %*% A
  Sigma <- A_tilde %*% t(A_tilde)

  # generate data
  X <- sim_SAR(n, A_tilde)
  Sigma_hat <- tpdm_estimate(X, k = k)
  
  # raw rho estimate
  rho_empirical_hat <- optim(par = 0, 
                             fn = function(x) frob_loss(Sigma_rho(x, W), Sigma_hat),
                             lower = 0, upper = 0.99 * rho_max, method = "L-BFGS-B")$par
  
  # oracle lambdas
  lambda_lw_oracle <- optim(par = 0,
                            fn = function(x) frob_loss(Sigma, shrink_cov(Sigma_hat, "lw", x)),
                            lower = 0, upper = 1, method = "L-BFGS-B")$par
  
  lambda_soft_oracle <- optim(par = 0,
                              fn = function(x) frob_loss(Sigma, shrink_cov(Sigma_hat, "soft", x)),
                              lower = 0, upper = 0.2, method = "L-BFGS-B")$par
  
  lambda_adlasso_oracle <- optim(par = 0,
                              fn = function(x) frob_loss(Sigma, shrink_cov(Sigma_hat, "adlasso", x, eta = 2)),
                              lower = 0, upper = 0.4, method = "L-BFGS-B")$par
  
  # oracle rhos
  Sigma_lw_oracle <- shrink_cov(Sigma_hat, "lw", lambda_lw_oracle)
  rho_lw_oracle <- optim(par = 0,
                         fn = function(x) frob_loss(Sigma_rho(x, W), Sigma_lw_oracle),
                         lower = 0, upper = 0.99 * rho_max, method = "L-BFGS-B")$par
  
  Sigma_soft_oracle <- shrink_cov(Sigma_hat, "soft", lambda_soft_oracle)
  rho_soft_oracle <- optim(par = 0,
                           fn = function(x) frob_loss(Sigma_rho(x, W), Sigma_soft_oracle),
                           lower = 0, upper = 0.99 * rho_max, method = "L-BFGS-B")$par
  
  Sigma_adlasso_oracle <- shrink_cov(Sigma_hat, "adlasso", lambda_adlasso_oracle, eta = 2)
  rho_adlasso_oracle <- optim(par = 0,
                              fn = function(x) frob_loss(Sigma_rho(x, W), shrink_cov(Sigma_hat, "adlasso", lambda_adlasso_oracle, eta = 2)),
                              lower = 0, upper = 0.99 * rho_max, method = "L-BFGS-B")$par
  
  # Fix et al. soft-threshold parameter
  lambda_soft_hat <- lambda_soft_fix(dist_matrix, Sigma_hat)
  Sigma_soft_hat <- shrink_cov(Sigma_hat, "soft", lambda_soft_hat)
  rho_soft_hat <- optim(par = 0,
                        fn = function(x) frob_loss(Sigma_rho(x, W), Sigma_soft_hat),
                        lower = 0, upper = 0.99 * rho_max, method = "L-BFGS-B")$par
  
  # LW ratio method
  lambda_lw_hat <- lw_ratio_method(X, k_fun, kappa, B = B) %>% pull(lambda_hat)
  Sigma_lw_hat <- shrink_cov(Sigma_hat, "lw", lambda_lw_hat)
  rho_lw_hat <- optim(par = 0,
                      fn = function(x) frob_loss(Sigma_rho(x, W), Sigma_lw_hat),
                      lower = 0, upper = 0.99 * rho_max, method = "L-BFGS-B")$par
  
  return(tibble("n" = n,
                "k" = k,
                "d" = d,
                "rho" = rho,
                "B" = B,
                "rep" = rep,
                "kappa" = kappa,
                "k_kappa" = k_fun(kappa * n),
                "lambda_lw_oracle" = lambda_lw_oracle, 
                "lambda_lw_hat" = lambda_lw_hat, 
                "lambda_soft_oracle" = lambda_soft_oracle, 
                "lambda_soft_hat" = lambda_soft_hat, 
                "lambda_adlasso_oracle" = lambda_adlasso_oracle, 
                "rho_empirical_hat" = rho_empirical_hat, 
                "rho_lw_oracle" = rho_lw_oracle, 
                "rho_lw_hat" = rho_lw_hat, 
                "rho_soft_oracle" = rho_soft_oracle, 
                "rho_soft_hat" = rho_soft_hat, 
                "rho_adlasso_oracle" = rho_adlasso_oracle, 
                "fl_empirical_hat" = frob_loss(Sigma, Sigma_hat), 
                "fl_lw_oracle" = frob_loss(Sigma, Sigma_lw_oracle), 
                "fl_lw_hat" = frob_loss(Sigma, Sigma_lw_hat), 
                "fl_soft_oracle" = frob_loss(Sigma, Sigma_soft_oracle), 
                "fl_soft_hat" = frob_loss(Sigma, Sigma_soft_hat), 
                "fl_adlasso_oracle" = frob_loss(Sigma, Sigma_adlasso_oracle)))
    
}) %>% bind_rows() 

res %>% saveRDS("scripts/shrinkage-tpdm/results/extremal-sar-ledoit-wolf.RDS")


p1 <- res %>% pivot_longer(
  cols = -c(n, k, d, rho, B, rep, kappa, k_kappa),
  names_to = c("param", "tpdm_estimator", "type"), 
  names_sep = "_",
  values_to = "value") %>%
  filter(param == "fl") %>%
  ggplot(aes(x = n, y = value, colour = tpdm_estimator, linetype = type)) +
  stat_summary(fun = mean, geom = "line") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)), labels = label_number(scale = 1/1000)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), transform = log10_trans()) +
  scale_colour_manual(values = c("darkgreen", "grey", "red", "blue"), labels = c(expression("Adaptive lasso," ~ eta == 1), "Empirical", "Ledoit-Wolf", "Soft")) +
  scale_linetype_manual(values = 1:2, labels = c("Estimated", "Oracle")) +
  labs(x = expression(n %*% 10^{-3}),
       y = expression(sum((tilde(sigma)[ij](n) - sigma[ij])^2, "", "")), 
       colour = "TPDM estimator",
       linetype = "Type") +
  theme_light()

p2 <- res %>% pivot_longer(
  cols = -c(n, k, d, rho, B, rep, kappa, k_kappa),
  names_to = c("param", "tpdm_estimator", "type"), 
  names_sep = "_",
  values_to = "value") %>%
  filter(param == "lambda") %>%
  ggplot(aes(x = n, y = value, colour = tpdm_estimator, linetype = type)) +
  stat_summary(fun = mean, geom = "line") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)), labels = label_number(scale = 1/1000)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
  scale_colour_manual(values = c("darkgreen", "red", "blue"), labels = c(expression("Adaptive lasso," ~ eta == 1), "Ledoit-Wolf", "Soft")) +
  scale_linetype_manual(values = 1:2, labels = c("Estimated", "Oracle")) +
  labs(x = expression(n %*% 10^{-3}),
       y = expression(hat(lambda)), 
       colour = "TPDM estimator",
       linetype = "Type") +
  theme_light()

p3 <- res %>% pivot_longer(
  cols = -c(n, k, d, rho, B, rep, kappa, k_kappa),
  names_to = c("param", "tpdm_estimator", "type"), 
  names_sep = "_",
  values_to = "value") %>%
  filter(param == "rho") %>%
  ggplot(aes(x = n, y = value - rho, colour = tpdm_estimator, linetype = type)) +
  stat_summary(fun = mean, geom = "line") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)), labels = label_number(scale = 1/1000)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_colour_manual(values = c("darkgreen", "grey", "red", "blue"), labels = c(expression("Adaptive lasso," ~ eta == 1), "Empirical", "Ledoit-Wolf", "Soft")) +
  scale_linetype_manual(values = 1:2, labels = c("Estimated", "Oracle")) +
  labs(x = expression(n %*% 10^{-3}),
       y = expression(hat(rho) - rho), 
       colour = "TPDM estimator",
       linetype = "Type") +
  theme_light()

ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = TRUE, legend = "top")

# ind <- which.min(res$frob_loss_Sigma_tilde)
# list(res$Sigma[[ind]], res$Sigma_hat[[ind]], res$Sigma_fix[[ind]], res$Sigma_tilde[[ind]]) %>%
#   lapply(function(S) {
#     rownames(S) <- 1:nrow(S)
#     colnames(S) <- 1:ncol(S)
#     return(S)
#   }) %>%
#   set_names(c("True TPDM", "Empirical TPDM", "Fix et al. (2021)", "Piecewise LW")) %>%
#   plot_tpdm()


# p3 <- res %>%
#   filter(n == 30000) %>%
#   pivot_longer(starts_with("lambda"), names_to = "lambda_type", values_to = "lambda_value") %>%
#   filter(lambda_type != "lambda_fix") %>%
#   mutate(lambda_type = fct_reorder(lambda_type, -lambda_value)) %>%
#   ggplot(aes(x = lambda_type, y = lambda_value)) +
#   geom_boxplot(outliers = FALSE) +
#   scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
#   scale_x_discrete(labels = c(expression(q == 1), expression(q == 2), expression(q == 3), expression(q == 4))) +
#   labs(y = expression(lambda),
#        x = expression(hat(lambda)^{(q)})) +
#   theme_light()
# 
# p4 <- lapply(c("Sigma", "Sigma_hat", "Sigma_tilde", "Sigma_fix"), function(S) {
#   s <- res %>%
#     filter(n == 30000) %>%
#     slice_min(frob_loss_Sigma_tilde) %>%
#     pull(S) %>%
#     Reduce("+", .) %>%
#     vecu()
#   return(data.frame("ij" = seq_len(length(s)), "sigma_ij" = s, "sigma_type" = S, "dist" = vecu(dist_matrix)))
#   }) %>% 
#   bind_rows() %>%
#   ggplot(aes(x = dist, y = sigma_ij, colour = sigma_type)) +
#   geom_jitter(alpha = 0.05, width = 0.1, shape = 20) +
#   scale_colour_manual(values = c("black", "blue", "grey", "red"), labels = c("True", "Fix et al. (2021)", "Empirical", "Piecewise Ledoit-Wolf")) +
#   scale_y_continuous(expand = expansion(mult = c(0.01, 0.03))) +
#   labs(x = "Distance",
#        y = "Dependence strength") +
#   theme_light()



res %>%
  filter(kappa == 0.75) %>%
  pivot_longer(c(lambda_soft, lambda_fix, lambda_lw), names_to = "lambda_type", values_to = "lambda_value") %>%
  group_by(n, kappa, rho, lambda_type) %>%
  summarise(lambda_value = median(lambda_value))
  ggplot(aes(x = n, y = lambda_value, colour = lambda_type)) +
  stat_summary(fun = median, geom = "line") +
  facet_grid(as.factor(rho) ~ ., scales = "free") +
  theme_light()



