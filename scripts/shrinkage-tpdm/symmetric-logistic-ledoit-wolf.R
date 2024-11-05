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

set.seed(1)

k_fun <- function(n) floor(4 * sqrt(n))
sl_tpdm <- Vectorize(sl_tpdm)

kappa_vals <- c(0.5, 0.75, 0.9)

params <- expand_grid(gamma = c(0.5, 0.8, 0.9),
            d = 3,
            n = c(seq(from = 5000, to = 15000, by = 5000), 
                  seq(from = 20000, to = 50000, by = 10000),
                  seq(from = 60000, to = 200000, by = 20000)),
            B = 10,
            rep = seq_len(100))

res <- pblapply(seq_len(nrow(params)), function(ind) {
  
  gamma <- params$gamma[ind]
  d <- params$d[ind]
  n <- params$n[ind]
  kappa <- params$kappa[ind]
  B <- params$B[ind]
  rep <- params$rep[ind]
  
  # compute true TPDM
  Sigma <- matrix(sl_tpdm(gamma), d, d)
  
  # generate data
  X <- rmev(n = n, d = d, par = gamma, model = "log") %>% sqrt()
  
  # compute empirical TPDM
  k <- k_fun(n)
  Sigma_hat <- tpdm_estimate(X, k)
  
  # bootstrap data of size n and compute empirical TPDM
  Sigma_hat_boot <- lapply(seq_len(B), function(b) {
    X[sample(1:n, size = n, replace = TRUE), ] %>% tpdm_estimate(k = k)
  })
  
  # compute average bootstrapped empirical TPDM
  Sigma_hat_boot_mean <- Reduce("+", Sigma_hat_boot) / length(Sigma_hat_boot)
  
  # compute lambda star
  lambda_star <- optim(par = 0, 
                       fn = function(x) frob_loss(Sigma, shrink_cov(Sigma_hat, "lw", x)),
                       lower = 0, upper = 1, method = "L-BFGS-B")$par
  
  # run the rest of the method for each kappa
  tmp <- lapply(kappa_vals, function(kappa) {
    
    k_kappa <- k_fun(floor(kappa * n))
    
    # bootstrap data of size kappa * n and compute empirical TPDM
    Sigma_hat_small_boot <- lapply(seq_len(B), function(b) {
      X[sample(1:n, size = floor(kappa * n), replace = TRUE), ] %>% tpdm_estimate(k = k_kappa)
    })
    
    # compute average bootstrapped empirical TPDM
    Sigma_hat_small_boot_mean <- Reduce("+", Sigma_hat_small_boot) / length(Sigma_hat_small_boot)
    
    # estimate ratio zeta
    zeta_hat_num <- k * sum(vecu(Sigma_hat_boot_mean)^2)
    zeta_hat_den <- k_kappa * sum(vecu(Sigma_hat_small_boot_mean)^2)
    zeta_hat <- zeta_hat_num / zeta_hat_den
    
    # compute lambda_hat
    lambda_hat <- optim(par = 0, 
                        fn = function(x) frob_loss(shrink_cov(Sigma_hat_boot_mean, "lw", x),
                                                   shrink_cov(Sigma_hat_small_boot_mean, "lw", x * zeta_hat)),
                        lower = 0, upper = min(1, 1 / zeta_hat), method = "L-BFGS-B")$par
    
    return(data.frame("n" = n,
                      "k" = k,
                      "d" = d,
                      "gamma" = gamma,
                      "sigma" = sl_tpdm(gamma),
                      "B" = B,
                      "rep" = rep,
                      "sigma_hat_ij_mean" = mean(vecu(Sigma_hat)),
                      "E_sigma_hat_norm" = sum(vecu(Sigma_hat_boot_mean)^2),
                      "lambda_star" = lambda_star,
                      "kappa" = kappa,
                      "k_kappa" = k_kappa,
                      "E_sigma_hat_small_norm" = sum(vecu(Sigma_hat_small_boot_mean)^2),
                      "zeta_hat" = zeta_hat,
                      "lambda_hat" = lambda_hat))
    
    }) %>% bind_rows()
  return(tmp) 
  }) %>% 
  bind_rows() %>%
  saveRDS("scripts/shrinkage-tpdm/results/symmetric-logistic-ledoit-wolf.RDS")

res <- readRDS("scripts/shrinkage-tpdm/results/symmetric-logistic-ledoit-wolf.RDS")

# lambda_true and lambda_hat against n
p1 <- ggplot(res, aes(x = n, y = lambda_hat, colour = as.factor(gamma), linetype = as.factor(kappa))) +
  stat_summary(aes(y = lambda_star, colour = as.factor(gamma)), fun = median, geom = "line", linetype = 1) +
  stat_summary(fun = median, geom = "line") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)), labels = label_number(scale = 1/1000)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.04))) +
  scale_colour_manual(values = c("blue", "red", "darkgreen")) +
  scale_linetype_manual(values = c(2, 3, 6)) +
  labs(x = expression(n %*% 10^{-3}),
       y = expression(hat(lambda)(n)),
       colour = expression(gamma),
       linetype = expression(kappa)) +
  theme_light()

# mean entry of Sigma_hat against n
p2 <- ggplot(res, aes(x = n, y = E_sigma_hat_norm - (choose(d, 2) * sl_tpdm(gamma)^2), colour = as.factor(gamma))) +
  stat_summary(fun = mean, geom = "line") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)), labels = label_number(scale = 1/1000)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
  scale_colour_manual(values = c("blue", "red", "darkgreen")) +
  scale_linetype_manual(values = c(2, 3, 6)) +
  labs(x = expression(n %*% 10^{-3}),
       y = expression(sum((hat(sigma)[ij](n) - sigma[ij](n))^2, "", "")), 
       colour = expression(gamma)) +
  theme_light()

# zeta_hat against n
p3 <- ggplot(res, aes(x = n, y = zeta_hat, colour = as.factor(gamma), linetype = as.factor(kappa))) +
  stat_summary(fun = median, geom = "line") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)), labels = label_number(scale = 1/1000)) +
  scale_linetype_manual(values = c(2, 3, 6)) +
  scale_colour_manual(values = c("blue", "red", "darkgreen")) +
  labs(x = expression(n %*% 10^{-3}),
       y = expression(hat(zeta)(n,kappa)),
       colour = expression(gamma),
       linetype = expression(kappa)) +
  theme_light()

# k * norm(Sigma_hat) * lambda_true against n
p4 <- ggplot(res, aes(x = n, y = k * E_sigma_hat_norm * lambda_star, colour = as.factor(gamma))) +
  stat_summary(fun = median, geom = "line") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)), labels = label_number(scale = 1/1000)) +
  scale_colour_manual(values = c("blue", "red", "darkgreen")) +
  scale_linetype_manual(values = c(2, 3, 6)) +
  labs(x = expression(n %*% 10^{-3}),
       y = expression(k(n) ~ lambda(n) ~ sum(hat(sigma)[ij](n)^2, "", "")),
       colour = expression(gamma)) +
  theme_light()

ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, legend = "right", common.legend = TRUE)



  