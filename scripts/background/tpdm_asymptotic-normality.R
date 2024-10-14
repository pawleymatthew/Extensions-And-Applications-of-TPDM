
# sl_h1 <- function(x, gamma) {
#   term1 <- (1 - gamma) / (2 * gamma)
#   term2 <- (x * (1 - x))^(1 / gamma - 2)
#   term3 <- (x^(1 / gamma) + (1 - x)^(1 / gamma))^(gamma - 2)
#   h1 <- term1 * term2 * term3
#   return(h1)
# }
# 
# sigma <- function(h1, par) {
#   integrate(f = function(x) 2 * sqrt(x * (1 - x)) * h1(x, par), lower = 0, upper = 1) %>%
#     extract2("value")
# }
# 
# nu_squared <- function(h1, par) {
#   term1 <- integrate(f = function(x) 2^2 * x * (1 - x) * h1(x, par),lower = 0, upper = 1)%>%
#     extract2("value")
#   term2 <- sigma(h1, par)^2
#   return(term1 - term2)
# }
# 
# sigma_ci <- function(h1, par, k, ci_level = 0.05) {
#   z <- qnorm(1 - ci_level / 2)
#   mean <- sigma(h1, par)
#   sd <- sqrt(nu_squared(h1, par) / k)
#   lower <- mean - z * sd
#   upper <- mean + z * sd
#   return(c(lower, upper))
# }
# 
# 
# gamma_vals <- c(0.6, 0.8)
# nrep <- 250
# n <- 10^6
# k <- sqrt(n)
# i <- 1
# j <- 2
# l <- 2
# m <- 3
# 
# res <- lapply(gamma_vals, function(gamma) {
#   sigma_hat_vals <- pbsapply(seq_len(nrep), function(rep) {
#     X <- rmev(n = n, d = 3, par = 1 / gamma, model = "log") %>% sqrt()
#     Sigma_hat <- tpdm_estimate(X = X, k = k)
#     sigma_hat_ij <- Sigma_hat[i,j]
#     sigma_hat_lm <- Sigma_hat[l,m]
#     return(c(sigma_hat_ij, sigma_hat_lm))
#   })
#   data.frame("gamma" = gamma,
#              "rep" = seq_len(nrep),
#              "sigma_hat_ij" = sigma_hat_vals[1, ],
#              "sigma_hat_lm" = sigma_hat_vals[2, ])
# }) %>%
#   bind_rows() %>% 
#   rowwise() %>%
#   mutate(sigma_true = sigma(sl_h1, gamma),
#          asy_sd = sqrt(nu_squared(sl_h1, par = gamma) / k))
# 
# ggplot(res, aes(x = sigma_hat_ij, y = sigma_hat_lm, colour = as.factor(gamma), fill = as.factor(gamma))) +
#   geom_point() +
#   geom_xsidehistogram(aes(y = after_stat(density)), alpha = 0.3, bins = 50, colour = "black", linewidth = 0.05) +
#   geom_ysidehistogram(aes(x = after_stat(density)), alpha = 0.3, bins = 50, colour = "black", linewidth = 0.05) +
#   stat_xsidefunction(fun = dnorm, args = list(mean = sigma(sl_h1, par = gamma_vals[1]), 
#                                               sd = sqrt(nu_squared(sl_h1, par = gamma_vals[1]) / k)),
#                      colour = "red", xlim = c(0.65, 0.85)) +
#   stat_xsidefunction(fun = dnorm, args = list(mean = sigma(sl_h1, par = gamma_vals[2]), 
#                                               sd = sqrt(nu_squared(sl_h1, par = gamma_vals[2]) / k)),
#                      colour = "blue", xlim = c(0.38, 0.65)) +
#   stat_ysidefunction(fun = dnorm, args = list(mean = sigma(sl_h1, par = gamma_vals[1]), 
#                                               sd = sqrt(nu_squared(sl_h1, par = gamma_vals[1]) / k)),
#                      colour = "red", ylim = c(0.65, 0.85)) +
#   stat_ysidefunction(fun = dnorm, args = list(mean = sigma(sl_h1, par = gamma_vals[2]), 
#                                               sd = sqrt(nu_squared(sl_h1, par = gamma_vals[2]) / k)),
#                      colour = "blue", ylim = c(0.38, 0.65)) +
#   geom_hline(aes(yintercept = sigma_true, colour = as.factor(gamma)), linetype = "dashed") +
#   geom_vline(aes(xintercept = sigma_true, colour = as.factor(gamma)), linetype = "dashed") +
#   scale_y_continuous(expand = c(0.01, 0.01)) +
#   scale_x_continuous(expand = c(0.01, 0.01)) +
#   scale_fill_manual(values = c("red", "blue")) +
#   scale_colour_manual(values = c("red", "blue")) +
#   theme_light() +
#   labs(colour = expression(gamma),
#        fill = expression(gamma),
#        x = expression(hat(sigma)[ij]),
#        y = expression(hat(sigma)[lm])) +
#   theme(legend.position = "top",
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         ggside.panel.scale.x = 0.4,
#         ggside.panel.scale.y = 0.25)
#   
# 
# gamma <- 0.6
# sigma(sl_h1, par = gamma)
# nu_squared(sl_h1, par = gamma)
# sigma(sl_h1, par = gamma) - qnorm(1 - 0.05 / 2) * sqrt(nu_squared(sl_h1, par = gamma) / k)
# sigma(sl_h1, par = gamma) + qnorm(1 - 0.05 / 2) * sqrt(nu_squared(sl_h1, par = gamma) / k)
# res %>% group_by(gamma) %>% summarise(cor_sigma = cor(sigma_hat_ij, sigma_hat_lm))
# 
# res <- pbsapply(seq_len(nrep), function(i) {
#   X <- rmev(n = n, d = 3, par = 1 / gamma, model = "log") %>% sqrt()
#   sigma <- tpdm_estimate(X = X, k = k) %>% as.numeric() %>% extract(2)
#   return(sigma)
# })
# 
# data <- data.frame("rep" = seq_len(nrep), "sigma_hat" = res)
# 
# ggplot(data, aes(x = sigma_hat)) +
#   geom_histogram(aes(y = after_stat(density)), alpha = 0.3, colour = "darkgrey") +
#   stat_function(fun = dnorm, args = list(mean = sigma(sl_h1, par = gamma), 
#                                          sd = sqrt(nu_squared(sl_h1, par = gamma) / k))) +
#   geom_vline(xintercept = sigma(sl_h1, par = gamma), colour = "red", linetype = "dashed") +
#   scale_x_continuous(expand = expansion(mult = c(0, 0))) +
#   scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
#   theme_light() +
#   labs(colour = "Measure",
#        x = expression(hat(sigma)[ij]),
#        y = "Density") +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank())
# 
# hist(res, breaks = 50, prob = TRUE)
# abline(v = sigma(sl_h1, par = 0.5), col = "red")
# abline(v = lower, col = "blue")
# abline(v = upper, col = "blue")
# curve(dnorm(x, 
#             mean = sigma(sl_h1, par = gamma), 
#             sd = sqrt(nu_squared(sl_h1, par = gamma) / k)),
#       add = TRUE)
# 
# mean(res > lower & res < upper)
# 
# res <- pbsapply(seq_len(200), function(i) {
#   X <- rmev(n = n, d = 3, par = 1 / gamma, model = "log") %>% sqrt()
#   sigma <- tpdm_estimate(X = X, k = k)
#   return(c(sigma[1, 2], sigma[2, 3]))
# }) %>%
#   t()
# 
# cor(res[, 1], res[, 2])


# Covariance of TPDM entries using max-linear model ---------------------------------

library(tidyverse)
library(magrittr)
library(mev)
library(pbapply)

sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/background", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

set.seed(3)
alpha <- 1
q <- 12
d <- 4
A <- matrix(rexp(d * q, rate = 0.01), nrow = d, ncol = q) 
A <- A / rowSums(A)
rownames(A) <- paste0("X[", 1:d, "]")
colnames(A) <- paste0("bold(a)[", 1:q, "]")

Sigma <- A^(alpha / 2) %*% t(A^(alpha / 2))
rownames(Sigma) <- paste0("X[", 1:d, "]")
colnames(Sigma) <- rownames(Sigma)

V <- asy_V(A)

n <- 10^4
k <- sqrt(n)
res <- pbsapply(seq_len(1000), function(rep) {
  X <- SpatialExtremes::rmaxlin(n = n, dsgn.mat = A)
  R <- apply(X, 1, pracma::Norm, p = alpha)
  ext_ind <- which(R > quantile(R, probs = 1 - k/n))
  Theta_ext <- X[ext_ind, ] / R[ext_ind]
  W <- Theta_ext^(alpha / 2)
  Sigma_hat <- d * t(W) %*% W / k
  Sigma_hat <- t(Sigma_hat)[lower.tri(t(Sigma_hat))] # upper half vectorisation
  return(Sigma_hat)
}) %>%
  t() %>%
  set_colnames(colnames(V))

cov(res) # empirical V
asy_V(A) / k # theoretical V

saveRDS(list("A" = A,
             "Sigma" = Sigma,
             "V" = V, 
             "n" = n,
             "k" = k,
             "emp_V" = res), 
        file.path("scripts", "background", "results", "max-linear-asymptotic-tpdm.RDS"))


