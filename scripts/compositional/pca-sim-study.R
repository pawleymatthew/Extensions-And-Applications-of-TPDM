rm(list = ls())

options(encoding = "utf-8")

sapply(list.files(path = "R/compositional", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

library(tidyverse)
library(magrittr)
library(SpatialExtremes)
library(compositions)
library(Ternary)
library(pracma)
library(pbapply)
library(ggh4x)
library(scales)
library(ggpubr)
library(mev)
library(Matrix)


# Max-linear exact compositional example --------------------------------------------


# Generate A matrix -----------------------------------------------------------------

set.seed(1)

tau_vals <- runif(n = 50, -4, 4)
pc_dir <- c(0.12, 0.58, 0.3) # c(0.2, 0.5, 0.3) very curved; c(0.25, 0.5, 0.25) very straight

A <- matrix(NA, nrow = length(pc_dir), ncol = length(tau_vals))
for (j in 1:ncol(A)) {
  A[, j] <- power.acomp(pc_dir, tau_vals[j])
}
A <- A / rowSums(A)

saveRDS(A, file = file.path("scripts", "compositional", "results", "sim-pca-az-A-matrix.RDS"))

## Run simulations -------------------------------------------------------------------

set.seed(1)

expand_grid(
  n_train = c(5000, 50000),
  n_test = 5000,
  k_frac = c(0.01, 0.05),
  maxlin_model = c("maxlin", "cooley"),
  pca_method = c("ds", "coda"),
  n_pcs = 1:2,
  rep = 1:50
) %>%
  mutate(param_index = row_number()) %>%
  rowwise() %>%
  mutate(pca = list(extreme_pca(X_train = rmaxlinearfactor(n = n_train, A = A, type = maxlin_model), 
                                X_test = rmaxlinearfactor(n = n_test, A = A, type = maxlin_model), 
                                k = floor(n_train * k_frac), 
                                method = pca_method, 
                                n_pcs = n_pcs, 
                                plot = FALSE,
                                failure_u = 1))) %>%
  unnest_wider(pca) %>%
  select(-c(Theta_ext_train:lambda)) %>%
  mutate(p_max = sum(apply(A, 2, function(a) max(a / 1))),
         p_min = sum(apply(A, 2, function(a) min(a / 1)))) %>%
  saveRDS(file = file.path("scripts", "compositional", "results", "sim-pca-az.RDS"))


# Husler-Reiss 3D example --------------------------------------------------------------

# beta = 0.5
# CoDA-PCA > DS-PCA -- seed = 5, 12
# CoDA-PCA = DS-PCA -- seed = 3 (very slight curve), 9, 10, 11

# generate variograms
# seed = 5: marked curvature
# seed = 3: very slight curvature
# seed = 11: straight line from vertex to opposite edge
seed_vals <- c(5, 3, 11)
Gamma_list <- lapply(seq_along(seed_vals), function(i) {
  set.seed(seed_vals[i])
  make_hr_variogram(d = 3, beta = 0.5)
})
saveRDS(Gamma_list, file = file.path("scripts", "compositional", "results", "sim-pca-hr-threedim-Gamma.RDS"))

# compute failure probabilities based on u=50 empirically
p_emp <- pbsapply(seq_along(seed_vals), function(i) {
  set.seed(seed_vals[i])
  X_emp <- mev::rmev(n = 10^6, d = 3, sigma = Gamma_list[[i]], model = "hr")
  p_min <- X_emp %>% 
    as.data.frame() %>% 
    filter(rowSums(.) > 100) %>% 
    apply(1, function(x) min(x) > 100) %>% 
    mean()
  p_max <- X_emp %>% 
    as.data.frame() %>% 
    filter(rowSums(.) > 100) %>% 
    apply(1, function(x) max(x) > 100) %>% 
    mean()
  return(c(p_min, p_max))
}) %>%
  set_colnames(seq_along(seed_vals)) %>%
  set_rownames(c("p_min", "p_max"))

# simulate results
set.seed(1)
expand_grid(
  n_train = c(5000),
  n_test = 5000,
  k_frac = c(0.01),
  vario_seed_index = seq_along(seed_vals),
  pca_method = c("ds", "coda"),
  n_pcs = 1:2,
  rep = 1:50
) %>%
  mutate(param_index = row_number()) %>%
  rowwise() %>%
  mutate(pca = list(extreme_pca(X_train = set_colnames(mev::rmev(n = n_train, d = 3, sigma = Gamma_list[[vario_seed_index]], model = "hr"), paste0("X", 1:3)), 
                                X_test = set_colnames(mev::rmev(n = n_test, d = 3, sigma = Gamma_list[[vario_seed_index]], model = "hr"), paste0("X", 1:3)), 
                                k = floor(n_train * k_frac), 
                                method = pca_method, 
                                n_pcs = n_pcs, 
                                plot = FALSE,
                                failure_u = 1))) %>%
  unnest_wider(pca) %>%
  select(-c(Theta_ext_train:lambda)) %>%
  mutate(p_max = p_emp["p_max", vario_seed_index],
         p_min = p_emp["p_min", vario_seed_index]) %>%
  saveRDS(file = file.path("scripts", "compositional", "results", "sim-pca-hr-threedim.RDS"))

# Husler-Reiss high-dimensional example --------------------------------------------------------------

set.seed(3)

d <- 10

Gamma <- bdiag(make_hr_variogram(d = 4, beta = 0.05), 
               make_hr_variogram(d = 4, beta = 2),
               make_hr_variogram(d = 2, beta = 0.03)) %>% 
  as.matrix() %>%
  set_colnames(paste0("X", 1:d)) %>% 
  set_rownames(paste0("X", 1:d))
Gamma[Gamma == 0] <- 10
diag(Gamma) <- 0

X_emp <- mev::rmev(n = 10^6, d = d, sigma = Gamma, model = "hr")
p_min_emp <- X_emp %>% 
  as.data.frame() %>% 
  filter(rowSums(.) > 100) %>% 
  apply(1, function(x) min(x) > 100) %>% 
  mean()
p_max_emp <- X_emp %>% 
  as.data.frame() %>% 
  filter(rowSums(.) > 100) %>% 
  apply(1, function(x) max(x) > 100) %>% 
  mean()

# compute loadings matrix
k <- 200
R_train <- apply(X_emp, 1, pracma::Norm, p = 1)
R_thresh <- Rfast::nth(R_train, k = k + 1, descending = TRUE)
ext_train <- which(R_train > R_thresh)
Theta_ext_train <- X_emp[ext_train, ] / R_train[ext_train]
Theta_ext_train <- acomp(Theta_ext_train)
pc <- princomp(Theta_ext_train)

expand_grid(
  n_train = 10^4,
  n_test = 10^4,
  k_frac = 0.02,
  pca_method = c("ds", "coda"),
  n_pcs = 1:(d-1),
  rep = 1:10) %>%
  filter(!(n_pcs == d & pca_method == "coda")) %>%
  mutate(param_index = row_number()) %>%
  rowwise() %>%
  mutate(pca = list(extreme_pca(X_train = set_colnames(mev::rmev(n = n_train, d = d, sigma = Gamma, model = "hr"), paste0("X", 1:d)), 
                                X_test = set_colnames(mev::rmev(n = n_test, d = d, sigma = Gamma, model = "hr"), paste0("X", 1:d)), 
                                k = floor(n_train * k_frac), 
                                method = pca_method, 
                                n_pcs = n_pcs, 
                                plot = FALSE,
                                failure_u = 1))) %>%
  unnest_wider(pca) %>%
  select(-c(Theta_ext_train:k_test)) %>%
  rowwise() %>% 
  mutate(prop_var = sum(lambda[1:n_pcs]) / sum(lambda)) %>%
  mutate(p_min = p_min_emp,
         p_max = p_max_emp) %>%
  saveRDS(file = file.path("scripts", "compositional", "results", "sim-pca-hr-highdim.RDS"))


saveRDS(pc, file = file.path("scripts", "compositional", "results", "sim-pca-hr-highdim-loadings.RDS"))
saveRDS(Gamma, file = file.path("scripts", "compositional", "results", "sim-pca-hr-highdim-Gamma.RDS"))
# saveRDS(X_emp, file = file.path("scripts", "compositional", "results", "sim-pca-hr-highdim-Xemp.RDS"))



# High-dimensional A example --------------------------------------------------------

set.seed(1)

d <- 10
q <- 50
tau1 <- runif(n = q, -3, 3)
tau2 <- runif(n = q, -3, 3)
u1 <- rbeta(n = d, shape1 = 2, shape2 = 2) %>% divide_by(sum(.))
u2 <- rbeta(n = d, shape1 = 2, shape2 = 5) %>% divide_by(sum(.))

A <- matrix(NA, nrow = d, ncol = q)
for (j in 1:ncol(A)) {
  A[, j] <- perturbe(power.acomp(u1, tau1[j]), power.acomp(u2, tau2[j]))
}
A <- A / rowSums(A)

set.seed(1)

expand_grid(
  n_train = 10^4,
  n_test = 10^4,
  k_frac = 0.02,
  maxlin_model = c("maxlin", "cooley"),
  pca_method = c("ds", "coda"),
  n_pcs = 1:(d-1),
  rep = 1:5
) %>%
  mutate(param_index = row_number()) %>%
  rowwise() %>%
  mutate(pca = list(extreme_pca(X_train = rmaxlinearfactor(n = n_train, A = A, type = maxlin_model), 
                                X_test = rmaxlinearfactor(n = n_test, A = A, type = maxlin_model), 
                                k = floor(n_train * k_frac), 
                                method = pca_method, 
                                n_pcs = n_pcs, 
                                plot = FALSE,
                                failure_u = 1))) %>%
  unnest_wider(pca) %>%
  select(-c(Theta_ext_train:k_test)) %>%
  rowwise() %>% 
  mutate(prop_var = sum(lambda[1:n_pcs]) / sum(lambda)) %>%
  mutate(p_max = sum(apply(A, 2, function(a) max(a / 1))),
         p_min = sum(apply(A, 2, function(a) min(a / 1)))) %>%
  saveRDS(file = file.path("scripts", "compositional", "results", "sim-pca-az-highdim.RDS"))

saveRDS(A, file = file.path("scripts", "compositional", "results", "sim-pca-az-highdim-A.RDS"))

# 
# X_emp <- rmaxlinearfactor(n = 10^6, A = A, type = "maxlin")
# k <- 200
# R_train <- apply(X_emp, 1, pracma::Norm, p = 1)
# R_thresh <- Rfast::nth(R_train, k = k + 1, descending = TRUE)
# ext_train <- which(R_train > R_thresh)
# Theta_ext_train <- X_emp[ext_train, ] / R_train[ext_train]
# Theta_ext_train <- acomp(Theta_ext_train)
# pc <- princomp(Theta_ext_train)
# 
# clrInv(pc$loadings[1:10, 1])
# 
# par(mfrow = c(1, 3))
# for (i in c(1, 4, 7)) {
#   extreme_pca(X_emp[, i:(i+2)], X_emp[, i:(i+2)], k = 200, method = "coda", n_pcs = 1:2, plot = TRUE)
# }
# 
# extreme_pca()