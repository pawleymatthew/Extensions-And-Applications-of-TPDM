rm(list = ls())

library(cluster)
library(tibble)
library(tidyr)
library(ggplot2)
library(scales)
library(evd)
library(dplyr)
library(purrr)
library(lattice)
library(gridExtra)

# Directory -------------------------------------------------------------------------

data_dir <- file.path("data/eva-data-challenge")
fn_dir <- file.path("R/eva-data-challenge")

# Source functions and load data ----------------------------------------------------

fn_files <- list.files(path = fn_dir, pattern = "*.R$", recursive = TRUE, full.names = TRUE)
sapply(fn_files, source)

sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/shrinkage-tpdm", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

# Load data -------------------------------------------------------------------------

UtopulaU1 <- read.csv(file = file.path(data_dir, "UtopulaU1.csv"), header = TRUE)
UtopulaU2 <- read.csv(file = file.path(data_dir, "UtopulaU2.csv"), header = TRUE)

colnames(UtopulaU1) <- paste(colnames(UtopulaU1), "U1", sep = ".")
colnames(UtopulaU2) <- paste(colnames(UtopulaU2), "U2", sep = ".")

Utopula <- cbind(UtopulaU1, UtopulaU2)
Utopula <- as.matrix(Utopula)

n <- nrow(Utopula)
d <- ncol(Utopula)

Utopula <- exp(Utopula / 2) # Frechet(2) margins


# Thresholds ------------------------------------------------------------------------

u1 <- function(varnames) dplyr::case_when(endsWith(varnames, "U1") ~ evd::qfrechet(299/300, shape = 2),
                                          endsWith(varnames, "U2") ~ evd::qfrechet(288/300, shape = 2))
u2 <- evd::qfrechet(299/300, shape = 2)

# C4 --------------------------------------------------------------------------------

## Clustering ------------------------------------------------------------------------

K <- 5

F_Utopula <- evd::pfrechet(Utopula, shape = 2)
D <- dist(t(F_Utopula), method = "manhattan") / (2 * n)
Utopula_cluster <- cluster::pam(x = D, k = K, diss = TRUE, cluster.only = TRUE)
cluster_dims <- table(Utopula_cluster)

X_clusters <- lapply(1:K, function(l) Utopula[, Utopula_cluster == l])


# Strategy for various k ------------------------------------------------------------

set.seed(1)

k_frac_vals <- seq(from = 0.02, to = 0.07, by = 0.005)
k_vals <- floor(k_frac_vals) * n 
nrep <- 10

res <- pblapply(k_frac_vals, function(k_frac) {
  
  k_fun_tmp <- function(n) floor(100 * k_frac * sqrt(n))
  k_tmp <- k_fun_tmp(n)
  
  # STRATEGY 1: empirical estimation approach
  A_hat_clusters <- lapply(X_clusters, function(X_l) max_linear_A_emp(X_l, k_tmp, alpha = 2))
  p_u1_emp_clusters <- lapply(A_hat_clusters, function(A) max_linear_p_hat(A = A, u = u1(rownames(A)), alpha = 2)) 
  p_u1_emp <- p_u1_emp_clusters %>% purrr::reduce(prod)
  p_u2_emp_clusters <- lapply(A_hat_clusters, function(A) max_linear_p_hat(A = A, u = u2, alpha = 2)) 
  p_u2_emp <- p_u2_emp_clusters %>% purrr::reduce(prod)
  
  # STRATEGY 2: sparse empirical estimation approach
  Utopula1 <- Utopula^2 # Frechet shape 1
  
  X1_clusters <- lapply(1:K, function(l) Utopula1[, Utopula_cluster == l])
  A_hatstar_clusters <- lapply(X1_clusters, function(X_l) {
    R_l <- apply(X_l, 1, pracma::Norm, p = 1)
    R_l_kplus1 <- Rfast::nth(R_l, k = k_tmp + 1, descending = TRUE)
    ext_ind <- which(R_l > R_l_kplus1)
    piTheta_ext_l <- (X_l[ext_ind, ] / R_l_kplus1) %>% apply(1, euc_proj) %>% t()
    d_l <- ncol(X_l)
    A_hatstar_l <- (d_l / k_tmp) *  t(piTheta_ext_l)
    return(A_hatstar_l)
  })
  
  p_u1_spemp_clusters <- lapply(A_hatstar_clusters, function(A) max_linear_p_hat(A = A, u = u1(rownames(A))^2, alpha = 1)) 
  p_u1_spemp <- p_u1_spemp_clusters %>% purrr::reduce(prod)
  p_u2_spemp_clusters <- lapply(A_hatstar_clusters, function(A) max_linear_p_hat(A = A, u = u2^2, alpha = 1)) 
  p_u2_spemp <- p_u2_spemp_clusters %>% purrr::reduce(prod)
  
  # STRATEGY 3: CP estimation approach with empirical TPDM

  Sigma_hat_clusters <- lapply(X_clusters, function(X_l) tpdm_estimate(X_l, k_tmp))
  
  Theta_hat_kz <- tibble("cluster" = rep(1:K, cluster_dims), "path_start" = sequence(cluster_dims)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate("cp_decomp_kz" = list(cp_decomp_kz_search(Sigma_hat_clusters[[cluster]], path_start = path_start, max_paths = 100))) %>%
    tidyr::unnest_wider(cp_decomp_kz) %>%
    dplyr::rowwise() %>%
    dplyr::mutate("p_u1" = max_linear_p_hat(A = A, beta = 1:nrow(A), u = u1(names(which(Utopula_cluster == cluster))), alpha = 2)) %>%
    dplyr::mutate("p_u2" = max_linear_p_hat(A = A, u = u2, alpha = 2)) %>%
    dplyr::ungroup()
  
  p_u1_cp <- Theta_hat_kz %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::reframe("p_u1" = list(unique(p_u1))) %>% 
    dplyr::pull(p_u1) %>% 
    expand.grid() %>% 
    apply(1, prod)
  
  p_u2_cp <- Theta_hat_kz %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::reframe("p_u2" = list(unique(p_u2))) %>% 
    dplyr::pull(p_u2) %>% 
    expand.grid() %>% 
    apply(1, prod)
  
  # STRATEGY 4: CP estimation approach with Ledoit-Wolf TPDM
  
  lambda_hat_clusters <- replicate(n = nrep, expr = {
    lapply(X_clusters, function(X_l) lw_ratio_method(X_l, k_fun = k_fun_tmp, kappa_vals = 0.75, B = 5)) %>%
      bind_rows() %>%
      pull(lambda_hat)
  })
    lambda_hat_tmp <- lambda_hat_clusters %>%
    t() %>%
    apply(2, function(z) quantile(z, 0.9))

  # lambda_hat_clusters <- lapply(X_clusters, function(X_l) lw_ratio_method(X_l, k_fun = k_fun_tmp, kappa_vals = 0.75, B = 5)) %>%
  #     bind_rows() %>%
  #     pull(lambda_hat)
  # 
  # lambda_hat_tmp <- lambda_hat_clusters
  
  lambda_hat_tmp[lambda_hat_tmp > 0.999] <- 0.999
  
  Sigma_tilde_clusters <- lapply(1:K, function(l) {
    shrink_cov(Sigma_hat_clusters[[l]], "lw", lambda = lambda_hat_tmp[l])
  })
  
  Theta_tilde_kz <- tibble("cluster" = rep(1:K, cluster_dims), "path_start" = sequence(cluster_dims)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate("cp_decomp_kz" = list(cp_decomp_kz_search(Sigma_tilde_clusters[[cluster]], path_start = path_start, max_paths = 100))) %>%
    tidyr::unnest_wider(cp_decomp_kz) %>%
    dplyr::rowwise() %>%
    dplyr::mutate("p_u1" = max_linear_p_hat(A = A, beta = 1:nrow(A), u = u1(names(which(Utopula_cluster == cluster))), alpha = 2)) %>%
    dplyr::mutate("p_u2" = max_linear_p_hat(A = A, u = u2, alpha = 2)) %>%
    dplyr::ungroup()
  
  p_u1_tilde <- Theta_tilde_kz %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::reframe("p_u1" = list(unique(p_u1))) %>% 
    dplyr::pull(p_u1) %>% 
    expand.grid() %>% 
    apply(1, prod)
  
  p_u2_tilde <- Theta_tilde_kz %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::reframe("p_u2" = list(unique(p_u2))) %>% 
    dplyr::pull(p_u2) %>% 
    expand.grid() %>% 
    apply(1, prod)
  
  # information about the shrinkage
  Sigma_summary <- Sigma_hat_clusters %>% 
    lapply(FUN = function(S) unlist(S[upper.tri(S, diag = F)])) %>% 
    enframe(name = "cluster", value = "sigma") %>%
    unnest(sigma) %>%
    group_by(cluster) %>%
    summarise(min_sigma = min(sigma), median_sigma = median(sigma), max_sigma = max(sigma)) 
  
  return(tibble("k_frac" = k_frac,
                "empiricalstandard_p1" = p_u1_emp,
                "empiricalstandard_p2" = p_u2_emp,
                "empiricalsparse_p1" = p_u1_spemp,
                "empiricalsparse_p2" = p_u2_spemp,
                "cpempirical_p1" = median(p_u1_cp),
                "cpempirical_p2" = median(p_u2_cp),
                "cplw_p1" = median(p_u1_tilde),
                "cplw_p2" = median(p_u2_tilde),
                "lambda_hat" = list(lambda_hat_tmp),
                "min_sigma" = list(pull(Sigma_summary, min_sigma)),
                "median_sigma" = list(pull(Sigma_summary, median_sigma)),
                "max_sigma" = list(pull(Sigma_summary, max_sigma))))
  
}) %>% bind_rows()

res %>% saveRDS("scripts/shrinkage-tpdm/results/eva-c4-ledoit-wolf.RDS")


fancy_scientific_short <- function(l) {
  l %>%
    format(scientific = TRUE) %>%
    gsub("^(.*)e", "e", .) %>%
    gsub("e", "10^", .) %>%
    parse(text = .)
}

res %>% 
  pivot_longer(
    cols = empiricalstandard_p1:cplw_p2,
    names_to = c("method", "threshold"), 
    names_sep = "_",
    values_to = "estimate") %>%
  mutate(true_prob = case_when(threshold == "p1" ~ 8.4e-23,
                               threshold == "p2" ~ 5.4e-25)) %>%
  mutate(threshold = case_when(threshold == "p1" ~ "p[1]",
                               threshold == "p2" ~ "p[2]")) %>%
  ggplot(aes(x = k_frac, y = estimate, colour = method)) +
  geom_smooth(se = FALSE) +
  geom_hline(aes(yintercept = true_prob), linetype = "dashed") +
  geom_vline(xintercept = 0.025, linetype = "dashed") +
  facet_grid(. ~ threshold, labeller = label_parsed) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)), labels = label_percent()) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)), transform = log10_trans(), breaks = breaks_log(n = 7), labels = fancy_scientific_short) +
  scale_colour_manual(values = c("grey", "red", "blue", "darkgreen"), labels = c("CP empirical TPDM", "CP Ledoit-Wolf TPDM", "Sparse empirical", "Empirical")) +
  scale_linetype_manual(values = 1:2, labels = c("p1", "p2")) +
  labs(x = expression(k/n),
       y = expression(hat(p)), 
       colour = "Method",
       linetype = "Threshold") +
  theme_light() +
  theme(legend.position = "top")

res %>% 
  unnest_longer(lambda_hat, indices_to = "cluster", values_to = "lambda_hat_cluster") %>%
  select(k_frac, cluster, lambda_hat_cluster) %>%
  mutate(cluster = as.factor(cluster)) %>%
  ggplot(aes(x = k_frac, y = lambda_hat_cluster, colour = cluster)) +
  geom_smooth(se = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)), labels = label_percent()) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.01))) +
  scale_colour_manual(values = c("grey", "red", "blue", "darkgreen", "orange"), labels = 1:5) +
  scale_linetype_manual(values = 1:2, labels = c("p1", "p2")) +
  labs(x = expression(k/n),
       y = expression(hat(lambda)), 
       colour = "Cluster") +
  theme_light()


# Tables and figures for Rmd document -----------------------------------------------

# table of summary statistics
library(kableExtra)
Sigma_hat_clusters %>% 
  lapply(FUN = function(S) unlist(S[upper.tri(S, diag = F)])) %>% 
  enframe(name = "cluster", value = "sigma") %>%
  unnest(sigma) %>%
  group_by(cluster) %>%
  summarise(min_sigma = min(sigma), median_sigma = median(sigma), max_sigma = max(sigma)) %>%
  mutate(dl = as.numeric(table(Utopula_cluster))) %>%
  mutate(U1_sites = as.numeric(table(Utopula_cluster[grepl("U1$", names(Utopula_cluster))]))) %>%
  mutate(U2_sites = as.numeric(table(Utopula_cluster[grepl("U2$", names(Utopula_cluster))]))) %>%
  mutate(lambda_hat_l = lambda_hat_tmp) %>%
  relocate(dl, U1_sites, U2_sites, lambda_hat_l, .after = cluster) %>%
  kbl(col.names = c("Cluster", "Size", "U1 sites", "U2 sites", "Shrinkage", "Min.", "Median", "Max."),
      booktabs = TRUE, digits = c(rep(0, 5), rep(2, 3)), escape = FALSE) %>%
  add_header_above(c(' ' = 5, "$\\\\{\\\\sigma_{ij}:i\\\\neq j\\\\}$" = 3), line_sep = 1, escape = FALSE) %>%
  kable_styling(latex_options = "striped")

colpal <- colorspace::sequential_hcl(n = 100, "Viridis")
breaks <- lattice::do.breaks(c(0, 1.3), length(colpal))

library(Matrix)
Sigma_hat_bdiag <- bdiag(Sigma_hat_clusters) %>% as.matrix()
S1 <- t(Sigma_hat_bdiag[nrow(Sigma_hat_bdiag):1, ])
p1 <- lattice::levelplot(S1, col.regions = colpal, at = breaks,
                         xlab = "", ylab = "",
                         colorkey = TRUE,
                         scales = list(x = list(draw = FALSE),
                                       y = list(draw = FALSE)))

Sigma_tilde_bdiag <- bdiag(Sigma_tilde_clusters) %>% as.matrix()
S2 <- t(Sigma_tilde_bdiag[nrow(Sigma_tilde_bdiag):1, ])
p2 <- lattice::levelplot(S2, col.regions = colpal, at = breaks,
                         xlab = "", ylab = "",
                         colorkey = TRUE,
                         scales = list(x = list(draw = FALSE),
                                       y = list(draw = FALSE)))
grid.arrange(p1, p2, ncol = 2)

# table of estimates and true values
data.frame("method" = c("Truth", 
                        "CP-decomposition of empirical TPDM",
                        "CP-decomposition of Ledoit-Wolf TPDM",
                        "Empirical estimation",
                        "Sparse empirical estimation"),
           "p1" = c(8.4e-23, median(p_u1_cp),  median(p_u1_tilde), p_u1_emp, p_u1_spemp),
           "p2" = c(5.4e-25, median(p_u2_cp),  median(p_u2_tilde), p_u2_emp, p_u2_spemp)) %>%
  kbl(col.names = c("Method", "$p_1$", "$p_2$"),
      booktabs = TRUE, format.args = list(scientific = TRUE, digits = 4), escape = FALSE) %>%
  add_header_above(c(' ' = 5, "$\\\\{\\\\sigma_{ij}:i\\\\neq j\\\\}$" = 3), line_sep = 1, escape = FALSE) %>%
  kable_styling(latex_options = "striped")

