sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/shrinkage-tpdm", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

library(tidyverse)
library(tidyr)
library(magrittr)
library(scales)
library(ggh4x)
library(ggpubr)
library(colorspace)
library(kableExtra)
library(reshape2)
library(mev)
library(GGally)
library(ggforce)

X <- load_red_sea_temp(alpha = 2)$X 
colnames(X) <- paste0("X", "[", 1:ncol(X), "]")

Sigma_hat <- tpdm_estimate(X, k = 0.15 * 1617)
max_sigma <- max(Sigma_hat[upper.tri(Sigma_hat, diag = FALSE)])


lambda_vals <- seq(from = 0, to = 0.5, by = 0.1)

hard_Sigma <- lapply(lambda_vals, function(lambda) shrink_cov(Sigma_hat, shrinkage_type = "hard", lambda = lambda)) %>%
  set_names(lambda_vals)

soft_Sigma <- lapply(lambda_vals, function(lambda) shrink_cov(Sigma_hat, shrinkage_type = "soft", lambda = lambda)) %>%
  set_names(lambda_vals)

adlasso_Sigma <- lapply(lambda_vals, function(lambda) shrink_cov(Sigma_hat, shrinkage_type = "adlasso", lambda = lambda, eta = 2)) %>%
  set_names(lambda_vals)

lw_Sigma <- lapply(lambda_vals, function(lambda) shrink_cov(Sigma_hat, shrinkage_type = "lw", lambda = lambda)) %>%
  set_names(lambda_vals)

saveRDS(list("hard_Sigma" = hard_Sigma, 
             "soft_Sigma" = soft_Sigma,
             "adlasso_Sigma" = adlasso_Sigma,
             "lw_Sigma" = lw_Sigma),
        "scripts/shrinkage-tpdm/results/red-sea-tpdm.RDS")


set.seed(1)
d <- 5
Lambda <- make_hr_variogram(d = d, beta = 3)
X <- rmev(n = 5000, d = d, sigma = Lambda, model = "hr")
Sigma_hat <- tpdm_estimate(X, k = 500)
Sigma <- matrix(hr_tpdm(Lambda), d, d)
lambda_vals <- seq(from = 0, to = 1, by = 0.005)

lapply(lambda_vals, function(lambda) {
  hard <- frob_loss(Sigma, shrink_cov(Sigma_hat, "hard", lambda))  
  soft <- frob_loss(Sigma, shrink_cov(Sigma_hat, "soft", lambda))
  adlasso <- frob_loss(Sigma, shrink_cov(Sigma_hat, "adlasso", lambda, eta = 2))
  lw <- frob_loss(Sigma, shrink_cov(Sigma_hat, "lw", lambda)) 
  return(data.frame("lambda" = lambda,
                    "hard" = hard,
                    "soft" = soft,
                    "adlasso" = adlasso,
                    "lw" = lw))
}) %>% 
  bind_rows() %>% 
  pivot_longer(cols = -lambda, names_to = "method", values_to = "risk") %>%
  group_by(method) %>%
  mutate(lambda_star = lambda[which.min(risk)], min_risk = min(risk)) %>%
  saveRDS("scripts/shrinkage-tpdm/results/frob-risk-example.RDS")

