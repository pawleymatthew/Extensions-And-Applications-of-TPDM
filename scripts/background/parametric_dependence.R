library(tidyverse)
library(magrittr)
library(mev)
library(pbapply)

sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/background", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

# Empirical estimates of TPDM under logistic and HR models ------------------------------------------------------

n <- 5 * 10^5
k_tpdm <- 100
k_chi <- 250
nrep <- 2
a_vals <- seq(from = 0.025, to = 0.975, by = 0.025)
lambda_vals <- seq(from = 0.1, to = 3, by = 0.1)

sl_res <- pblapply(a_vals, function(a) {
  sapply(seq_len(nrep), function(i) {
    X <- rmev(n = n, d = 2, par = 1 / a, model = "log") %>% sqrt()
    sigma <- tpdm_estimate(X = X, k = k_tpdm) %>% as.numeric() %>% extract(2)
    chi <- taildep(X, u = 1 - k_chi/n, meas = "chi")$chi[, "coef"] %>% as.numeric()
    return(c("model" = "log","dep_par" = a, "rep" = i, "chi" = chi, "sigma" = sigma))
  })
}) %>%
  lapply(t) %>%
  lapply(as.data.frame) %>%
  bind_rows()

hr_res <- pblapply(lambda_vals, function(lambda) {
  sapply(seq_len(nrep), function(i) {
    X <- rmev(n = n, d = 2, sigma = rbind(c(0, lambda^2), c(lambda^2, 0)), model = "hr") %>% sqrt()
    sigma <- tpdm_estimate(X = X, k = k_tpdm) %>% as.numeric() %>% extract(2)
    chi <- taildep(X, u = 1 - k_chi/n, meas = "chi")$chi[, "coef"] %>% as.numeric()
    return(c("model" = "hr","dep_par" = lambda, "rep" = i, "chi" = chi, "sigma" = sigma))
  })
}) %>%
  lapply(t) %>%
  lapply(as.data.frame) %>%
  bind_rows()

rbind(sl_res, hr_res) %>%
  mutate(dep_par = as.numeric(dep_par), chi = as.numeric(chi), sigma = as.numeric(sigma)) %>%
  saveRDS(file = file.path("scripts", "background", "results", "parametric_chi_tpdm_empirical.RDS"))

# Empirical estimates of TPDM under logistic and HR models with small n ------------------------------------------------------

n <- 5000
k_tpdm <- 250
k_chi <- 250
nrep <- 10
a_vals <- seq(from = 0.025, to = 0.975, by = 0.025)
lambda_vals <- seq(from = 0.1, to = 3, by = 0.1)


sl_res <- pblapply(a_vals, function(a) {
  sapply(seq_len(nrep), function(i) {
    X <- rmev(n = n, d = 2, par = 1 / a, model = "log") %>% sqrt()
    sigma <- tpdm_estimate(X = X, k = k_tpdm) %>% as.numeric() %>% extract(2)
    chi <- taildep(X, u = 1 - k_chi/n, meas = "chi")$chi[, "coef"] %>% as.numeric()
    return(c("model" = "log","dep_par" = a, "rep" = i, "chi" = chi, "sigma" = sigma))
  })
}) %>%
  lapply(t) %>%
  lapply(as.data.frame) %>%
  bind_rows()

hr_res <- pblapply(lambda_vals, function(lambda) {
  sapply(seq_len(nrep), function(i) {
    X <- rmev(n = n, d = 2, sigma = rbind(c(0, lambda^2), c(lambda^2, 0)), model = "hr") %>% sqrt()
    sigma <- tpdm_estimate(X = X, k = k_tpdm) %>% as.numeric() %>% extract(2)
    chi <- taildep(X, u = 1 - k_chi/n, meas = "chi")$chi[, "coef"] %>% as.numeric()
    return(c("model" = "hr","dep_par" = lambda, "rep" = i, "chi" = chi, "sigma" = sigma))
  })
}) %>%
  lapply(t) %>%
  lapply(as.data.frame) %>%
  bind_rows()

rbind(sl_res, hr_res) %>%
  mutate(dep_par = as.numeric(dep_par), chi = as.numeric(chi), sigma = as.numeric(sigma)) %>%
  saveRDS(file = file.path("scripts", "background", "results", "parametric_chi_tpdm_empirical_smalln.RDS"))


