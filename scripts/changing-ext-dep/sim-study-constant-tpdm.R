rm(list = ls())

sapply(list.files(path = "R/changing-ext-dep", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)


set.seed(1)

# Find A matrices ----------------------------------------------------------

nrep <- 50000
q <- 20

ml_params <- tibble(A = lapply(1:nrep,
                               function(i) matrix(runif(2 * q), nrow = 2, ncol = q) %>%
                                 multiply_by(1 / sqrt(rowSums(.^2))))) %>%
  rowwise() %>%
  mutate(sigma = (A %*% t(A))[1, 2]) %>%
  mutate(v = 2 * sum(apply(A, 2, function(x) prod(x^2) / sum(x^2))) - sigma^2) %>%
  ungroup()

A <- ml_params %>%
  filter(abs(sigma - 0.8) < 0.0005 & abs(v - 0.06) < 0.00015) %>%
  pull(A) 

# Generate dataset ------------------------------------------------------------------

set.seed(1)
n <- 10000
d <- 2
X <- rbind(SpatialExtremes::rmaxlin(n = n/2, dsgn.mat = A[[1]]),
           SpatialExtremes::rmaxlin(n = n/2, dsgn.mat = A[[2]]))
colnames(X) <- paste0("X", 1:2)
X <- as.data.frame(X)

# Do tests --------------------------------------------------------------------------

test_drees(X, b = 400, k = 40, return_all = TRUE) %>%
  saveRDS(file = file.path("scripts", "changing-ext-dep", "results", "constant-tpdm-drees-example.RDS"))

test_pawley(X, b = 400, k = 40, return_all = TRUE) %>%
  saveRDS(file = file.path("scripts", "changing-ext-dep", "results", "constant-tpdm-pawley-example.RDS"))

# Empirical power of tests --------------------------------------------------------------------------

bb_L2 <- readRDS(file = file.path("scripts", "changing-ext-dep", "results", "bb_L2.RDS"))

ks_cv <- sapply(1:100, function(d) CPAT:::qkolmogorov(0.95^(2/(d*(d-1)))))
cm_cv <- sapply(1:100, function(d) quantile(bb_L2, probs = 0.95^(2/(d*(d-1))))) %>% as.numeric()


drees_res <- pbreplicate(n = 1000, expr = {
  X <- rbind(SpatialExtremes::rmaxlin(n = n/2, dsgn.mat = A[[1]]),
             SpatialExtremes::rmaxlin(n = n/2, dsgn.mat = A[[2]]))
  colnames(X) <- paste0("X", 1:2)
  X <- as.data.frame(X)
  test_drees(X, b = 400, k = 40)
}) %>%
  t() %>%
  as_tibble() %>%
  mutate(test_method = "drees") %>%
  relocate(test_method) %>%
  unnest(c(ks, cm, elapsed_time)) %>%
  pivot_longer(c(ks, cm), names_to = "test_type", values_to = "test_stat") %>%
  mutate(crit_val = case_when(
    test_method == "pawley" & test_type == "ks" ~ ks_cv[d],
    test_method == "pawley" & test_type == "cm" ~ cm_cv[d],
    test_method == "drees" & test_type == "ks" & d == 2 ~ 0.8135, # see Table 1 in Drees (2023)
    test_method == "drees" & test_type == "cm" & d == 2 ~ 0.1939, # see Table 1 in Drees (2023)
    .default = NA_real_
  ))

pawley_res <- pbreplicate(n = 1000, expr = {
  X <- rbind(SpatialExtremes::rmaxlin(n = n/2, dsgn.mat = A[[1]]),
             SpatialExtremes::rmaxlin(n = n/2, dsgn.mat = A[[2]]))
  colnames(X) <- paste0("X", 1:2)
  X <- as.data.frame(X)
  test_pawley(X, b = 400, k = 40)
}) %>%
  t() %>%
  as_tibble() %>%
  mutate(test_method = "pawley") %>%
  relocate(test_method) %>%
  unnest(c(ks, cm, elapsed_time)) %>%
  pivot_longer(c(ks, cm), names_to = "test_type", values_to = "test_stat") %>%
  mutate(crit_val = case_when(
    test_method == "pawley" & test_type == "ks" ~ ks_cv[d],
    test_method == "pawley" & test_type == "cm" ~ cm_cv[d],
    test_method == "drees" & test_type == "ks" & d == 2 ~ 0.8135, # see Table 1 in Drees (2023)
    test_method == "drees" & test_type == "cm" & d == 2 ~ 0.1939, # see Table 1 in Drees (2023)
    .default = NA_real_
  ))

rbind(drees_res, pawley_res) %>%
  group_by(test_method, test_type) %>%
  summarise(power = 100 * mean(test_stat > crit_val)) %>% 
  saveRDS(file = file.path("scripts", "changing-ext-dep", "results", "constant-tpdm-power.RDS"))

# Drees CM: 100%
# Drees KS: 99.8%
# Pawley CM: 5.5%
# Pawley KS: 3.5%

# Coverage of asymptotic CI ---------------------------------------------------------

set.seed(1)

lower_ci <- 0.8 - qnorm(0.975) * sqrt(0.060/40)
upper_ci <- 0.8 + qnorm(0.975) * sqrt(0.060/40)

sigmas <- pbreplicate(n = 1000, expr = {
  X <- rbind(SpatialExtremes::rmaxlin(n = n/2, dsgn.mat = A[[1]]),
             SpatialExtremes::rmaxlin(n = n/2, dsgn.mat = A[[2]]))
  colnames(X) <- paste0("X", 1:2)
  X <- as.data.frame(X)
  res <- test_pawley(X, b = 400, k = 40, return_all = TRUE)
  pull(res$data, sigma)
}) %>%
  as.numeric()

emp_cov <- mean(sigmas > lower_ci & sigmas < upper_ci) # 93.56%

