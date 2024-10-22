rm(list = ls())

sapply(list.files(path = "R/changing-ext-dep", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

# Compute CM critical values --------------------------------------------------------

bb_L2 <- readRDS(file = file.path("scripts", "changing-ext-dep", "results", "bb_L2.RDS"))

ks_cv <- sapply(1:100, function(d) CPAT:::qkolmogorov(0.95^(2/(d*(d-1)))))
cm_cv <- sapply(1:100, function(d) quantile(bb_L2, probs = 0.95^(2/(d*(d-1))))) %>% as.numeric()


# Type I error and power ------------------------------------------------------------

## Set up parameter combinations -----------------------------------------------------

params_small_sample <- expand_grid(n = c(2500, 5000, 10000),
                      d = c(2, 5, 10, 25)) %>%
  cross_join(rbind(expand_grid(model = c("log"), 
                               param0 = c(2),
                               param1 = c(2, 3, 4)),
                   expand_grid(model = c("hr"),
                               param0 = c(1),
                               param1 = c(0.5, 1, 2)))) %>%
  mutate(change_type = case_when(param0 == param1 ~ list(c("none")),
                                 .default = list(c("jump", "linear")))) %>%
  unnest(change_type) %>%
  cross_join(tibble(test_method = c("pawley", "drees"))) %>%
  filter(d == 2 | change_type == "none") %>%
  filter(!(d > 2 & test_method == "drees")) %>%
  rowwise() %>%
  mutate(n_blocks = list(c(25, 50))) %>%
  unnest(n_blocks) %>%
  mutate(b = n / n_blocks) %>%
  rowwise() %>%
  mutate(k_frac = list(c(0.05, 0.1, 0.15))) %>%
  unnest(k_frac) %>%
  ungroup() %>%
  mutate(k = floor(k_frac * b)) %>%
  filter(k * n / b > d * (d - 1) / 2) %>% # V rank condition
  filter(k > d) # TPDM rank condition

params_large_sample <- expand_grid(n = c(10^6),
                      d = c(2, 5)) %>%
  cross_join(rbind(expand_grid(model = c("log"), 
                               param0 = c(2),
                               param1 = c(2, 2.5)),
                   expand_grid(model = c("hr"),
                               param0 = c(1),
                               param1 = c(1, 1.5)))) %>%
  mutate(change_type = case_when(param0 == param1 ~ list(c("none")),
                                 .default = list(c("jump")))) %>%
  unnest(change_type) %>%
  cross_join(tibble(test_method = c("pawley", "drees"))) %>%
  filter(!(d > 2 & test_method == "drees")) %>%
  rowwise() %>%
  mutate(n_blocks = list(c(500))) %>%
  unnest(n_blocks) %>%
  mutate(b = n / n_blocks) %>%
  rowwise() %>%
  mutate(k_frac = list(c(0.025))) %>%
  unnest(k_frac) %>%
  ungroup() %>%
  mutate(k = floor(k_frac * b)) %>%
  filter(k * n / b > d * (d - 1) / 2) %>% # V rank condition
  filter(k > d)

params <- rbind(params_small_sample, params_large_sample) %>%
  mutate(param_index = row_number())

saveRDS(params, file = file.path("scripts", "changing-ext-dep", "results", "sim-params.RDS"))

## Run simulations and tests ---------------------------------------------------------

set.seed(8)
pblapply(seq_len(nrow(params)), function(i) {
  
  # get parameter values
  n <- params$n[i]
  d <- params$d[i]
  model <- params$model[i] 
  param0 <- params$param0[i] 
  param1 <- params$param1[i] 
  change_type <- params$change_type[i]
  test_method <- params$test_method[i]
  b <- params$b[i]
  k <- params$k[i]
  test_method <- params$test_method[i]
  
  # set number of repetitions -- 
  # seed1: (40,20,20) took 02:38:12
  # seed2: (100,50,50) took 06:46:22
  # seed3: (110,30,30) took 06:13:41
  # seed4: (250,50,50) took 13:52:00
  # seed5: (150,50,50) took 08:54:15
  # seed6: (150,50,50) took 09:21:38
  # seed7: (150,50,50) took 09:03:25
  # seed8: (50,100,100) took 09:03:25
  nreps <- case_when(
    n < 10^6 & d < 10 ~ 150,
    n < 10^6 & d >= 10 ~ 50, 
    n == 10^6 ~ 50
  )

  # perform simulations and run tests
  out <- replicate(n = nreps, expr = {
    X <- sim_X_changing_dep(n, d, model, param0, param1, change_type)
    if (test_method == "pawley") {
      test_pawley(X, b, k, return_all = FALSE)
    } else {
      test_drees(X, b, k, return_all = FALSE)
    }
  })
  
}) %>%
  lapply(t) %>%
  lapply(function(M) as.data.frame(M) %>% mutate(rep = seq_len(nrow(M)))) %>%
  bind_rows(.id = "param_index") %>%
  mutate(param_index = as.integer(param_index)) %>%
  full_join(params, by = "param_index") %>%
  unnest(c(ks, cm, elapsed_time)) %>%
  pivot_longer(c(ks, cm), names_to = "test_type", values_to = "test_stat") %>%
  mutate(crit_val = case_when(
    test_method == "pawley" & test_type == "ks" ~ ks_cv[d],
    test_method == "pawley" & test_type == "cm" ~ cm_cv[d],
    test_method == "drees" & test_type == "ks" & d == 2 ~ 0.8135, # see Table 1 in Drees (2023)
    test_method == "drees" & test_type == "cm" & d == 2 ~ 0.1939, # see Table 1 in Drees (2023)
    .default = NA_real_
  )) %>%
  mutate(reject_H0 = case_when(test_stat > crit_val ~ TRUE, .default = FALSE)) %>%
  saveRDS(file = file.path("scripts", "changing-ext-dep", "results", "sim-test-results-seed8.RDS"))


# Extra -----------------------------------------------------------------------------

list(file.path("scripts", "changing-ext-dep", "results", "sim-test-results-seed1.RDS"),
     file.path("scripts", "changing-ext-dep", "results", "sim-test-results-seed2.RDS"),
     file.path("scripts", "changing-ext-dep", "results", "sim-test-results-seed3.RDS"),
     file.path("scripts", "changing-ext-dep", "results", "sim-test-results-seed4.RDS")) %>%
  lapply(readRDS) %>%
  bind_rows() %>%
  filter(change_type == "none") %>%
  group_by(n, d, n_blocks, k_frac, test_method, test_type) %>%
  summarise(empirical_type1 = 100 * mean(reject_H0)) %>%
  arrange(n, d, n_blocks, k_frac) %>%
  mutate(k_frac = 100 * k_frac) %>%
  pivot_wider(names_from = c(d, test_method, test_type), values_from = empirical_type1) %>%
  kbl(col.names = c("$n$", "No. blocks, $n/b$", "$k/b$ (%)", rep(c("CM", "KS"), 5)),
      booktabs = TRUE, digits = 1, escape = FALSE) %>%
  add_header_above(c(' ' = 3, 'Drees' = 2, 'Pawley' = 2, 'Pawley' = 2, 'Pawley' = 2, 'Pawley' = 2), line_sep = 1) %>%
  add_header_above(c(' ' = 3, "$d=2$" = 4, "$d=5$" = 2, "$d=10$" = 2, "$d=25$" = 2), line_sep = 1, bold = TRUE, escape = FALSE) %>%
  collapse_rows(columns = 1:2, latex_hline = "major", row_group_label_position = "first") %>%
  kable_styling(latex_options = "striped")


list(file.path("scripts", "changing-ext-dep", "results", "sim-test-results-seed1.RDS"),
     file.path("scripts", "changing-ext-dep", "results", "sim-test-results-seed2.RDS")) %>%
  lapply(readRDS) %>%
  bind_rows() %>%
  group_by(n, d, test_method, test_type, k_frac, n_blocks) %>%
  summarise(average_time = mean(elapsed_time)) %>%
  ggplot(aes(x = d, y = average_time, colour = as.factor(n_blocks), linetype = as.factor(k_frac))) +
  geom_line() + 
  facet_grid(~ n)
