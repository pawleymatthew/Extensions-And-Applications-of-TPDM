rm(list = ls())

sapply(list.files(path = "R/changing-ext-dep", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

# Compute CM critical values --------------------------------------------------------

bb_L2 <- readRDS(file = file.path("scripts", "changing-ext-dep", "results", "bb_L2.RDS"))

ks_cv <- sapply(1:100, function(d) CPAT:::qkolmogorov(0.95^(2/(d*(d-1)))))
cm_cv <- sapply(1:100, function(d) quantile(bb_L2, probs = 0.95^(2/(d*(d-1))))) %>% as.numeric()


n <- 3000
k_total_fracs <- seq(from = 0.05, to = 0.2, by = 0.05)
b_vals <- c(50, 100, 150, 200, 250, 300, 500, 600)

res_alt <- pbsapply(k_total_fracs, function(k_total_frac) {
  sapply(b_vals, function(b) {
    replicate(n = 100, expr = {
      k <- floor(k_total_frac * b)
      X <- sim_X_changing_dep(n = n, d = 2, model = "hr", param0 = 1, param1 = 1.5, change_type = "jump")
      test_pawley(X, b, k, return_all = FALSE)$cm > 0.460
    }) %>% mean()
  }) 
}) %>%
  melt(res, value.name = "power", varnames = c("b", "k_frac")) 

res_null <- pbsapply(k_total_fracs, function(k_total_frac) {
  sapply(b_vals, function(b) {
    replicate(n = 100, expr = {
      k <- floor(k_total_frac * b)
      X <- sim_X_changing_dep(n = n, d = 2, model = "hr", param0 = 1, param1 = 1, change_type = "none")
      test_pawley(X, b, k, return_all = FALSE)$cm > 0.460
    }) %>% mean()
  }) 
}) %>%
  melt(res, value.name = "power", varnames = c("b", "k_frac")) 

list(res_alt, res_null) %>%
  bind_rows(.id = "scenario") %>%
  mutate(scenario = case_when(scenario == 1 ~ "alt", .default = "null")) %>%
  mutate(b = b_vals[b], k_frac = k_total_fracs[k_frac]) %>%
  ggplot(aes(x = b, y = power, colour = scenario, linetype = scenario)) +
  geom_path() +
  geom_point() + 
  facet_grid(~ as.factor(k_frac)) + 
  xlab("Block size") +
  ylab("Rejection rate") +
  labs(colour = "Scenario", linetype = "Scenario") +
  scale_color_manual(labels = c("Alternative", "Null"), values = c("red", "blue")) +
  scale_linetype_manual(labels = c("Alternative", "Null"), values = c(1, 2)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme_light()
