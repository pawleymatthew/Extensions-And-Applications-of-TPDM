library(tidyverse)
library(tidyr)
library(magrittr)
#library(SimCop)
library(pracma)
library(expm)
library(sn)
library(scales)
library(ggh4x)
#library(CPAT)
library(maps)
library(patchwork)
library(ggpubr)
library(colorspace)
library(pbapply)
library(kableExtra)
library(reshape2)
library(e1071) # simulate Brownian bridge
library(tictoc)

# CPAT not on CRAN anymore, run this to install CPAT from archive if needed
#install.packages("CPAT_0.1.0.tar", repos = NULL, type = "source")


rm(list = ls())

sapply(list.files(path = "R/changing-ext-dep", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/background", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

# new version of this function
# for convenience, just overwrite the one created in R/changing_ext_dep/sim_data_changing_dep.R
sim_X_changing_dep <- function(n, d, model, param0, param1, change_type = "jump", Xonly = TRUE) {
  
  # make base variogram if HR model
  if (model == "hr") {
    Gamma0 <- make_hr_variogram(d = d)
  }
  
  # helper function to simulate n samples from model with given parameter
  sim_X <- function(n, param) {
    if (model == "hr") {
      X <- mev::rmev(n = n, d = d, sigma = param * Gamma0, model = "hr")
    } else {
      X <- mev::rmev(n = n, d = d, param = param, model = model)
    }
    return(as.data.frame(X))
  }
  # simulate data according to change_type
  if (param0 == param1) { # no change in dependence
    X <- sim_X(n = n, param = param0)
  } else if (change_type == "jump") { # jump change in dependence
    n0 <- floor(0.5 * n)
    n1 <- n - n0
    X <- rbind(sim_X(n = n0, param = param0), sim_X(n = n1, param = param1))
  } else if (change_type == "linear") { # linear change in dependence
    X <- lapply((1:n)/n, function(t) {
      X_t <- sim_X(n = 1, param = param0 + t * (param1 - param0))
    }) %>%
      bind_rows()
  } else {
    stop("Specify valid dependence change type.")
  }
  X <- sqrt(X)
  if (Xonly | model != "hr") {
    return(X)
  } else {
    return(list("X" = X, "Gamma0" = Gamma0))
  }
}

# Load critical values and data ----------------------------------------------------------------------

bb_L2 <- readRDS(file = file.path("scripts", "changing-ext-dep", "results", "bb_L2.RDS"))

ks_cv <- sapply(1:100, function(d) CPAT:::qkolmogorov(0.95^(2/(d*(d-1)))))
cm_cv <- sapply(1:100, function(d) quantile(bb_L2, probs = 0.95^(2/(d*(d-1))))) %>% as.numeric()


# Simulate data ---------------------------------------------------------------------

set.seed(1)

#redsea <- load_red_sea_temp(alpha = 2)
param0 = 1
param1 = 4
d = 40
dat <- sim_X_changing_dep(n = 5e3, d = d, model = "hr", param0 = param0, param1 = param1, change_type = "linear", Xonly = FALSE)

X = dat$X
Gamma_t0 = param0 * dat$Gamma0
Gamma_t1 = param1 * dat$Gamma0

# Resample sites and run tests ------------------------------------------------------

set.seed(1)

d_vals <- c(2, 5, 10)

data <- lapply(seq_along(d_vals), function(d_ind) {
  
  d <- d_vals[d_ind]
  nreps <- 1000
  
  pbreplicate(n = nreps, expr = {
    
    sample_sites <- sort(sample(names(X), size = d, replace = FALSE))
    
    X[, sample_sites] %>%
      as_tibble() %>%
      test_pawley(b = 200, k = 20, return_all = FALSE)
    
  }) %>%
    t() %>%
    as_tibble() %>%
    mutate(d = d, region = "north") %>%
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
    )) %>%
    rowwise() %>%
    mutate(p_value = case_when(
      test_type == "ks" ~ 1 - CPAT:::pkolmogorov(test_stat)^(d*(d-1)/2), 
      test_type == "cm" ~ 1 - mean(bb_L2 < test_stat)^(d*(d-1)/2), 
    ))
}) %>%
  bind_rows()

label_data <- data %>%
  group_by(d, region, test_type) %>%
  summarise(rejection_rate = mean(p_value < 0.05)) %>%
  ungroup() %>%
  mutate(facet_var = interaction(region, test_type, d))

pdf("figures/paper-fig-sim-study-p-values.pdf", width = 6, height = 3.5)
data %>%
  ggplot(aes(x = p_value)) +
  geom_histogram(breaks = seq(0, 1, 0.025), color = "darkgrey", fill = "grey", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, colour = "darkblue", linetype = "dashed", linewidth = 0.3) +
  facet_nested(test_type ~ fct_inorder(paste0("d = ", d)), scales = "free", labeller = labeller(region = str_to_title, test_type = toupper), nest_line = element_line(colour = "white")) +
  geom_text(data = label_data, aes(label = paste0(round(100 * rejection_rate, 3), "%"), group = facet_var), x = 0.5, y = Inf, vjust = 3, size = 3.5) +
  scale_x_continuous(expand = c(0, 0), breaks = breaks_extended(n = 6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  xlab("p-value") +
  ylab("Frequency density") +
  theme_light() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
dev.off()

# plot the gamma matrices, not the TPDMs
# log transformed the off-diagonals to compress the scale a bit

pdf("figures/paper-fig-sim-study-hr-params.pdf", width = 6, height = 4)
plot_tpdm(list("t = 0" = log1p(Gamma_t0), "t = 1" = log1p(Gamma_t1)), x_labels = FALSE, y_labels = FALSE)
dev.off()


