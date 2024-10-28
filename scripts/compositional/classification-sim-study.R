rm(list = ls())

sapply(list.files(path = "R/compositional", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

library(tidyverse)
library(magrittr)
library(pbapply)
library(caret)
library(Compositional)
library(CompositionalML)
library(ggh4x)
library(scales)
library(ggpubr)



# Define dependence models ----------------------------------------------------------

# define dependence models to generate data from
# see sim_X_classify for interpretation of dependence parameters

model_params <- expand_grid(
  model = c("log", "bilog", "neglog"),
  param1_frac = c(1.5, 3)) %>%
  mutate(param0 = case_when(model == "log" ~ 0.9,
                            model == "bilog" ~ 0.3,
                            model == "neglog" ~ 0.5)) %>%
  mutate(param1 = param0 * param1_frac,
         d = case_when(model == "bilog" ~ 5, .default = 5)) %>%
  mutate(nice_model = case_when(model == "log" ~ "Symmetric",
                           model == "bilog" ~ "Bi",
                           model == "neglog" ~ "Negative")) %>%
  mutate(model_index = as.character(row_number())) %>%
  relocate(model_index, model, nice_model, d, param0, param1, param1_frac)

# Simulate validation sets ----------------------------------------------------------

# test sets from the exact angular measure
# cols: X1 ... Xd class
# data are angles on the L1-simplex

test_data <- pblapply(seq_len(nrow(model_params)), function(i) {
  sim_X_classify(n = 10^5,
                 d = model_params$d[i],
                 param0 = model_params$param0[i],
                 param1 = model_params$param1[i],
                 model = model_params$model[i],
                 type = "rmevspec")
})


# Ternary plots ---------------------------------------------------------------------

library(ggtern)

test_data %>%
  bind_rows(.id = "model_index") %>%
  full_join(model_params, by = "model_index") %>%
  mutate(param1_frac = paste0("vartheta[1]/vartheta[0] == ", param1_frac)) %>%
  ggtern(aes(X1, X2, X3, colour = class)) +
  geom_point(alpha = 0.4, shape = 4) +
  facet_grid(param1_frac ~ nice_model, labeller = label_parsed) +
  scale_colour_manual(values = c("red", "blue"), labels = c(expression(vartheta[0]), expression(vartheta[1]))) +
  xlab("") +
  ylab("") +
  labs(colour = "Class dependence parameter") +
  Llab(expression(X[1])) +
  Tlab(expression(X[2])) +
  Rlab(expression(X[3])) +
  theme_hidegrid() +
  theme_hidelabels() +
  theme(legend.position = "bottom")

detach(package:ggtern, unload = TRUE)


# Run simulations -------------------------------------------------------------------

set.seed(2)

n_train_vals <- 5000
k_frac_vals <- 0.05
nrep <- 30

knn_a_vals <- seq(from = 0, to = 1, by = 0.1)
svm_a_vals <- seq(from = 0, to = 1, by = 0.2)
rf_a_vals <- seq(from = 0, to = 1, by = 0.1)
knn_k_vals <- 1:25

params <- model_params %>%
  cross_join(expand_grid(n_train = n_train_vals, 
                         k_frac = k_frac_vals, 
                         rep = seq_len(nrep)))

res <- pblapply(seq_len(nrow(params)), function(i) {
  
  # generate training data
  train_data <- sim_X_classify(n = params$n_train[i], 
                               d = params$d[i],
                               param0 = params$param0[i],
                               param1 = params$param1[i],
                               model = params$model[i],
                               type = "rmev",
                               k_frac = params$k_frac[i])
  
  x <- train_data %>% select(starts_with("X"))
  y <- train_data %>% pull(class)
  xnew <- test_data %>%
    extract2(as.numeric(params$model_index[i])) %>%
    select(starts_with("X"))
  ynew <- test_data %>%
    extract2(as.numeric(params$model_index[i])) %>%
    pull(class)

  # tune classifiers
  knn_tune <- alfaknn.tune(x = x, ina = y, a = knn_a_vals, k = knn_k_vals) %>% 
    extract2("ela") %>%
    set_rownames(gsub("alpha=", "", rownames(.))) %>%
    set_colnames(gsub("k=", "", colnames(.))) %>%          
    reshape2::melt(varnames = c("alpha", "k"), value.name = "train_error") %>%
    mutate(train_error = 1 - train_error) %>% # alfaknn.tune returns classification success rate
    group_by(alpha) %>%
    slice_min(train_error, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(classifier = "knn") %>%
    relocate(classifier)
  
  svm_tune <- alfasvm.tune(y = y, x = x, a = svm_a_vals) %>%
    extract2("per") %>%
    as_tibble() %>%
    mutate(alpha = svm_a_vals) %>%
    mutate(classifier = "svm") %>%
    relocate(classifier, alpha) %>%
    rename(train_error = performance) %>%
    mutate(train_error = 1 - train_error) # alfasvm.tune returns classification success rate

  rf_tune <- alfa.rf(xnew = x, y = y, x = x, a = rf_a_vals) %>%
    lapply(function(tmp) {
      config <- tmp$config
      train_error <- tmp$mod %>% 
        lapply(function(item) extract(item, "prediction.error")) %>%
        unlist() %>%
        as.numeric()
      out <- cbind(config, train_error) %>% as_tibble()
      return(out)
    }) %>%
    bind_rows(.id = "alpha") %>%
    mutate(alpha = as.numeric(gsub("alpha=", "", alpha))) %>%
    mutate(classifier = "rf") %>%
    relocate(classifier, alpha) %>%
    group_by(alpha) %>%
    slice_min(train_error, with_ties = FALSE) %>%
    ungroup()
  
  knn_tune$test_error <- sapply(seq_len(nrow(knn_tune)), function(i) {
    alfa.knn(xnew = xnew, x = x, ina = y, 
             a = knn_tune$alpha[i], 
             k = knn_tune$k[i]) %>%
      ModelMetrics::ce(ynew, as.factor(.))
  })
  
  svm_tune$test_error <- sapply(seq_len(nrow(svm_tune)), function(i) {
    alfa.svm(xnew = xnew, y = y, x = x,
             a = svm_tune$alpha[i],
             gamma = svm_tune$gamma[i],
             cost = svm_tune$cost[i]) %>%
      extract2(1) %>%
      extract2("est") %>%
      ModelMetrics::ce(ynew, as.factor(.))
  })
  
  rf_tune$test_error <- sapply(seq_len(nrow(rf_tune)), function(i) {
    alfa.rf(xnew = xnew, y = y, x = x,
            a = rf_tune$alpha[i],
            size = rf_tune$size[i],
            depth = rf_tune$depth[i],
            splits = rf_tune$splits[i],
            R = rf_tune$R[i]) %>%
      extract2(1) %>%
      extract2("est") %>%
      ModelMetrics::ce(ynew, as.factor(.))
  })
  
  tune <- bind_rows(knn_tune, svm_tune, rf_tune) %>%
    mutate(model_index = params$model_index[i], 
           rep = params$rep[i])
  
}) %>%
  bind_rows() %>%
  full_join(params, by = c("model_index", "rep")) %>%
  saveRDS(file = file.path("scripts", "compositional", "results", "classification-sim-results-knn-svm-rf-2.RDS"))


# Ternary plots ---------------------------------------------------------------------

# generate data from each model to display in ternary plots

set.seed(1)

lapply(seq_len(nrow(model_params)), function(i) {
  sim_X_classify(n = 5 * 10^3,
                 d = model_params$d[i],
                 param0 = model_params$param0[i],
                 param1 = model_params$param1[i],
                 model = model_params$model[i],
                 type = "rmev",
                 k_frac = 0.05)
}) %>%
  bind_rows(.id = "model_index") %>%
  full_join(model_params, by = "model_index") %>%
  saveRDS(file = file.path("scripts", "compositional", "results", "classification-sim-ternary-data.RDS"))



# Hyperparameter tuning examples ----------------------------------------------------

# plot loss against (gamma, cost) for SVM with fixed alpha
# plot loss against (k, alpha) for KNN

# generate training data
i <- 2
train_data <- sim_X_classify(n = 10000, 
                             d = model_params$d[i],
                             param0 = model_params$param0[i],
                             param1 = model_params$param1[i],
                             model = model_params$model[i],
                             type = "rmev",
                             k_frac = 0.05)

x <- train_data %>% select(starts_with("X"))
y <- train_data %>% pull(class)

knn_tune <- alfaknn.tune(x = x, ina = y, a = knn_a_vals, k = knn_k_vals) %>% 
  extract2("ela") %>%
  set_rownames(gsub("alpha=", "", rownames(.))) %>%
  set_colnames(gsub("k=", "", colnames(.)))

n_colpal <- 13
colpal_seq <- rev(colorspace::sequential_hcl(n = n_colpal, "Viridis"))
knn_tune <- 1 - knn_tune
lattice::levelplot(knn_tune,
                   col.regions = colpal_seq,
                   at = lattice::do.breaks(c(min(knn_tune), max(knn_tune)), length(colpal_seq)),
                   xlab = expression("Transformation parameter," ~ alpha), ylab = expression("Neighbours," ~ k),
                   colorkey = TRUE,
                   scales = list(x = list(at = seq_along(knn_a_vals),
                                          labels = knn_a_vals),
                                 y = list(at = seq_along(knn_k_vals),
                                          labels = knn_k_vals))) %>%
  update(aspect = 0.8)

rf_tune <- alfa.rf(xnew = x, y = y, x = x, a = 0) %>%
  lapply(function(tmp) {
    config <- tmp$config
    train_error <- tmp$mod %>% 
      lapply(function(item) extract(item, "prediction.error")) %>%
      unlist() %>%
      as.numeric()
    out <- cbind(config, train_error) %>% as_tibble()
    return(out)
  }) %>%
  bind_rows(.id = "alpha") %>%
  mutate(alpha = as.numeric(gsub("alpha=", "", alpha)))

ggplot(rf_tune, aes(x = size, y = train_error, colour = as.factor(splits))) + 
  geom_line() +
  facet_grid(as.factor(depth) ~ ., scales = "free_y") +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(labels = label_percent()) +
  scale_colour_manual(values = c("red", "blue", "darkgreen", "black")) +
  labs(colour = "Splits", x = "Size", y = "Empirical risk") +
  theme_light() +
  theme(legend.position = "top", 
        panel.grid.major = element_blank()) 

