rm(list = ls())

sapply(list.files(path = "R/changing-ext-dep", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)

# Load critical values and data ----------------------------------------------------------------------

bb_L2 <- readRDS(file = file.path("scripts", "changing-ext-dep", "results", "bb_L2.RDS"))

ks_cv <- sapply(1:100, function(d) CPAT:::qkolmogorov(0.95^(2/(d*(d-1)))))
cm_cv <- sapply(1:100, function(d) quantile(bb_L2, probs = 0.95^(2/(d*(d-1))))) %>% as.numeric()

redsea <- load_red_sea_temp(alpha = 2)


# Resample sites and run tests ------------------------------------------------------

set.seed(1)

d_vals <- c(2, 5, 10)

lapply(seq_along(d_vals), function(d_ind) {
  
  d <- d_vals[d_ind]
  nreps <- 1000
  
  pbreplicate(n = nreps, expr = {
    
    north_sites <- redsea$coord %>%
      mutate(loc_index = row_number()) %>%
      filter(region == "north") %>%
      slice_sample(n = d) %>%
      pull(loc_index)
    
    redsea$X[1:1605, north_sites] %>%
      as_tibble() %>%
      test_pawley(b = 107, k = 20, return_all = FALSE)
    
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
  bind_rows() %>%
  saveRDS(file = file.path("scripts", "changing-ext-dep", "results", "redsea-north-p-values.RDS"))


lapply(seq_along(d_vals), function(d_ind) {
  
  d <- d_vals[d_ind]
  nreps <- 1000
  
  pbreplicate(n = nreps, expr = {
    
    south_sites <- redsea$coord %>%
      mutate(loc_index = row_number()) %>%
      filter(region == "south") %>%
      slice_sample(n = d) %>%
      pull(loc_index)
    
    redsea$X[1:1605, south_sites] %>%
      as_tibble() %>%
      test_pawley(b = 107, k = 20, return_all = FALSE)
    
  }) %>%
    t() %>%
    as_tibble() %>%
    mutate(d = d, region = "south") %>%
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
  bind_rows() %>%
  saveRDS(file = file.path("scripts", "changing-ext-dep", "results", "redsea-south-p-values.RDS"))


