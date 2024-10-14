
# sim_X_classify <- function(n, d, param0, param1, model, p = 0.5) {
#   
#   # "log" : param should be in (0, 1). 
#   # "neglog" : param must be positive
#   # "bilog" : param must be in [0, 1]
#   
#   n0 <- floor(n * p)
#   n1 <- n - n0
#   
#   if (model == "log") {
#     X0 <- mev::rmev(n = n0, d = d, param = 1 / param0, model = model) 
#     X1 <- mev::rmev(n = n1, d = d, param = 1 / param1, model = model) 
#   } else if (model == "neglog") {
#     X0 <- mev::rmev(n = n0, d = d, param = param0, model = model) 
#     X1 <- mev::rmev(n = n1, d = d, param = param1, model = model)  
#   } 
#   else if (model == "bilog") {
#     X0 <- mev::rmev(n = n0, d = d, param = rep(param0 / seq_len(ceiling(d / 2)), length.out = d), model = model) 
#     X1 <- mev::rmev(n = n1, d = d, param = rep(param1 / seq_len(ceiling(d / 2)), length.out = d), model = model) 
#   } else {
#     stop("Specify valid dependence model.")
#   }
#   
#   X0 <- X0 %>%
#     set_colnames(paste0("X", seq_len(ncol(.)))) %>%
#     as_tibble() %>%
#     mutate(class = "1")
#   
#   X1 <- X1 %>%
#     set_colnames(paste0("X", seq_len(ncol(.)))) %>%
#     as_tibble() %>%
#     mutate(class = "2")
#   
#   X <- bind_rows(X0, X1) %>%
#     mutate(class = as.factor(class)) %>%
#     sample_frac(1)
#   
#   return(X)
# }


sim_X_classify <- function(n, d, param0, param1, model, type = "rmev", k_frac = 0.05, p = 0.5) {
  
  # "log" : param should be in (0, 1). 
  # "neglog" : param must be positive
  # "bilog" : param must be in [0, 1]
  
  # map parameters to inputs for mev::rmev / mev::rmevspec
  if (model == "log") {
    param0 <- 1 / param0
    param1 <- 1 / param1
  } else if (model == "bilog") {
    param0 <- rep(param0 / seq_len(ceiling(d / 2)), length.out = d)
    param1 <- rep(param1 / seq_len(ceiling(d / 2)), length.out = d)
  }
  
  n0 <- floor(n * p)
  n1 <- n - n0
  
  if (type == "rmev") {
    X0 <- mev::rmev(n = n0, d = d, param = param0, model = model)
    X1 <- mev::rmev(n = n1, d = d, param = param1, model = model)
  } else {
    X0 <- mev::rmevspec(n = n0, d = d, param = param0, model = model)
    X1 <- mev::rmevspec(n = n1, d = d, param = param1, model = model)
  }
  
  X0 <- X0 %>% 
    set_colnames(paste0("X", seq_len(ncol(.)))) %>%
    as_tibble() %>% 
    mutate(class = "1")
  X1 <- X1 %>% 
    set_colnames(paste0("X", seq_len(ncol(.)))) %>%
    as_tibble() %>% 
    mutate(class = "2")
  
  X <- bind_rows(X0, X1) %>%
    mutate(class = as.factor(class)) %>%
    sample_frac(1)
  
  if (type == "rmev") {
    X <- X %>%
      mutate(R = rowSums(select(., -class))) %>%
      slice_max(R, prop = k_frac) %>%
      mutate(across(.cols = starts_with("X"),
                    .fns = ~ .x / R)) %>%
      select(-R)
  }
  
  return(X)
}




