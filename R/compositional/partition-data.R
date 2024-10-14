extreme_data_partition <- function(data, train_frac, k_frac) {
  
  # response variable must be called "class"
  # assumes all other variables are the predictors
  
  train_index <- createDataPartition(data$class, p = train_frac, list = FALSE)
  train_data <- data[train_index, ]
  test_data <- data[-train_index, ]
  
  train_ext_data <- train_data %>%
    mutate(R = rowSums(select(., -class))) %>%
    slice_max(R, prop = k_frac)
  
  R_thresh <- min(train_ext_data$R)
  
  train_ext_data <- train_ext_data %>%
    mutate(across(.cols = -c("class", "R"),
                  .fns = ~ .x / R,
                  .names = "theta_{.col}")) %>%
    select(-R)
  
  test_ext_data <- test_data %>%
    mutate(R = rowSums(select(., -class))) %>%
    filter(R >= R_thresh) %>%
    mutate(across(.cols = -c("class", "R"),
                  .fns = ~ .x / R)) %>%
    rename_with(.cols = -c("class", "R"),
                .fn = ~ paste("theta", .x, sep = "_")) %>%
    arrange(desc(R)) %>%
    select(-R)
  
  return(list("train_ext_data" = train_ext_data,
              "test_ext_data" = test_ext_data,
              "threshold" = R_thresh))
}