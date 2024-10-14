test_drees <- function(X, b, k, return_all = FALSE, y_vals = seq(from = 0.01, to = 0.99, by = 0.01)) {
  
  tic()
  
  n <- nrow(X)
  d <- ncol(X)
  k <- floor(k) # ensures k is an integer
  
  block_angles <- X %>% 
    mutate(block = ceiling(row_number() / b)) %>%
    relocate(block) %>%
    mutate(R = sqrt(rowSums(select(., -block)^2))) %>%
    group_by(block) %>%
    slice_max(R, n = k) %>%
    mutate(across(
      .cols = !R,
      .fns = ~ .x / R,
      .names = "theta_{.col}")) %>%
    group_by(block) %>%
    group_split(.keep = FALSE) %>%
    lapply(function(tmp) as.matrix(select(tmp, starts_with("theta"))))
  
  # compute integrated TPDM at (end) of each block
  block_Sy <- block_angles %>%
    lapply(function(Theta_ext) sapply(y_vals, function(y)  mean(Theta_ext[, 1] < y))) %>%
    do.call(rbind, .) %>%
    set_colnames(paste0("y_", seq_along(y_vals)))
  
  block_ISy <- block_Sy %>%
    apply(2, function(x) cumsum(x / (n / b)))
  
  t <- seq(from = b / n, to = 1, by = b / n)
  ISy_1 <- matrix(tail(block_ISy, n = 1), byrow = TRUE, nrow = nrow(block_ISy), ncol = ncol(block_ISy))
  Z <- sqrt(k * n / b) * t(block_ISy - t * ISy_1) %>% t()
  colnames(Z) <- colnames(block_ISy)
  
  # test statistics
  ks <- max(abs(Z))
  cm <- apply(Z, 2, function(z) (b / n) * sum(z^2)) %>% max()
  
  elapsed_time <- toc(quiet = TRUE)
  elapsed_time <- (elapsed_time$toc - elapsed_time$tic) %>% as.numeric()
  
  # return objects
  if (!return_all) {
    return(list("ks" = ks, "cm" = cm, "elapsed_time" = elapsed_time))
  } else {
    data <- full_join(
      melt(block_Sy, varnames = c("block", "set"), value.name = "S"),
      full_join(
        melt(block_ISy, varnames = c("block", "set"), value.name = "IS"),
        melt(Z, varnames = c("block", "set"), value.name = "z"),
        by = c("block", "set")),
      by = c("block", "set")) %>%
      arrange(block) %>%
      group_by(set) %>%
      mutate(ks = cummax(abs(z)),
             cm = (b / n) * cumsum(z^2)) %>%
      ungroup()
    return(list("data" = data, "ks" = ks, "cm" = cm, "elapsed_time" = elapsed_time))
  }
}
