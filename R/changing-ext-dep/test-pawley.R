vecu <- function(M) {
  m <- t(M)[lower.tri(t(M))]
  names(m) <- gtools::combinations(n = ncol(M), r = 2, v = colnames(M)) %>% 
    apply(1, str_flatten, collapse = "_")
  return(m)
}

test_pawley <- function(X, b, k, return_all = FALSE) {
  
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
  block_sigma <- block_angles %>%
    lapply(function(Theta_ext) d / k * t(Theta_ext) %*% Theta_ext) %>%
    lapply(vecu) %>%
    do.call(rbind, .) %>%
    set_colnames(gsub("theta_", "", colnames(.)))
  
  block_psi <- block_sigma %>%
    apply(2, function(x) cumsum(x / (n / b)))
  
  if (d > 2) {
    V_sqrtinv <- block_angles %>%
      do.call(rbind, .) %>%
      apply(1, function(theta) d * combn(theta, 2, prod)) %>%
      t() %>%
      cov() %>%
      pracma::inv() %>% # chol %>% chol2inv
      expm::sqrtm()
  } else {
    V_sqrtinv <- block_angles %>%
      lapply(function(W) d * apply(W, 1, prod)) %>%
      lapply(var) %>%
      unlist() %>%
      mean() %>% 
      raise_to_power(-1/2)
  }
  
  t <- seq(from = b / n, to = 1, by = b / n)
  psi_1 <- matrix(tail(block_psi, n = 1), byrow = TRUE, nrow = nrow(block_psi), ncol = ncol(block_psi))
  Z <- sqrt(k * n / b) * V_sqrtinv %*% t(block_psi - t * psi_1) %>% t()
  colnames(Z) <- colnames(block_psi)
  
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
      melt(block_sigma, varnames = c("block", "variables"), value.name = "sigma"),
      full_join(
        melt(block_psi, varnames = c("block", "variables"), value.name = "psi"),
        melt(Z, varnames = c("block", "variables"), value.name = "z"),
        by = c("block", "variables")),
      by = c("block", "variables")) %>%
      arrange(block) %>%
      group_by(variables) %>%
      mutate(ks = cummax(abs(z)),
             cm = (b / n) * cumsum(z^2)) %>%
      ungroup() %>%
      mutate(variables = gsub("_", " & ", variables))
    return(list("data" = data, "ks" = ks, "cm" = cm, "elapsed_time" = elapsed_time))
  }
}

