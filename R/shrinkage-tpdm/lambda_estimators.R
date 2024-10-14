lw_ratio_method <- function(X, k_fun, kappa_vals, B = 5, sub_ind = NULL, lambda_min = 0, lambda_max = 1) {
  
  if (lambda_max == 0) {return(data.frame("n" = n,
                                          "k" = k,
                                          "d" = d,
                                          "B" = B,
                                          "sigma_hat_ij_mean" = NA,
                                          "E_sigma_hat_norm" = NA,
                                          "kappa" = NA,
                                          "k_kappa" = NA,
                                          "E_sigma_hat_small_norm" = NA,
                                          "zeta_hat" = NA,
                                          "lambda_hat" = 0))}
  
  d <- ncol(X)
  n <- nrow(X)
  
  # compute empirical TPDM
  k <- k_fun(n)
  Sigma_hat <- tpdm_estimate(X, k)
  
  # bootstrap data of size n and compute empirical TPDM
  Sigma_hat_boot <- lapply(seq_len(B), function(b) {
    S <- X[sample(1:n, size = n, replace = TRUE), ] %>% tpdm_estimate(k = k)
    if (!is.null(sub_ind)) {S <- S[sub_ind[, 1], sub_ind[, 2]]}
    return(S)
  })
  
  # compute average bootstrapped empirical TPDM
  Sigma_hat_boot_mean <- Reduce("+", Sigma_hat_boot) / length(Sigma_hat_boot)
  
  # run the rest of the method for each kappa
  tmp <- lapply(kappa_vals, function(kappa) {
    
    k_kappa <- k_fun(floor(kappa * n))
    
    # bootstrap data of size kappa * n and compute empirical TPDM
    Sigma_hat_small_boot <- lapply(seq_len(B), function(b) {
      S <- X[sample(1:n, size = floor(kappa * n), replace = TRUE), ] %>% tpdm_estimate(k = k_kappa)
      if (!is.null(sub_ind)) {S <- S[sub_ind[, 1], sub_ind[, 2]]}
      return(S)
    })
    
    # compute average bootstrapped empirical TPDM
    Sigma_hat_small_boot_mean <- Reduce("+", Sigma_hat_small_boot) / length(Sigma_hat_small_boot)
    
    # estimate ratio zeta
    zeta_hat_num <- k * sum(vecu(Sigma_hat_boot_mean)^2)
    zeta_hat_den <- k_kappa * sum(vecu(Sigma_hat_small_boot_mean)^2)
    zeta_hat <- zeta_hat_num / zeta_hat_den
    
    # compute lambda_hat
    if (lambda_min >= min(1, 1 / zeta_hat, lambda_max)) {
      lambda_hat <- lambda_min
    } else {
      lambda_hat <- optim(par = lambda_min,
                          fn = function(x) frob_loss(shrink_cov(Sigma_hat_boot_mean, "lw", x),
                                                     shrink_cov(Sigma_hat_small_boot_mean, "lw", x * zeta_hat)),
                          lower = lambda_min, upper = min(1, 1 / zeta_hat, lambda_max), method = "L-BFGS-B")$par
    }
    
    return(data.frame("n" = n,
                      "k" = k,
                      "d" = d,
                      "B" = B,
                      "sigma_hat_ij_mean" = mean(vecu(Sigma_hat)),
                      "E_sigma_hat_norm" = sum(vecu(Sigma_hat_boot_mean)^2),
                      "kappa" = kappa,
                      "k_kappa" = k_kappa,
                      "E_sigma_hat_small_norm" = sum(vecu(Sigma_hat_small_boot_mean)^2),
                      "zeta_hat" = zeta_hat,
                      "lambda_hat" = lambda_hat))
    
  }) %>% bind_rows()
  return(tmp)
  
}


sl_oracle_lambda <- function(Sigma_hat, gamma, shrinkage_type, eta = 1) {
  
  sigma_bar <- mean(vecu(Sigma_hat))
  sigma <- sl_tpdm(gamma)
  if (sigma_bar < sigma) return(0)
  
  if (shrinkage_type == "hard") {
    optimal_lambda <- ifelse(abs(sigma_bar - sigma) > sigma, 1, 0)
  } else if (shrinkage_type == "soft") {
    optimal_lambda <- sigma_bar - sigma
    optimal_lambda <- max(0, optimal_lambda)
  } else if (shrinkage_type == "adlasso") {
    optimal_lambda <- (sigma_bar^eta * (sigma_bar - sigma))^(1 / (eta + 1))
    optimal_lambda <- max(0, optimal_lambda)
  } else if (shrinkage_type == "lw") {
    optimal_lambda <- 1 - sigma / sigma_bar
    optimal_lambda <- min(1, max(optimal_lambda, 0))
  } else {
    stop("Input a valid shrinkage type.")
  }
  return(optimal_lambda)
}