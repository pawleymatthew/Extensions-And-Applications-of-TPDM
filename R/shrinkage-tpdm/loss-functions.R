vecu <- function(M, add_names = FALSE) {
  m <- t(M)[lower.tri(t(M))]
  if (add_names) {
    names(m) <- gtools::combinations(n = ncol(M), r = 2, v = seq_len(ncol(M))) %>% 
      apply(1, str_flatten, collapse = "_")
  }
  return(m)
}

frob_loss <- function(Sigma, Sigma_tilde) {
  sum((vecu(Sigma) - vecu(Sigma_tilde))^2)
}

krzanowski_measure <- function(Sigma, Sigma_tilde) {
  
  Sigma_tilde <- nearPD(Sigma_tilde, keepDiag = TRUE)$mat %>% as.matrix()
  U_tilde <- eigen(Sigma_tilde)$vectors
  U <- eigen(Sigma)$vectors
  
  K_p <- sapply(1:d, function(p) {
    K <- t(U_tilde[, 1:p, drop = FALSE]) %*% U[, 1:p, drop = FALSE] 
    return(sum(K^2))
  })
  
  return(K_p)
  
}