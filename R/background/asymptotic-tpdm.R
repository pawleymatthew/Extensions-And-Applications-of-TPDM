library(gtools)
library(pracma)

vecu <- function(M, add_names = FALSE) {
  m <- t(M)[lower.tri(t(M))]
  if (add_names) {
    names(m) <- gtools::combinations(n = ncol(M), r = 2, v = seq_len(ncol(M))) %>% 
      apply(1, str_flatten, collapse = "_")
  }
  return(m)
}

# asymptotic covariance matrix V
asy_V <- function(A, alpha = 1) {
  d <- nrow(A)
  q <- ncol(A)
  B <- matrix(NA, nrow = choose(d, 2), ncol = q)
  ij_pairs <- gtools::combinations(n = d, r = 2, repeats.allowed = FALSE) %>% set_colnames(c("i", "j"))
  for (ind in seq_len(nrow(ij_pairs))) {
    ind_i <- ij_pairs[ind, "i"]
    ind_j <- ij_pairs[ind, "j"]
    for (s in 1:q) {
      B[ind, s] <- A[ind_i, s]^(alpha / 2) * A[ind_j, s]^(alpha / 2) / pracma::Norm(A[, s], p = alpha)^(alpha / 2)
    }
  }
  V <- matrix(NA, nrow = choose(d, 2), ncol = choose(d, 2))
  for (ind1 in seq_len(nrow(ij_pairs))) {
    for (ind2 in seq_len(nrow(ij_pairs))) {
      ind_i <- ij_pairs[ind1, "i"]
      ind_j <- ij_pairs[ind1, "j"]
      ind_l <- ij_pairs[ind2, "i"]
      ind_m <- ij_pairs[ind2, "j"]
      V[ind1, ind2] <- d * sum(B[ind1, ] * B[ind2, ]) - Sigma[ind_i, ind_j] * Sigma[ind_l, ind_m]
    }
  }
  
  rownames(V) <- paste0("hat(sigma)[", ij_pairs[, "i"], ij_pairs[, "j"], "]")
  colnames(V) <- rownames(V)
  
  return(V)
}
