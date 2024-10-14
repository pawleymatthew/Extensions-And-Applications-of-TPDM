ppolar <- function(X, p = 2) {
  R <- apply(X, 1, pracma::Norm, p = 2)
  Theta <- X / R
  order_stat <- order(R, decreasing = TRUE)
  return(list("R" = R, "Theta" = Theta, "order_stat" = order_stat))
}

A_hat <- function(X, k) {
  m <- ncol(X)
  pX <- ppolar(X, p = 2)
  Ak_hat <- lapply(k, function(k) {
    A_k <- sqrt(m/k) * t(pX$Theta[pX$order_stat[1:k], ])
    rownames(A_k) <- colnames(X)
    return(A_k)
  })
  names(Ak_hat) <- paste("k", k, sep = "=")
  return(Ak_hat)
} 

tpdm_estimate <- function(X, k) {
  Sigma_hat <- lapply(A_hat(X, k), function(A) A %*% t(A))
  if (length(k) == 1) Sigma_hat <- Sigma_hat[[1]]
  return(Sigma_hat)
}
