s_hard <- function(x, lambda) x * (abs(x) > lambda)

s_soft <- function(x, lambda) sign(x) * pmax(abs(x) - lambda, 0)

s_adlasso <- function(x, lambda, eta = 1) sign(x) * pmax(abs(x) - lambda^(eta + 1) / abs(x)^eta, 0)

shrink_cov <- function(M, shrinkage_type, lambda, shrink_diag = FALSE, ...) {
  
  if (shrinkage_type == "hard") {
    sM <- s_hard(M, lambda = lambda)
  } else if (shrinkage_type == "soft") {
    sM <- s_soft(M, lambda = lambda)
  } else if (shrinkage_type == "adlasso") {
    sM <- s_adlasso(M, lambda = lambda, ...)
  } else if (shrinkage_type == "lw") {
    if (lambda < 0 | lambda > 1) stop("LW shrinkage intensity must be between 0 and 1.")
    sM <- lambda * diag(ncol(M)) + (1 - lambda) * M
    shrink_diag <- FALSE
  } else {
    stop("Input a valid shrinkage type.")
  }
  if(!shrink_diag) diag(sM) <- diag(M)
  
  return(sM)
}