# Marginal estimation and transformation --------------------------------------------

# Hill's estimator of tail index
# if X is a dxn data matrix: `apply(X, 2, hill_estimate, k = 250)`
hill_estimate <- function(X, k){
  apply(X, 2, function(x) {
    x_order <- sort(x, decreasing = TRUE)
    alpha <- sum(log(x_order[1:k])) / k - log(x_order[k]) %>% as.numeric()
    return(1 / alpha)
  })
}

margins_to_frechet <- function(X, alpha = 1, type = "empirical") {
  if (type == "empirical") {
    Y <- ExtremalDep::trans2UFrechet(data = X, type = "Empirical")
  } else if (type == "GEV") {
    pars <- apply(X, 2, function(x) evd::fgev(x, std.err = FALSE)$estimate) %>% t()
    Y <- ExtremalDep::trans2UFrechet(data = X, pars = pars, type = "GEV")
  } else {
    stop('type must be "empirical" or "GEV".')
  }
  Y <- Y^(1 / alpha)
  dimnames(Y) <- dimnames(X)
  return(Y)
}
