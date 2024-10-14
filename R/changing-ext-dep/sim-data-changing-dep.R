make_hr_variogram <- function(d, beta = 3) {
  H <- matrix(data = EnvStats::rpareto(d^2, location = 1, shape = 2.5), nrow = d, ncol = d)
  Gamma <- (beta / d) * dist(H, upper = TRUE, diag = TRUE)^2 %>% as.matrix()
  return(Gamma)
}

sim_X_changing_dep <- function(n, d, model, param0, param1, change_type = "jump") {
  
  # make base variogram if HR model
  if (model == "hr") {
    Gamma0 <- make_hr_variogram(d = d)
  }
  
  # helper function to simulate n samples from model with given parameter
  sim_X <- function(n, param) {
    if (model == "hr") {
      X <- mev::rmev(n = n, d = d, sigma = param * Gamma0, model = "hr")
    } else {
      X <- mev::rmev(n = n, d = d, param = param, model = model)
    }
    return(as.data.frame(X))
  }
  # simulate data according to change_type
  if (param0 == param1) { # no change in dependence
    X <- sim_X(n = n, param = param0)
  } else if (change_type == "jump") { # jump change in dependence
    n0 <- floor(0.5 * n)
    n1 <- n - n0
    X <- rbind(sim_X(n = n0, param = param0), sim_X(n = n1, param = param1))
  } else if (change_type == "linear") { # linear change in dependence
    X <- lapply((1:n)/n, function(t) {
      X_t <- sim_X(n = 1, param = param0 + t * (param1 - param0))
      }) %>%
      bind_rows()
  } else {
    stop("Specify valid dependence change type.")
  }
  X <- sqrt(X)
  return(X)
}
