sl_tpdm <- function(gamma) {
  sl_integrand <- function(x, gamma)
    (1 - gamma) / gamma * (x * (1 - x)) ^ (1 / gamma - 3 / 2) * ((1 - x) ^ (1 / gamma) + x ^
                                                                   (1 / gamma)) ^ (gamma - 2)
  if (gamma == 0) {
    1
  } else if (gamma == 1) {
    0
  } else {
    integrate(
      function(x)
        sl_integrand(x, gamma = gamma),
      lower = 0,
      upper = 1,
      subdivisions = 500
    )$value
  }
} 

hr_tpdm <- function(lambda) {
  hr_integrand <- function(x, a) {
    g1 <- a / 2 - 1 / a * log(x ^ 2 / (1 - x ^ 2))
    g2 <- a / 2 - 1 / a * log((1 - x ^ 2) / x ^ 2)
    term1 <- 2 / (a * sqrt(1 - x ^ 2))
    term2 <- ((1 - g1 / a) / x ^ 2) * dnorm(g1)
    term3 <- ((1 - g2 / a) / (1 - x ^ 2)) * dnorm(g2)
    out <- term1 * (term2 + term3)
    return(out)
  }
  ifelse(
    lambda == 0,
    1,
    integrate(
      function(x)
        hr_integrand(x, a = 2 * lambda),
      lower = 0,
      upper = 1,
      subdivisions = 500
    )$value
  )
  
}

hr_tpdm <- Vectorize(hr_tpdm)
sl_tpdm <- Vectorize(sl_tpdm)

make_hr_variogram <- function(d, beta = 3) {
  H <- matrix(data = EnvStats::rpareto(d^2, location = 1, shape = 2.5), nrow = d, ncol = d)
  Gamma <- (beta / d) * dist(H, upper = TRUE, diag = TRUE)^2 %>% as.matrix()
  return(Gamma)
}
