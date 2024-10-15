sl_tpdm <- function(a) { 
  sl_integrand <- function(x, a) (1 - a) / a * (x * (1-x))^(1/a - 3/2) * ((1-x)^(1/a) + x^(1/a))^(a-2) 
  if (a == 0) { 
    1 
  } else if (a == 1) { 
    0 
  } else { 
    integrate(function(x) sl_integrand(x, a = a), lower = 0, upper = 1, subdivisions = 500)$value 
  } 
} 
sl_tpdm <- Vectorize(sl_tpdm)

hr_tpdm <- function(lambda) { 
  hr_integrand <- function(x, a) { 
    g1 <- a/2 - 1/a * log(x^2/(1-x^2)) 
    g2 <- a/2 - 1/a * log((1-x^2)/x^2) 
    term1 <- 2 / (a * sqrt(1 - x^2)) 
    term2 <- ((1 - g1/a) / x^2) * dnorm(g1) 
    term3 <- ((1 - g2/a) / (1-x^2)) * dnorm(g2) 
    out <- term1 * (term2 + term3) 
    return(out) 
  } 
  ifelse(lambda == 0, 
         1, 
         integrate(function(x) hr_integrand(x, a = 2 * lambda), lower = 0, upper = 1, 
                   subdivisions = 500)$value)   
  
} 

hr_tpdm <- Vectorize(hr_tpdm)

sl_chi <- function(x) {2 - 2^x} 
hr_chi <- function(lambda) {2 * (1 - pnorm(lambda))}