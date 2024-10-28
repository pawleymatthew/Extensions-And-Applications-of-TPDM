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

