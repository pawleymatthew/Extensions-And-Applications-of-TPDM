# MEV copulae -----------------------------------------------------------------------

rgumbelcop <- function(n, lambda, d = 2, alpha = 2) {
  cop <- SimCop::NewMEVGumbelCopula(r = lambda)
  X <- cop %>%
    SimCop::GetApprox(dim = d) %>%
    SimCop::GenerateRV(n = n) %>%
    set_colnames(paste0("X", 1:ncol(.))) %>%
    margins_to_frechet(alpha = alpha)
  return(list("X" = X,
              "n" = n,
              "d" = d,
              "description" = glue::glue("Gumbel copula with dependence parameter {lambda} on Fréchet({alpha}) margins. Simulated using SimCop package.")))
}

# # simulate ndata datasets of size n with time-varying lambda
# rgumbelcop_timedep <- function(ndata, n, lambda_fun, ...) {
#   t_vals <- (1:n) / n
#   lambda_vals <- sapply(t_vals, lambda_fun)
#   data <- pblapply(lambda_vals, function(lambda_t) rgumbelcop(n = ndata, lambda = lambda_t, ...)$X)
#   lapply(data, function(Xt) {
#     Xt %>%
#       as.data.frame() %>%
#       group_by(ind = row_number(), .add = FALSE)
#   }) %>%
#     bind_rows() %>%
#     group_split(.keep = FALSE) %>%
#     lapply(as.matrix) %>%
#     set_names()
# }

