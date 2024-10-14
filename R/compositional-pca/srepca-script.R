# source("Code/srepca-functions.R")
# 
# options(warn = -1, dplyr.summarise.inform = FALSE)
# 
# set.seed(1)
# block.d.vals <- c(3, 3) 
# block.m.vals <- c(2, 2)
# base.A <- create.block.A(block.d.vals, block.m.vals)
# n.ext <- 100
# n.vals <- c(1000)
# nrep <- 20
# sigma.noise.vals <- c(0, 0.6)
# lambda.vals <- 10^seq(from = -5, to = 0, by = 0.25)
# pen.types <- c("test", "p.value")
# 
# all.data <- expand_grid("sigma.noise" = sigma.noise.vals, "n" = n.vals, "rep" = 1:nrep) %>%
#   mutate(A = list(base.A)) %>%
#   relocate(A, .before = sigma.noise) %>%
#   rowwise() %>%
#   mutate(snr = tlfm.snr(A = A, sigma = sigma.noise)) %>%
#   mutate(data = list(sim.data.A(A = A, n = n, sigma = sigma.noise, univar.qthresh = 1 - n.ext / n))) %>%
#   ungroup()
# 
# tpdm.est.params <- expand_grid(method = "lasso", pen.type = pen.types) %>%
#   full_join(tibble(method = "standard"))
# 
# test.params <- expand_grid(all.data, tpdm.est.params) %>%
#   rowwise() %>%
#   mutate(Sigma.true = list(tlfm.tpdm(A)), n.ext = n.ext, rt.qthresh = 1 - n.ext / n) %>%
#   relocate(Sigma.true, .after = A) %>% 
#   ungroup()
# 
# est.results <- test.params %>%
#   rowwise() %>%
#   mutate(tpdm.est = list(tryCatch(
#     {tpdm.est(data = data, rt.qthresh = rt.qthresh, method = method, lambda = lambda.vals, pen.type = pen.type)},
#     error = function(e) {}))) %>%
#   unnest(tpdm.est) %>%
#   rename(Sigma.hat = Sigma) %>%
#   rowwise() %>%
#   mutate(Sigma.error = tpdm.error(tpdm.est = Sigma.hat, tpdm.true = Sigma.true)) %>%
#   mutate(Sigma.hat.eigen = list(tryCatch({
#     Sigma.hat %>% Sigma.dat.to.mat(value.var = "sigma.hat.ij") %>% tpdm.eigen()
#   }, error = function(e) {}))) %>%
#   mutate(Sigma.hat.evals = list(Sigma.hat.eigen$values)) %>%
#   mutate(Sigma.hat.evecs = list(Sigma.hat.eigen$vectors)) %>%
#   mutate(Sigma.true.eigen = list(tryCatch({
#     Sigma.true %>% Sigma.dat.to.mat(value.var = "sigma.ij") %>% tpdm.eigen()
#   }, error = function(e) {}))) %>%
#   mutate(Sigma.true.evals = list(Sigma.true.eigen$values)) %>%
#   mutate(Sigma.true.evecs = list(Sigma.true.eigen$vectors)) %>%
#   ungroup() %>%
#   dplyr::select(-Sigma.hat.eigen) %>%
#   dplyr::select(-Sigma.true.eigen)
# 
# saveRDS(est.results, "Code/tpdm.est.results.rds")
# 
# base.A
# base.A %*% t(base.A)
