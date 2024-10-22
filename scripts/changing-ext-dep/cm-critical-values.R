rm(list = ls())

set.seed(1)

delta_t <- 1 / 400
nrep <- 5 * 10^6

bb_L2 <- pbreplicate(n = nrep, e1071::rbridge(end = 1, frequency = 1 / delta_t)) %>%
  apply(2, function(z) sum(z^2) * delta_t)

saveRDS(bb_L2, file = file.path("scripts", "changing-ext-dep", "results", "bb_L2.RDS"))