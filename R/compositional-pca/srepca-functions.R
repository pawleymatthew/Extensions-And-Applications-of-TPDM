# # General ---------------------------------------------------------------------------
# 
# rowNorms <- function(M) {
#   if (is.vector(M)) M <- matrix(M, nrow = 1)
#   r <- sqrt(rowSums(M^2))
#   return(r)
# }
# 
# colNorms <- function(M) rowNorms(t(M))
# 
# i.set <- function(data) data %>% distinct(i) %>% pull(i)
# 
# t.set <- function(data) data %>% distinct(t) %>% pull(t)
# 
# # list all combinations of two vectors (i.e. pairs up to symmetry, so AB=BA)
# make.i.pairs <- function(i.set) {
#   i.pairs <- gtools::combinations(n = length(i.set), r = 2, repeats.allowed = T, v = i.set) %>%
#     as_tibble() %>%
#     rename(i = V1, j = V2)
#   return(i.pairs)
# }
# 
# # Transformed linear operations -----------------------------------------------------
# 
# tau <- function(x, inverse = FALSE) {
#   # softplus function and its inverse
#   if (!inverse){
#     y <- -plogis(x, lower.tail = FALSE, log.p = TRUE) # log(1 + e^x)
#   } else {
#     y <- qlogis(-x, lower.tail = FALSE, log.p = TRUE) # log(e^x - 1)
#   }
#   return(y)
# }
# 
# vecadd <- function(u, v) {
#   tau(tau(u, inverse = TRUE) + tau(v, inverse = TRUE))
# }
# 
# scalmult <- function(a, v) {
#   tau(a * tau(v, inverse = TRUE))
# }
# 
# matvecmult <- function(M, v) {
#   tau(M %*% tau(v, inverse = TRUE))
# }
# 
# # Transformed linear factor model (TLFM) --------------------------------------------
# 
# create.block.A <- function(block.d.vals, block.m.vals) {
#   A <- lapply(seq_along(block.d.vals), FUN = function(i) {
#     di <- block.d.vals[i]
#     mi <- block.m.vals[i]
#     Ai <- matrix(data = runif(di * mi), nrow = di, ncol = mi)
#     return(Ai)
#   }) %>% 
#     bdiag() %>% 
#     as.matrix()
#   A <- A / rowNorms(A)
#   d <- nrow(A)
#   if (d <= 26) {
#     rownames(A) <- LETTERS[1:d]
#   } else {
#     rownames(A) <- paste0("Loc", 1:d)
#   }
#   return(A)
# }
# 
# add.noise.A <- function(A, sigma.noise) {
#   A <- A + (A == 0) * runif(n = nrow(A) * ncol(A), max = sigma.noise)
#   return(A)
# }
# 
# rtlfm <- function(n, A, sigma = 0, alpha = 2) {
#   d <- nrow(A)
#   m <- ncol(A)
#   X <- matrix(NA, nrow = n, ncol = d)
#   for (i in 1:n) {
#     Z <- evd::rgpd(n = m, loc = 0, scale = 1, shape = 1/alpha)
#     X[i, ] <- matvecmult(A, Z)
#     if (sigma > 0) {
#       N <- TruncatedNormal::rtnorm(n = d, mu = 0, sd = 1, lb = 0, ub = Inf)
#       eta <- evd::rgpd(n = 1, loc = 0, scale = 1, shape = 1/alpha)
#       epsilon <- scalmult(a = eta, v = N)
#       X[i, ] <- vecadd(u = X[i, ], v = scalmult(a = sigma, v = epsilon))
#     } 
#   }
#   return(X)
# }
# 
# tlfm.tpdm <- function(A) {
#   i.set <- rownames(A)
#   i.pairs <- make.i.pairs(i.set)
#   Sigma.true <- A %*% t(A)
#   Sigma <- i.pairs %>%
#     rowwise() %>%
#     mutate(sigma.ij = Sigma.true[i, j]) %>%
#     ungroup()
#   return(Sigma)
# }
# 
# tlfm.snr <- function(A, sigma) {
#   if (sigma == 0) {
#     return(1)
#   } else {
#     signal <- sum(colNorms(A))
#     d <- nrow(A)
#     EN <- sqrt(2) * gamma((d + 1) / 2) / gamma(d / 2) # see exp. value of chi random variable
#     noise <- sigma * EN
#     ratio <- signal / (signal + noise)
#     return(ratio)
#   }
# }
# 
# # EVT -------------------------------------------------------------------------------
# 
# ## Univariate ------------------------------------------------------------------------
# 
# fit.gpd <- function(x, qthresh) {
#   u <- as.numeric(quantile(x, probs = qthresh))
#   fit <- evmix::fgpd(x = x, u = u)
#   data.mle <- tibble(u = u, sigma = fit$mle[1], xi = fit$mle[2])
#   return(data.mle)
# }
# 
# frechet2.inv.cdf <- function(x) (-log(x))^(-1/2)
# 
# ## Create data set and preprocessing ---------------------------------------------------------------------
# 
# preproc <- function(data, univar.qthresh = 0.95, inv.cdf = frechet2.inv.cdf) {
#   data %<>%
#     arrange(t) %>%
#     group_by(i) %>%
#     mutate(fit.gpd(Xit.orig, qthresh = univar.qthresh)) %>% # creates three columns: u, sigma, xi (see fit.gpd output)
#     mutate(Fihat = case_when(
#       Xit.orig <= u ~ rank(Xit.orig) / (n() + 1),
#       Xit.orig > u ~ univar.qthresh + (1 - univar.qthresh) * evmix::pgpd(Xit.orig - u, u = 0, sigmau = sigma, xi = xi))) %>%
#     mutate(Xit = inv.cdf(Fihat)) %>%
#     ungroup() %>%
#     group_by(t) %>%
#     mutate(rt = sqrt(sum(Xit^2))) %>%
#     ungroup() %>%
#     mutate(rt.quantile = rank(rt) / (n() + 1), Wit = Xit / rt) 
#   return(data)
# }
# 
# sim.data.A <- function(A, n, sigma = 0, univar.qthresh = 0.95) {
#   i.names <- rownames(A)
#   t.vals <- 1:n
#   X.vals <- rtlfm(n = n, A = A, sigma = sigma, alpha = 2) # n x d matrix
#   data <- data.frame(i = rep(i.names, times = n),
#                      t = rep(t.vals, each = length(i.names)),
#                      Xit.orig = matrix(t(X.vals), ncol = 1))
#   data <- preproc(data, univar.qthresh = univar.qthresh)
#   return(data)
# }
# 
# sim.data.br <- function(loc, vario, n, sigma = 0, univar.qthresh = 0.95) {
#   i.names <- rownames(loc)
#   t.vals <- 1:n
#   X.vals <- sim.br(n, loc, vario, t.names = 1:n, sigma = sigma) # n x d matrix
#   data <- data.frame(i = rep(i.names, times = n),
#                      t = rep(t.vals, each = length(i.names)),
#                      Xit.orig = matrix(t(X.vals), ncol = 1))
#   data <- preproc(data, univar.qthresh = univar.qthresh)
#   return(data)
# }
# 
# 
# # TPDM estimation -------------------------------------------------------------------
# 
# Sigma.mat.to.dat <- function(Sigma, value.name = "sigma.ij") {
#   Sigma[lower.tri(Sigma)] <- NA
#   Sigma.dat <- Sigma %>%
#     reshape2::melt(na.rm = TRUE, value.name = value.name) %>%
#     rename(i = Var1, j = Var2) %>%
#     arrange(i)
#   return(Sigma.dat)
# }
# 
# Sigma.dat.to.mat <- function(Sigma, value.var = "sigma.ij") {
#   Sigma.mat <- Sigma %>%
#     dcast(i ~ j, value.var = value.var) %>% 
#     column_to_rownames(var = "i") %>% 
#     as.matrix() %>%
#     forceSymmetric() %>%
#     as.matrix()
#   return(Sigma.mat)
# }
# 
# tpdm.est.std <- function(data, rt.qthresh = 0.95) {
#   Wext <- data %>% 
#     filter(rt.quantile > rt.qthresh) %>% # extreme obs.
#     dplyr::select(t, i, Wit) %>%
#     pivot_wider(names_from = i, values_from = Wit) %>%
#     dplyr::select(-t) %>%
#     as.matrix()
#   Sigma <- ncol(Wext) / nrow(Wext) * t(Wext) %*% Wext # (d/k) \sum_{i:extreme} w_ii w_i^T
#   Sigma <- Sigma.mat.to.dat(Sigma, value.name = "sigma.hat.ij")
#   Sigma <- as.tibble(Sigma)
#   out <- tibble("lambda" = NA, 
#                 "Sigma" = list(Sigma),
#                 "penalty.info" = NULL)
#   return(out)
# }
# 
# # test each pairwise combination of locations for extremal dependence
# # returns tibble [i, j, p.value]
# # p.value < alpha => dependent
# # p.value > alpha => independent
# test.ext.indep <- function(data, i1, i2) {
#   if (i1 == i2) {
#     p.value <- 0
#   } else {
#     test.tail <- data %>%
#       filter(i %in% c(i1, i2)) %>%
#       dplyr::select(i, t, Xit) %>%
#       pivot_wider(values_from = Xit, names_from = i) %>%
#       dplyr::select(-t) %>%
#       as.matrix() %>%
#       POT::tailind.test()
#     p.value <- test.tail$stats["ChiSq", "p.values"]
#   }
#   return(p.value)
# }
# 
# tpdm.est.lasso <- function(data, rt.qthresh = 0.95, lambda = NULL, pen.type = NULL, alpha.test = 0.05) {
#   i.set <- i.set(data)
#   d <- length(i.set)
#   n <- t.set(data) %>% length()
#   i.pairs <- make.i.pairs(i.set)
#   n.pairs <- nrow(i.pairs)
#   # penalty information
#   pen.dat <- i.pairs %>%
#     mutate(pen.type = pen.type) %>%
#     rowwise() %>%
#     mutate(p.value = test.ext.indep(data = data, i1 = i, i2 = j)) %>%
#     mutate(pen.factor = case_when(pen.type == "test" ~ as.numeric(p.value > alpha.test),
#                                   pen.type == "p.value" ~ p.value)) %>%
#     ungroup()
#   if (sum(pen.dat$pen.factor) < 1e-2) {
#     out <- tpdm.est.std(data = data, rt.qthresh = rt.qthresh)
#     print("All penalty factors are zero - standard TPDM estimator applied.")
#     return(out)
#   } else {
#     # response variable for lasso regression
#     Y <- matrix(NA, nrow = n, ncol = n.pairs)
#     data.ext <- data %>% filter(rt.quantile > rt.qthresh)
#     for (l in 1:n.pairs) {
#       i1 <- i.pairs$i[l]
#       i2 <- i.pairs$j[l]
#       Y[, l] <- data.ext %>%
#         filter(i %in% i.pairs[l, c("i", "j")]) %>%
#         group_by(t) %>%
#         summarise(y = case_when(i1 == i2 ~ Wit^2, i1 != i2 ~ prod(Wit))) %>% # Wit * Wjt (only one row in the filter if i1==i2)
#         ungroup() %>%
#         pull(y)
#     }
#     # estimate the model coefficients
#     fit <- glmnet(x = Matrix::bdiag(rep(list(rep(1, n)), n.pairs)),
#                   y = pracma::Reshape(Y, n = length(Y), m = 1), # stack columns of Y
#                   lambda = lambda,
#                   family = "gaussian",
#                   alpha = 1,
#                   intercept = FALSE,
#                   standardize = FALSE,
#                   penalty.factor = pen.dat$pen.factor)
# 
#     Sigma.fit <- cbind(i.pairs, as.matrix(fit$beta) * d) %>%
#       rename_with(function(sname) fit$lambda[as.numeric(gsub("s", "", sname)) + 1], starts_with("s")) %>%
#       pivot_longer(cols = -c(i, j), names_to = "lambda", values_to = "sigma.hat.ij") %>%
#       group_by(lambda)
#     out <- tibble("lambda" = pull(group_keys(Sigma.fit)),
#                   "Sigma" = as.list(group_split(Sigma.fit, keep = FALSE)),
#                   "penalty.info" = list(pen.dat))
#     return(out)
#   }
# }
# 
# tpdm.est <- function(data, rt.qthresh = 0.95, method = "standard", lambda = NULL, pen.type = NULL, alpha.test = 0.05) {
#   if (method == "standard") {
#     out <- tpdm.est.std(data = data, rt.qthresh = rt.qthresh)
#   } else {
#     out <- tpdm.est.lasso(data = data, rt.qthresh = rt.qthresh, lambda = lambda, pen.type = pen.type, alpha.test = alpha.test)
#   }
#   return(out)
# }
# 
# tpdm.error <- function(tpdm.est, tpdm.true){
#   tpdm.mse <- full_join(tpdm.est, tpdm.true, by = c("i", "j")) %>%
#     mutate(sigma.ij.err = sigma.ij - sigma.hat.ij) %>%
#     summarise(mse = mean(sigma.ij.err^2)) %>%
#     pull()
#   return(tpdm.mse)
# }
# 
# tpdm.eigen <- function(Sigma, fix.evec.sign = TRUE) {
#   if (!isSymmetric(Sigma)) stop("Sigma must be symmetric.")
#   Sigma <- as.matrix(Matrix::nearPD(Sigma, keepDiag = TRUE)$mat)
#   Sigma.eigen <- eigen(Sigma, symmetric = TRUE)
#   rownames(Sigma.eigen$vectors) <- colnames(Sigma)
#   if (fix.evec.sign) { # make first component of each eigenvector positive
#     Sigma.eigen$vectors <- t(t(Sigma.eigen$vectors) * ifelse(Sigma.eigen$vectors[1, ] >= 0, 1, -1))
#   }
#   return(Sigma.eigen)
# }
# 
# # Cross validation ------------------------------------------------------------------
# 
# create.ext.folds <- function(data, rt.qthresh, K) {
#   i.set <- i.set(data)
#   d <- length(i.set)
#   t.set <- t.set(data)
#   n <- length(t.set)
#   ext.t <- unique(data$t[data$rt.quantile > rt.qthresh])
#   folds <- rep(NA, length = n)
#   folds[ext.t] <- createFolds(y = ext.t, k = K, list = FALSE)
#   folds <- rep(folds, each = d)
#   data %<>% mutate(test.fold = folds)
#   return(data)
# }
# 
# # Simulating extremes ---------------------------------------------------------------
# 
# # rewrite this
# ext.pca <- function(data, lambda = 0, rt.qthresh = 0.95, reconstruct = TRUE) {
#   # TPDM and its eigendecomposition
#   #Sigma <- tpdm.lasso.estimate(data, lambda, rt.qthresh)
#   Sigma <- tpdm.estimate(data, rt.qthresh)
#   Sigma.eigen <- tpdm.eigen(Sigma, fix.evec.sign = TRUE)
#   # compute ePCs
#   Xmat <- data %>% 
#     dplyr::select(t, i, Xit) %>% # extract angular components
#     pivot_wider(names_from = i, values_from = Xit) # put into matrix-like format
#   tvals <- Xmat %>% pull(t) # save the t values
#   Xmat %<>%
#     dplyr::select(-t) %>% # remove t values from X to just leave the data
#     as.matrix()
#   V <- tau(Xmat, inverse = TRUE) %*% Sigma.eigen$vectors # V_t = U^T tau^{-1}(X_t)
#   dimnames(V) <- list(NULL, paste0("ePC.", 1:ncol(V)))
#   # compute low-dim reconstructions of X
#   if (reconstruct) {
#     X.rec <- vector(mode = "list", length = ncol(V))
#     for (m in 1:ncol(V)) {
#       X.rec[[m]] <- t(tau(Sigma.eigen$vectors[, 1:m] %*% t(V[, 1:m])))
#       colnames(X.rec[[m]]) <- colnames(Xmat)
#       X.rec[[m]] %<>%
#         as_tibble() %>%
#         mutate(t = tvals) %>%
#         pivot_longer(cols = !t, names_to = "i", values_to = paste("Xit.rec", m, sep = ".")) %>%
#         relocate(i, t)
#     }
#     X.rec <- Reduce(function(x, y) {full_join(x, y, by = c("i", "t"))}, X.rec) %>%
#       mutate(Xit = data$Xit) %>%
#       relocate(Xit, .after = t)
#   } else { # is !reconstruct 
#     X.rec <- NULL
#   }
#   # tidy up V
#   V <- as_tibble(V) %>%
#     mutate(t = tvals) %>%
#     relocate(t) %>%
#     pivot_longer(cols = starts_with("ePC"), names_to = "ePC.i", values_to = "Vit") %>%
#     group_by(t) %>%
#     mutate(rVt = sqrt(sum(Vit^2))) %>%
#     ungroup() %>%
#     mutate(rVt.quantile = rank(rVt) / (n() + 1),
#            WVit = Vit / rVt)
#   return(list("Sigma" = Sigma, 
#               "U" = Sigma.eigen$vectors,
#               "D" = Sigma.eigen$values,
#               "V" = V, 
#               "X.rec" = X.rec))
# }
# 
# sim.events.rc <- function(epca, nsim = nrow(epca$V), mvals = 1:ncol(epca$Sigma), rt.qthresh.V = 0.95) {
#   d <- ncol(epca$Sigma)
#   WV <- epca$V %>%
#     filter(rVt.quantile > rt.qthresh.V) %>%
#     dplyr::select(t, ePC.i, WVit) %>% # extract angular components
#     pivot_wider(names_from = ePC.i, values_from = WVit) %>% # put into matrix-like format
#     dplyr::select(-t) %>% # remove t values to just leave the data
#     as.matrix()
#   Z.mvals <- vector(mode = "list", length = length(mvals))
#   for (l in seq_along(mvals)) {
#     m <- mvals[l]
#     if (m == d) {
#       Z.mvals[[l]] <- WV
#     } else {
#       Z <- matrix(NA, nrow = nrow(WV), ncol = m + 1)
#       Z[, 1:m] <- WV[, 1:m]
#       Z[, m + 1] <- ifelse(WV[, m + 1] >= 0, sqrt(1 - colSums(t(Z[, 1:m]^2))), -sqrt(1 - colSums(t(Z[, 1:m]^2))))
#       Z.mvals[[l]] <- Z
#     }
#   }
#   kappa.vals <- lapply(Z.mvals, FUN = function(Z) Directional::vmfkde.tune(Z)[[1]]) %>% unlist()
#   # create sample for each m value
#   X.sim <- vector(mode = "list", length = length(mvals))
#   rV.sim <- evd::rfrechet(n = nsim, loc = 0, scale = sqrt(d), shape = 2) # use same radii for each m
#   for (l in seq_along(mvals)) {
#     m <- mvals[l]
#     kappa <- kappa.vals[l]
#     prob <- rep(1 / nrow(Z), nrow(Z))
#     k <- rep(1 / kappa^2, nrow(Z))
#     Z.sim <- Directional::rmixvmf(n = nsim, prob, mu = Z, k)$x
#     q <- max.col(Z.sim %*% t(Z)) # indexes q=(q1,...,qn) of nearest n'bours among Z to each Z.sim[i, ]
#     # simulate angles
#     WV.sim <- matrix(NA, nrow = nsim, ncol = d)
#     WV.sim[, 1:m] <- Z.sim[, 1:m]
#     WV.sim[, -(1:m)] <- abs(Z.sim[, m + 1] / Z[q, m + 1]) * WV[q, -(1:m)]
#     V.sim <- rV.sim * WV.sim
#     X.sim[[l]] <- t(tau(epca$U %*% t(V.sim)))
#   }
#   X.sim <- lapply(seq_along(X.sim), function(l) {
#     X <- X.sim[[l]]
#     X %<>% 
#       as_tibble() %>% 
#       mutate(t = 1:nsim) %>% 
#       pivot_longer(cols = !t, names_to = "i", values_to = "Xit.sim") %>%
#       mutate(m = mvals[l])
#   })
#   X.sim %<>% 
#     purrr::reduce(full_join, by = c("i", "t", "m", "Xit.sim")) %>%
#     group_by(t, m) %>%
#     mutate(rt.sim = sqrt(sum(Xit.sim^2))) %>%
#     ungroup() %>%
#     group_by(m) %>%
#     mutate(rt.sim.quantile = rank(rt.sim) / (n() + 1),
#            Wit.sim = Xit.sim / rt.sim) %>%
#     ungroup() %>%
#     relocate(i, t, m)
#   return(X.sim)
# }
# 
# # Packages --------------------------------------------------------------------------
# 
# # general R
# library(dplyr)
# library(tidyr)
# library(magrittr)
# library(tibble)
# library(purrr)
# library(reshape2)
# library(forcats)
# # EVT
# library(evmix)
# library(evd)
# library(POT)
# # linear algebra
# library(pracma)
# library(lattice)
# library(Matrix)
# library(abind)
# # general statistics
# library(stats)
# library(Directional)
# library(caret)
# # regularised estimation
# library(glmnet)
# library(genlasso)
# # plotting
# library(colorspace)
# library(ggplot2)
# library(GGally)
# library(patchwork)
# library(RColorBrewer)
# library(scales)
# 
# # TPDM estimation with TLFM data ----------------------------------------------------------------------
# 
# # moved to "Code/screpca-script.R"
# 
# # Plots -----------------------------------------------------------------------------
# 
# # # plot TPDM errors for each method and lambda value
# # est.results %>%
# #   dplyr::filter(n == 1000, sigma.noise == 0) %>%
# #   ggplot(aes(x = interaction(method, pen.type), y = Sigma.error, fill = forcats::fct_reorder(lambda, as.numeric(lambda)))) +
# #   geom_boxplot() +
# #   scale_fill_discrete_diverging("Blue-Red", labels = scientific(c(lambda.vals, NA), digits = 3)) + 
# #   scale_x_discrete(labels = c("Lasso: p-value penalty", "Lasso: accept/reject penalty", "Standard")) + 
# #   labs(x = "Estimator",
# #        y = "Error (Frobenius norm)",
# #        fill = expression("Penalty term strength," ~lambda)) + 
# #   theme_bw()
# # est.results %>%
# #   dplyr::filter(n == 1000, sigma.noise == 0.5) %>%
# #   ggplot(aes(x = interaction(method, pen.type), y = Sigma.error, fill = forcats::fct_reorder(lambda, as.numeric(lambda)))) +
# #   geom_boxplot() +
# #   scale_fill_discrete_diverging("Blue-Red", labels = scientific(c(lambda.vals, NA), digits = 3)) + 
# #   scale_x_discrete(labels = c("Lasso: p-value penalty", "Lasso: accept/reject penalty", "Standard")) + 
# #   labs(x = "Estimator",
# #        y = "Error (Frobenius norm)",
# #        fill = expression("Penalty term strength," ~lambda)) + 
# #   theme_bw()
# # 
# # # plot the best TPDM estimates for each method with rep=1, n=1000, sigma=0
# # optimal.df1 <- est.results %>%
# #   dplyr::filter(rep == 1, n == 1000, sigma.noise == 0) %>%
# #   group_by(method, pen.type) %>%
# #   slice(which.min(Sigma.error))
# # Sigma.hat.optimal <- optimal.df1$Sigma.hat %>%
# #   lapply(FUN = Sigma.dat.to.mat, value.var = "sigma.hat.ij") %>%
# #   simplify2array()
# # dimnames(Sigma.hat.optimal)[[3]] <- paste(optimal.df1$pen.type, optimal.df1$lambda)
# # 
# # # plot the errors in the eigenvalue estimates
# # eval.df <- est.results %>%
# #   dplyr::select(n, sigma.noise, rep, method, pen.type, lambda, Sigma.hat.evals, Sigma.true.evals) %>%
# #   unnest(Sigma.hat.evals, Sigma.true.evals) %>%
# #   group_by(n, sigma.noise, rep, method, pen.type, lambda) %>% 
# #   mutate(eval.index = row_number()) %>%
# #   ungroup()
# # lambda.star <- optimal.df1$lambda[1]
# # eval.df %>%
# #   dplyr::filter(n == 1000, sigma.noise == 0, lambda %in% c(NA, lambda.star)) %>%
# #   ggplot() + 
# #   geom_boxplot(aes(x = as.factor(eval.index), y = Sigma.true.evals - Sigma.hat.evals, fill = interaction(method, pen.type))) + 
# #   geom_hline(yintercept = 0) + 
# #   scale_fill_manual(values = c(alpha("red", 0.5), alpha("blue", 0.5), alpha("grey", 0.5)), labels = c("Lasso: p-value penalty", "Lasso: accept/reject penalty", "Standard")) + 
# #   labs(x = "Index, i",
# #        y = expression("Error in eigenvalue estimate," ~hat(~lambda)[i] - ~lambda[i]),
# #        fill = "TPDM estimator") + 
# #   theme_bw()
# 
