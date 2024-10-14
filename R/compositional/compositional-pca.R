comp_line <- function(x0, x1, tau = seq(from = -10, to = 10, by = 0.1)) {
  dim <- length(x0)
  n_pts <- length(tau)
  pts <- matrix(NA, nrow = n_pts, ncol = dim)
  sapply(tau, function(t) perturbe(x0, power.acomp(x1, t)))
}

extreme_pca <- function(X_train, X_test, k, method = "coda", n_pcs = 1, plot = FALSE, failure_u = 100, ...) {
  
  if (method == "coda") {
    n_pcs <- n_pcs[n_pcs < ncol(X_train)]
  }
  
  # do PCA on training data
  R_train <- apply(X_train, 1, pracma::Norm, p = ifelse(method == "coda", 1, 2))
  R_thresh <- Rfast::nth(R_train, k = k + 1, descending = TRUE)
  ext_train <- which(R_train > R_thresh)
  Theta_ext_train <- X_train[ext_train, ] / R_train[ext_train]
  if (method == "coda") Theta_ext_train <- acomp(Theta_ext_train)
  pc <- princomp(Theta_ext_train)
  
  # reconstruct test data above threshold
  R_test <- apply(X_test, 1, pracma::Norm, p = ifelse(method == "coda", 1, 2))
  ext_test <- which(R_test > R_thresh)
  Theta_ext_test <- X_test[ext_test, ] / R_test[ext_test]
  if (method == "coda") Theta_ext_test <- acomp(Theta_ext_test)
  
  Theta_ext_train_hat <- lapply(seq_along(n_pcs), function(i) {
    m <- n_pcs[i]
    pred_train <- predict(pc, newdata = Theta_ext_train)
    if (method == "coda") {
      Theta_hat <- pc$Center
      for (j in 1:m) {
        Theta_hat <- perturbe(Theta_hat, power.acomp(clrInv(pc$loadings[, j]), pred_train[, j]))
      }
    } else {
      Theta_hat <- pc$center + (pred_train[, 1:m, drop = FALSE] %*% t(pc$loadings[, 1:m, drop = FALSE]))
    }
    return(Theta_hat)
  })
  
  Theta_ext_test_hat <- lapply(seq_along(n_pcs), function(i) {
    m <- n_pcs[i]
    pred_test <- predict(pc, newdata = Theta_ext_test)
    if (method == "coda") {
      Theta_hat <- pc$Center
      for (j in 1:m) {
        Theta_hat <- perturbe(Theta_hat, power.acomp(clrInv(pc$loadings[, j]), pred_test[, j]))
      }
    } else {
      Theta_hat <- pc$center + (pred_test[, 1:m, drop = FALSE] %*% t(pc$loadings[, 1:m, drop = FALSE]))
    }
    return(Theta_hat)
  })
  
  if (plot) {
    colpal <- c("red", "blue", "orange")
    # training data
    TernaryPlot(grid.lines = 0, ...)
    TernaryPoints(Theta_ext_train, pch = 20, cex = 1, col = "black")
    for (i in seq_along(n_pcs)) {
      TernaryPoints(Theta_ext_train_hat[[i]], cex = 2, pch = 4, col = alpha(colpal[n_pcs[i]], 0.5))
    }
    if (method == "coda") {
      TernaryLines(comp_line(x0 = pc$Center, # clrInv(princomp(Theta)$center)
                             x1 = clrInv(pc$loadings[, 1])),
                   col = colpal[1], lty = 2, lwd = 4)
      TernaryLines(comp_line(x0 = pc$Center, # clrInv(princomp(Theta)$center)
                             x1 = clrInv(pc$loadings[, 2])),
                   col = colpal[2], lty = 2, lwd = 4)
    }
    # # test data
    # TernaryPlot(grid.lines = 0, ...)
    # TernaryPoints(Theta_ext_test, pch = 20, cex = 1, col = "black")
    # for (i in seq_along(n_pcs)) {
    #   TernaryPoints(Theta_ext_test_hat[[i]], pch = 1, col = colpal[n_pcs[i]])
    # }
    # if (method == "coda") {
    #   TernaryLines(comp_line(x0 = pc$Center, # clrInv(princomp(Theta)$center)
    #                          x1 = clrInv(pc$loadings[, 1])),
    #                col = colpal[1], lty = 2, lwd = 2)
    #   TernaryLines(comp_line(x0 = pc$Center, # clrInv(princomp(Theta)$center)
    #                          x1 = clrInv(pc$loadings[, 2])),
    #                col = colpal[2], lty = 2, lwd = 2)
    # }
  }
  
  # mean squared Aitchison reconstruction error
  aitchison_err <- sapply(seq_along(n_pcs), function(i) {
    norm(acomp(Theta_ext_test) - acomp(Theta_ext_test_hat[[i]]))^2 %>% mean(na.rm = TRUE) 
  })
  
  euclidean_err <- sapply(seq_along(n_pcs), function(i) {
    apply(
      data.frame(acomp(Theta_ext_test)) - data.frame(acomp(Theta_ext_test_hat[[i]])), 
      1, 
      function(x) norm(x)^2) %>% 
      mean(na.rm = TRUE) 
  })
  
  p_max_hat <- sapply(seq_along(n_pcs), function(i) {
    A_hat <- Theta_ext_train_hat[[i]] %>% acomp() %>% data.frame() %>% as.matrix()
    A_hat <- A_hat[apply(A_hat, 1, function(x) all(x > 0)), ]
    A_hat <- ncol(A_hat) / nrow(A_hat) * A_hat
    apply(A_hat, 1, function(a) max(a / failure_u)) %>% sum()
  })
  
  p_min_hat <- sapply(seq_along(n_pcs), function(i) {
    A_hat <- Theta_ext_train_hat[[i]] %>% acomp() %>% data.frame() %>% as.matrix()
    A_hat <- A_hat[apply(A_hat, 1, function(x) all(x > 0)), ]
    A_hat <- ncol(A_hat) / nrow(A_hat) * A_hat
    apply(A_hat, 1, function(a) min(a / failure_u)) %>% sum()
  })
  
  tmp <- list(
    "Theta_ext_train" = Theta_ext_train,
    "Theta_ext_test" = Theta_ext_test,
    "Theta_ext_train_hat" = Theta_ext_train_hat,
    "Theta_ext_test_hat" = Theta_ext_test_hat,
    "threshold" = R_thresh,
    "k_test" = nrow(Theta_ext_test),
    "lambda" = pc$sdev^2,
    "test_loss_aitchison" = aitchison_err,
    "test_loss_euclidean" = euclidean_err,
    "p_max_hat" = p_max_hat,
    "p_min_hat" = p_min_hat
  )
  
  return(tmp)
}


# Sim AxZ ---------------------------------------------------------------------------

tau <- function(x, inverse = FALSE) {
  # softplus function and its inverse
  if (!inverse){
    y <- -plogis(x, lower.tail = FALSE, log.p = TRUE) # log(1 + e^x)
  } else {
    y <- qlogis(-x, lower.tail = FALSE, log.p = TRUE) # log(e^x - 1)
  }
  return(y)
}

rmaxlinearfactor <- function(n, A, type = "maxlin", alpha = 1) {
  
  if (type == "maxlin") {
    X <- rmaxlin(n = n, dsgn.mat = A)
  } else {
    X <- matrix(NA, nrow = n, ncol = nrow(A))
    for (i in 1:n) {
      Z <- evd::rgpd(n = ncol(A), loc = 0, scale = 1, shape = 1 / alpha)
      X[i, ] <- tau(A %*% tau(Z, inverse = TRUE))
    }
  }
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  return(X)
}
