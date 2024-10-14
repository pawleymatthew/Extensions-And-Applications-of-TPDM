library(pracma, quietly = TRUE)
library(evmix, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(purrr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(reshape2, quietly = TRUE)
library(tibble, quietly = TRUE)
library(Ternary, quietly = TRUE)
library(scales, quietly = TRUE)
library(compositions, quietly = TRUE)
library(robCompositions, quietly = TRUE)
library(lsa, quietly = TRUE)

# General ---------------------------------------------------------------------------

rowNorms <- function(M) {
  if (is.vector(M)) M <- matrix(M, nrow = 1)
  r <- sqrt(rowSums(M^2))
  return(r)
}

tau <- function(x, inverse = FALSE) {
  # softplus function and its inverse
  if (!inverse){
    y <- -plogis(x, lower.tail = FALSE, log.p = TRUE) # log(1 + e^x)
  } else {
    y <- qlogis(-x, lower.tail = FALSE, log.p = TRUE) # log(e^x - 1)
  }
  return(y)
}

# Brown-Resnick ---------------------------------------------------------------------

sim.br <- function(n, loc, vario, t.names = 1:n) {
  X <- mvPot::simulBrownResnick(n = n, loc = loc, vario = vario) %>%
    unlist() %>% 
    matrix(nrow = n, ncol = nrow(loc), byrow = TRUE, dimnames = list(t.names, rownames(loc)))
  return(X)
}

# see Cooley supplementary material
tpdm.br <- function(loc, vario) {
  # compute distance between the locations
  D <- dist(loc) %>% as.matrix()
  # function: map distance to extremal dependence measure (TPDM) under BR model
  # helper function
  h <- function(x, y, s) {
    a <- (2 * vario(s))^(1 / 2)
    f1 <- (a / 2) - (1 / a) * log(x^2 / y^2)
    f2 <- (a / 2) - (1 / a) * log(y^2 / x^2)
    hnum <- 2 * (y^2 * (1 - f1 / a) * dnorm(f1) + x^2 * (1 - f2 / a) * dnorm(f2))
    hden <- a * x^3 * y^3
    h <- hnum / hden
    return(h)
  }
  # function for computing elements of Sigma
  tpdm.fun <- function(d) {
    if (d == 0) {
      sigma <- 1
    } else {
      htilde <- function(theta) h(cos(theta), sin(theta), d)
      sigma <- integrate(f = function(theta) cos(theta) * sin(theta) * htilde(theta), lower = 0, upper = pi/2)$value
    }
    return(sigma)
  }
  # compute Sigma
  Sigma <- D %>% 
    sapply(tpdm.fun) %>%
    matrix(nrow = nrow(D), ncol = nrow(D), dimnames = list(rownames(loc), rownames(loc))) # sym. matrix with diag=0
  return(Sigma)
}

# Preprocessing ---------------------------------------------------------------------

fit.gpd <- function(x, qthresh = 0.95) {
  u <- quantile(x, probs = qthresh) %>% as.numeric()
  fit <- evmix::fgpd(x = x, u = u)
  data.mle <- tibble(u = u, sigma = fit$mle[1], xi = fit$mle[2])
  return(data.mle)
}

frechet2.inv.cdf <- function(x) (-log(x))^(-1/2)
frechet1.inv.cdf <- function(x) (-log(x))^(-1)

mvts.to.long <- function(X, value.name = deparse(substitute(X))) {
  # input: 
  # X = n x d matrix/dataframe of values
  # value.name = what the value column will be called (".it" is then appended to the end)
  # output: nd x 3 dataframe, [t; i; X.it]
  if (is.matrix(X)) X <- as.data.frame(X)
  n <- nrow(X)
  d <- ncol(X)
  X.long <- reshape2::melt(X, value.name = paste(value.name, "it", sep = "."), id.vars = NULL) %>%
    rename(i = variable) %>%
    mutate(t = rep(1:n, d)) %>%
    relocate(t, i)
  return(X.long)
}

std.margins <- function(X, univar.qthresh = 0.95, inv.cdf = frechet2.inv.cdf) {

  temp <- X %>% 
    mvts.to.long(value.name = "X") %>%
    dplyr::group_by(i) %>%
    dplyr::mutate(fit.gpd(X.it, qthresh = univar.qthresh)) %>%
    dplyr::mutate(Fhat.i = dplyr::case_when(
      X.it <= u ~ rank(X.it) / (n() + 1),
      X.it > u ~ univar.qthresh + (1 - univar.qthresh) * evmix::pgpd(X.it - u, u = 0, sigmau = sigma, xi = xi))) %>%
    dplyr::mutate(X.std.it = inv.cdf(Fhat.i)) %>%
    dplyr::ungroup()
  
  X.std <- temp %>%
    dplyr::select(i, t, X.std.it) %>%
    reshape2::acast(t ~ i, value.var = "X.std.it") %>%
    as.data.frame()
  
  out <- list("X.std" = X.std, "std.info" = temp)
  return(out)
}

euc.proj <- function(x) {
  b <- 1 # radius of the simplex
  u <- sort(x, decreasing = TRUE)
  sx <- cumsum(u)
  rho <- which(u > (sx - b) / (1:length(x)))
  theta <- max(0, (sx[rho] - b) / rho)
  w <- x - theta
  w[w < 0] = 0
  return(w)
}

ext.polar <- function(X, proj = "l1.selfnorm", R.qthresh = 0.95) {
  # only the Euclidean projection relies on defining a threshold t
  if (proj == "l1.selfnorm") {
    R <- rowSums(X)
    Theta <- X / R
    R.thresh <- as.numeric(quantile(R, R.qthresh))
  } else if (proj == "l2.selfnorm") {
    R <- rowNorms(X)
    Theta <- X / R
    R.thresh <- as.numeric(quantile(R, R.qthresh))
  } else if (proj == "pi1.euc") {
    R <- rowSums(X)
    R.thresh <- as.numeric(quantile(R, R.qthresh))
    Theta <- t(apply(X / R.thresh, 1, euc.proj))
  } else {
    return("Enter a valid projection method.")
  }
  ext.ind <- which(R > R.thresh)
  R.ext <- R[ext.ind]
  Theta.ext <- Theta[ext.ind, ]
  X.pp <- list("R" = R, 
               "Theta" = Theta,
               "ext.ind" = ext.ind,
               "R.ext" = R.ext, 
               "Theta.ext" = Theta.ext, 
               "proj" = proj, 
               "R.thresh" = R.thresh)
  return(X.pp)
}

tpdm.standard <- function(X.pp) {
  W <- as.matrix(X.pp$Theta.ext)
  Sigma <- ncol(W) / nrow(W) * t(W) %*% W
  return(Sigma)
}

# PCA -------------------------------------------------------------------------------

# input matrix X has standardised columns (alpha = 1)
ext.pca.comp <- function(X.pp) {
  # create compositional data matrix (class=acomp)
  Theta.ext <- X.pp$Theta.ext %>% acomp()
  # dimensions
  k <- nrow(Theta.ext) # number of extreme observations
  d <- ncol(Theta.ext) # number of dimensions
  # do compositional PCA
  xi <- geometricmeanCol(Theta.ext) %>% clo() %>% acomp() # data centre, class = acomp
  Theta.ext.ctd <- Theta.ext - xi # centring (subtraction in the simplex)
  clr.Theta.ext.ctd <- clr(Theta.ext.ctd) # map to CLR space
  Gamma <- cov(clr.Theta.ext.ctd) # centred log-ratio covariance matrix
  ssvd <- svd(clr.Theta.ext.ctd) # 'simplicial SVD' (in CLR space)
  ssvd$d <- c(as.numeric(ssvd$d), rep(0, d - length(ssvd$d)))
  # create outputs
  var.pc <- data.frame("PC" = 1:d,
                       "sing.value" = ssvd$d,
                       "eigenvalue" = ssvd$d^2 / (k - 1),
                       "prop.var" = ssvd$d^2 / sum(ssvd$d^2),
                       "cum.prop.var" = cumsum(ssvd$d^2) / sum(ssvd$d^2))
  Gamma <- cov(clr.Theta.ext.ctd) # centred log-ratio covariance matrix
  V <- matrix(ssvd$v, nrow = d, ncol = d, dimnames = list(colnames(Theta.ext), paste0("PC", 1:d))) # evectors (in CLR space)
  U <- matrix(ssvd$u, nrow = k, ncol = d) # left sing. vectors
  
  # outputs in list
  out <- list("xi" = xi,
              "Gamma" = Gamma,
              "var.pc" = var.pc,
              "V" = V, 
              "U" = U,
              "pca.method" = "compositional")
  return(out)
}

ext.pca.sabourin <- function(X.pp, subset = colnames(Theta.ext)) {
  # create data matrix of angles (class=matrix)
  Theta.ext <- as.matrix(X.pp$Theta.ext)
  Theta.ext <- Theta.ext[, subset] / rowNorms(Theta.ext[, subset])
  # dimensions
  k <- nrow(Theta.ext) # number of extreme observations
  d <- ncol(Theta.ext) # number of dimensions
  # do PCA
  Sigma <- d / k * t(Theta.ext) %*% Theta.ext
  eigen.Sigma <- eigen(Sigma)
  lambda <- eigen.Sigma$values
  # create outputs
  var.pc <- data.frame("PC" = 1:d,
                       "sing.value" = NA,
                       "eigenvalue" = lambda,
                       "prop.var" = lambda / sum(lambda),
                       "cum.prop.var" = cumsum(lambda) / sum(lambda))
  V <- matrix(eigen.Sigma$vectors, nrow = d, ncol = d, dimnames = list(colnames(Theta.ext), paste0("PC", 1:d))) # evectors
  U <- NULL
  # outputs in list
  out <- list("Sigma" = Sigma,
              "var.pc" = var.pc,
              "V" = V, 
              "U" = U,
              "pca.method" = "sabourin")
  return(out)
}

ext.pca.cooley <- function(X.pp) {
  # create directional data matrix (class=matrix)
  Theta.ext <- as.matrix(X.pp$Theta.ext)
  # dimensions
  k <- nrow(Theta.ext) # number of extreme observations
  d <- ncol(Theta.ext) # number of dimensions
  # do PCA
  Sigma <- d / k * t(Theta.ext) %*% Theta.ext
  eigen.Sigma <- eigen(Sigma)
  lambda <- eigen.Sigma$values
  # create outputs
  var.pc <- data.frame("PC" = 1:d,
                       "sing.value" = NA,
                       "eigenvalue" = lambda,
                       "prop.var" = lambda / sum(lambda),
                       "cum.prop.var" = cumsum(lambda) / sum(lambda))
  V <- matrix(eigen.Sigma$vectors, nrow = d, ncol = d, dimnames = list(colnames(Theta.ext), paste0("PC", 1:d))) # evectors
  U <- NULL
  ePC <- tau(as.matrix(X.pp$Theta) * as.numeric(X.pp$R), inverse = TRUE) %*% V # V_t = U^T tau^{-1}(X_t)
  dimnames(ePC) <- list(NULL, paste0("ePC", 1:ncol(V)))
  # outputs in list
  out <- list("Sigma" = Sigma,
              "var.pc" = var.pc,
              "V" = V,
              "U" = U,
              "ePC" = ePC,
              "pca.method" = "cooley")
  return(out)
}

# X.pp and ext.pca should be compositional type
event.rec.comp <- function(X.pp, ext.pca){
  # preallocate kxdxd array with the data centre xi in each slice
  k <- dim(X.pp$Theta.ext)[1]
  d <- dim(X.pp$Theta.ext)[2]
  X.ext.rec <- array(1, dim = c(k, d, d), dimnames = append(dimnames(X.pp$Theta.ext), list(paste0("m=", 0:(d-1)))))
  for (i in 1:k) {
    X.ext.rec[i, , 1] <- ext.pca$xi
  }
  for (m in 2:d) { # loop over reconstruction rank m=1,...,d
    for (i in 1:k) { # loop over event i = 1,...,k
      X.ext.rec[i, , m] <- ext.pca$xi
      for (j in 1:(m-1)) {
        svd.pertpow.ij <- (ext.pca$U[i, j] * ext.pca$var.pc$sing.value[j]) * clrInv(ext.pca$V[, j])
        X.ext.rec[i, , m] <- X.ext.rec[i, , m] + svd.pertpow.ij
        X.ext.rec[i, , m] <- as.numeric(X.pp$R.ext[i]) * as.numeric(X.ext.rec[i, , m])
      }
    }
  }
  return(X.ext.rec)
}

# X.pp and ext.pca should be Sabourin/L2 type
# outputs are reconstructions of the L2 angle as vectors in Rd
# for m<d, reconstruction could have negative values and rowNorms!=1
event.rec.sabourin <- function(X.pp, ext.pca){
  # preallocate kxdxd array with the data centre xi in each slice
  k <- dim(X.pp$Theta.ext)[1]
  d <- dim(X.pp$Theta.ext)[2]
  Theta.ext.rec <- array(1, dim = c(k, d, d), dimnames = append(dimnames(X.pp$Theta.ext), list(1:d)))
  for (m in 1:d) {
    proj.subspace <- ext.pca$V[, 1:m] # matrix whose columns span a subspace of some R^d
    proj.points <- X.pp$Theta.ext %>% as.matrix() %>% t() # matrix whose COLUMNS are to be projected
    proj.out <- pracma::linearproj(A = proj.subspace, B = proj.points)$Q %>% t() # matrix whose rows are the projected points in R^d
    #proj.out <- pmax(proj.out, min(proj.out[proj.out > 0])) # projection may have non-positive values; convert these to some small value
    #proj.out <- pmax(proj.out, 1e-16)
    Theta.ext.rec[, , m] <- proj.out
    #Theta.ext.rec[, , m] <- proj.out / rowNorms(proj.out) # matrix whose rows are the projected points on unit circle
  }
  return(Theta.ext.rec)
}

# X.pp and ext.pca should be Cooley/L2 type
event.rec.cooley <- function(X.pp, ext.pca){
  # preallocate kxdxd array with the data centre xi in each slice
  k <- dim(X.pp$Theta.ext)[1]
  d <- dim(X.pp$Theta.ext)[2]
  X.ext.rec <- array(1, dim = c(k, d, d), dimnames = append(dimnames(X.pp$Theta.ext), list(1:d)))
  for (m in 1:d) {
    X.ext.rec[, , m] <- t(tau(ext.pca$V[, 1:m] %*% t(ext.pca$ePC[X.pp$ext.ind, 1:m]))) # kxd recon of all extreme X's
  }
  return(X.ext.rec)
}

# error.type = "aitchison" (on simplex), "euclidean", "cosine" (on unit sphere)
rec.err <- function(X.ext.true, X.ext.rec, err.type = "euclidean") {
  k <- dim(X.ext.rec)[1]
  d <- dim(X.ext.rec)[3]
  rec.err <- matrix(NA, nrow = k, ncol = d) # ij entry = error in j-rank recon. of event i
  for (i in 1:k) {
    for (m in 1:d) {
        if (err.type == "euclidean") {
          rec.err[i, m] <- dist(rbind(X.ext.true[i, ], X.ext.rec[i, , m]))[1]
        } else if (err.type == "cosine") {
          rec.err[i, m] <- 1 - lsa::cosine(as.numeric(X.ext.rec[i, , m]), as.numeric(X.ext.true[i, ]))
        }
      }
    }
  return(rec.err)
}

# # legend("topright", legend = paste0("PC", 1:tern.pcs), col = pc.cols[1:tern.pcs], lwd = 2)