
vect_of_thresholds <- function(v){
  # This function gives a vector listing the ordered thresholds for which
  # a vector v has a different number of nonnull coordinates
  d <- length(v)
  u <- sort(v, TRUE)
  threshold <- rep(0,d)
  su <- cumsum(u)
  
  for(k in 1:(d-1)){
    threshold[k] <- su[k]-k*u[k+1]}
  threshold[d] <- su[d] # threshold above which the projection does not make sense
  
  return(threshold)
}

cones_function <- function(v, Gamma){
  # This function gives the subsets which contains the mass of Z
  # with Gamma by decreasing order
  d <- length(v)
  ord <- order(v)
  thres <- vect_of_thresholds(v)
  length_Gamma <- length(Gamma)
  
  J <- which(Gamma < thres[d], arr.ind = TRUE) # the coordinates of the useful Gamma
  length_J <- length(J)
  # the one that are too big provide a vector of NA
  
  cones <- matrix(NA, nrow = d, ncol = length_Gamma)
  if (length_J > 0){
    # Initialization
    j <- J[1]
    r1 <- sum(thres<Gamma[j]) + 1 # number of positive coordinates for the threshold Gamma[1]
    cones[ord[1:(d-r1)], j] <- 0 # the coordinates on which v is not projected
    cones[ord[seq2(d-r1+1, d)], j] <- 1 # the coordinates on which v is projected
    
    j <- J[2]
    r2 <- d+2 # just need to be above d+1 I think
  }
  
  if (length_J > 1){
    # the loop
    while ((r2 >= 1)&&(j<=length_Gamma)) {
      r2 <- sum(thres<Gamma[j]) + 1 # number of positive coordinates for the threshold Gamma[1]
      if (r2 <= d){
        cones[ , j] <- cones[ , j-1]
        cones[ord[seq2(d-r1+1, d-r2)], j] <- 0 # seq2 only works for increasing sequence which helps if r1 < r2
      } # we let the NA if the Gamma[1] is too huge
      j <- j + 1
      r1 <- r2
    }
  }
  return(cones)
}

occ_subsets <- function(testmtx){
  # This function takes a matrix and gives the number of occurence of each column
  nc <- ncol(testmtx)
  occ <- seq(nc)  
  for (i in seq(nc-1)) {
    dup <- colSums( abs( testmtx[,seq(i+1,nc),drop=F] - testmtx[,i] ) ) == 0
    occ[which(dup)+i] <- occ[i]
  }
  result <- rbind(as.matrix(testmtx[ ,unique(occ)]), table(occ))
  colnames(result) <- NULL
  return(result)
}

muscle_plot <- function(X, prop){
  # This function gives for each k the value of s_hat which minimizes the KL
  # and the value of this minimizer
  n <- ncol(X)
  d <- nrow(X)
  p <- length(prop) # prop should be given in an increasing order
  
  sum_norm <- apply(X, 2, sum)
  sort_sum_norm <- sort(sum_norm, decreasing = TRUE)
  thresholds <- sort_sum_norm[round(n*prop+1, 0)] # this corresponds to the list of thresholds u_n
  # we use the function 'round' since otherwise there are some approximations in n*Proportion
  X <- X[ , sum_norm > thresholds[p]] # this is done only to avoid computing too many NA's
  
  length_thresh <- length(thresholds)
  
  # Initialization: the "adjacency matrix"
  binary <- apply(X, 2, cones_function, thresholds)
  s_hat <- rep(NA, length_thresh)
  result <- matrix(NA, ncol = 4)
  colnames(result) <- c("k", "s_tilde", "hat_s", "minimizer")
  
  for (j in 1:length_thresh){
    binary_bis <- binary[( (j-1)*d + 1 ): (j*d) , !is.na(binary[(j-1)*d+1, ]) ] # we remove the columns of NA of the considered block
    k <- ncol(binary_bis) # it corresponds to k_n
    M <- occ_subsets(binary_bis)
    r <- ncol(M)
    
    # We order the matrix M in increasing order regarding the occurence of the columns
    ordered_subcones <- order(M[d+1, ], decreasing=TRUE)
    M <- as.matrix(M[ , ordered_subcones])
    T <- M[(d+1), ]
    
    # we now optimize in s_hat
    optim <- 1:(r-1) - lfactorial(k) + k*log(k) + sum(lfactorial(T)) - cumsum(T[-r]*(log(T)[-r])) - (k-cumsum(T[-r]))*log( (k-cumsum(T[-r])) / ((r-1):1) )
    s_hat <- which.min(optim)
    
    minimizer <- (optim[s_hat])/k + k/n
    result <- rbind(result, c(k, r, s_hat, minimizer))
  }
  
  result <- result[-1, , drop=F] # we have to remove the first row of NA's
  return(result)
}

##########################################################

# note: X is dxn, so use with X = t(X)
muscle_clusters <- function(X, prop){
  # This function uses the previous algo, chooses the optimal k,
  # and gives the associated subsets on which the angular vector Z puts mass
  n <- ncol(X)
  d <- nrow(X)
  p <- length(prop) # prop should be given in an increasing order
  
  sum_norm <- apply(X, 2, sum)
  sort_sum_norm <- sort(sum_norm, decreasing = TRUE)
  
  thresholds <- sort_sum_norm[round(n*prop+1, 0)] # this corresponds to the list of thresholds u_n
  # we use the function 'round' since otherwise there are some approximations in n*Proportion
  X <- X[ , sum_norm > thresholds[p]] # this is done only to avoid computing too many NA's
  
  length_thresh <- length(thresholds)
  
  # Initialization: the "adjacency matrix"
  binary <- apply(X, 2, cones_function, thresholds)
  minimizer <- Inf
  s_hat <- NA
  k_hat <- NA
  
  for ( j in 1:length_thresh ){
    binary_bis <- binary[( (j-1)*d + 1 ): (j*d) , !is.na(binary[(j-1)*d+1, ]) ] # we remove the columns of NA of the considered block
    k <- ncol(binary_bis) # it corresponds to k_n
    M <- occ_subsets(binary_bis)
    r <- ncol(M)
    
    # We order the matrix M in increasing order regarding the occurence of the columns
    ordered_subcones <- order(M[d+1, ], decreasing=TRUE)
    M <- as.matrix(M[ , ordered_subcones])
    T_val <- M[(d+1), ]
    
    # we now optimize in s_hat
    optim <- 1:(r-1) - lfactorial(k) + k*log(k) + sum(lfactorial(T_val)) - cumsum(T_val[-r]*(log(T_val)[-r])) - (k-cumsum(T_val[-r]))*log( (k-cumsum(T_val[-r])) / ((r-1):1) )
    min_loc <- optim[which.min(optim)]/k + k/n
    
    if (minimizer > min_loc){
      minimizer <- min_loc
      s_hat <- which.min(optim)
      k_hat <- k
      Mat <- M[, 1:s_hat]
    }
  }
  M <- as.matrix(Mat)
  u <- sort_sum_norm[round(k_hat + 1, 0)]
  weights <- M[(d+1), ]/ sum(M[(d+1), ])
  return(list("M" = M, "k_hat" = k_hat, "u" = u, "s_hat" = s_hat, "weights" = weights))
}

tpdm_estimate_bias <- function(X, k) {
  
  euc_proj <- function(x) {
    b <- 1 # radius of the simplex
    u <- sort(x, decreasing = TRUE)
    sx <- cumsum(u)
    rho <- which(u > (sx - b) / (1:length(x)))
    theta <- max(0, (sx[rho] - b) / rho)
    w <- x - theta
    w[w < 0] = 0
    return(w)
  }
  
  extreme_angles <- function(X, k, alpha = 1, type = "selfnorm", order = "time") {
    
    # compute and sort radii
    r <- apply(X, 1, pracma::Norm, p = alpha)
    sort_r <- sort(r, decreasing = TRUE)
    
    # radial threshold
    if (k > nrow(X)) stop("k must be <= n.")
    rk <- sort_r[k]
    ext_ind <- r >= rk
    X_ext <- X[ext_ind, , drop = FALSE]
    
    # compute angles
    if (type == "selfnorm") {
      W <- X_ext / r[ext_ind]
    } else if (type == "euclidean") {
      if (alpha != 1) stop("You are using Euclidean projection but alpha != 1.")
      if (k > nrow(X) - 1) stop("k must be <= n-1 for Euclidean projection.")
      rk1 <- sort_r[k + 1]
      W <- apply(X_ext / rk1, 1, euc_proj) %>% t()
    } else {
      stop('type must be "selfnorm" or "euclidean".')
    }
    dimnames(W) <- dimnames(X_ext)
    if (order == "norm") W <- W[names(sort_r[1:k]), , drop = FALSE]
    return(W)
  }
  
  extreme_directions <- function(X, k) {
    
    W <- extreme_angles(X, k, alpha = 1, type = "euclidean")
    W_beta <- apply(W > 0, 1, function(w) paste(w * 1, collapse = ""))
    df <- W_beta %>%
      table() %>% 
      sort(decreasing = TRUE) %>%
      tibble::enframe(name = "beta", value = "freq") %>%
      dplyr::mutate("beta_size" = stringr::str_count(beta, "1")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate("t" = list(which(W_beta == beta))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate("s" = seq_len(nrow(.)))
    
    # log-likelihood (Eq 4.3) and model selection criterion (Eq 4.6) in Meyer and Wintenberger (2023)
    r <- nrow(df)
    l_constant <- lfactorial(k) - sum(lfactorial(df$freq))
    df$loglkhood_s <- sapply(df$s, function(s) {
      l_constant + cumsum(df$freq * log(df$freq / k))[s] + sum(df$freq[(s+1):r]) * log((1 - cumsum(df$freq / k)[s]) / (r - s))
    })
    df$meyer_criterion <- (k / nrow(X)) - ((df$s - df$loglkhood_s) / k)
    
    df$Sigma <- lapply(seq_along(df$s), function(s) {
      inds <- df$t[df$s <= s] %>% unlist()
      Sigma_hat <- (ncol(X) / length(inds)) * t(W[inds, ])^(1/2) %*% W[inds, ]^(1/2)
      return(Sigma_hat)
    })
    
    return(df)
  }
  
  d <- ncol(X)
  muscle_out <- muscle_clusters(t(X), prop = k / nrow(X)) 
  W <- extreme_angles(X, k = muscle_out$k_hat, alpha = 1, type = "euclidean")
  
  M_beta <- lapply(seq_len(ncol(muscle_out$M)), function(i) muscle_out$M[1:d, i])
  W_beta <- 1 * (W > 0)
  
  # APPROACH 1: only keep angles that are contained \hat{\mathcal{S}^\star}
  
  # sparse estimate
  keep_angles <- which(apply(W_beta, 1, function(wb) any(M_beta %in% list(as.numeric(wb)))) > 0) %>% as.numeric()
  A_keep <- d / length(keep_angles) * t((W[keep_angles, , drop = FALSE]))
  A_keep_alpha <- A_keep^(1/2)
  S <- A_keep_alpha %*% t(A_keep_alpha)
  
  dimnames(S) <- list(colnames(X), colnames(X))
  print(paste0("k_hat = ", muscle_out$k_hat, ", s_hat = ",  muscle_out$s_hat, ", No. keep angles = ", length(keep_angles)))
  return(S)
}
