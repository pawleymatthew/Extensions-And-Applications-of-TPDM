rm(list = ls())

sapply(list.files(path = "R/general", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
sapply(list.files(path = "R/changing-ext-dep", pattern = "*.R", full.names = TRUE, recursive = TRUE), source)
set.seed(1)


# Logistic type models --------------------------------------------------------------

set.seed(1)

n <- 5 * 10^4
d <- 3
k <- sqrt(n)

# symmetric logistic
par <- 0.3

# asymmetric logistic
asy <- list(0, .0, 0, 
            c(.3,.1), # 1 2
            c(.7,.5), # 1 3
            c(.9,.5), # 2 3
            c(0,0,0)) # 1 2 3 

# Hüsler-Reiss
Lambda <- cbind(c(0, 0.5, 0.5), 
                c(0.5, 0, 0.05), 
                c(0.5, 0.05, 0))


# max-linear
A <- cbind(c(0.2, 0.2, 0), 
           c(0.5, 0, 0),
           c(0.2, 0.2, 0.2),
           c(0.1, 0.6, 0.8))

X_sl <- rmev(n, d, param = par, model = "log")
X_al <- rmev(n, d, param = par, asy = asy, model = "alog")
X_hr <- rmev(n, d, sigma = Lambda, model = "hr")
X_ml <- SpatialExtremes::rmaxlin(n = n, dsgn.mat = A)

list(X_sl, X_al, X_hr, X_ml) %>%
  lapply(function(X) {
    X <- X %>%
      set_colnames(paste0("X", 1:d)) %>% 
      as.data.frame() %>%
      mutate(R = X1 + X2 + X3) %>%
      mutate(is_extreme = (R > quantile(R, 1 - k/n)))
    
    p <- ggtern(X, aes(X1, X2, X3)) +
      geom_mask() +
      geom_point(data = filter(X, !is_extreme), colour = "grey", size = 0.4, alpha = 0.2) +
      geom_point(data = filter(X, is_extreme), colour = "red", size = 0.4, alpha = 0.7) +
      xlab("") +
      ylab("") +
      Llab(expression(X[1])) +
      Tlab(expression(X[2])) +
      Rlab(expression(X[3])) +
      theme_light() +
      theme_hidegrid() +
      theme_hidelabels()
    return(p)
  }) %>%
  ggtern::grid.arrange(grobs = ., ncol = 2, nrow = 2)







