rm(list = ls())

# library(dplyr)
# library(lattice)
# library(gridExtra)
# library(Ternary)
# library(Rfast)
# library(stringr)

# Directory -------------------------------------------------------------------------

data_dir <- file.path("data/eva-data-challenge")
fn_dir <- file.path("R/eva-data-challenge")

# Source functions and load data ----------------------------------------------------

fn_files <- list.files(path = fn_dir, pattern = "*.R$", recursive = TRUE, full.names = TRUE)
sapply(fn_files, source)

# Load data -------------------------------------------------------------------------

Coputopia <- read.csv(file = file.path(data_dir, "Coputopia.csv"), header = TRUE)

X <- Coputopia %>% dplyr::select(Y1, Y2, Y3) %>% as.matrix() %>% exp() # Frechet(1) margins
colnames(X) <- paste0("X", seq_len(ncol(X)))
d <- ncol(X)

# C3 --------------------------------------------------------------------------------

k <- 500 # value used for submission

R <- apply(X, 1, pracma::Norm, p = 1)
R_kplus1 <- Rfast::nth(R, k = k + 1, descending = TRUE)
ext_ind <- which(R > R_kplus1)
piTheta_ext <- (X[ext_ind, ] / R_kplus1) %>% apply(1, euc_proj) %>% t()

face <- apply(piTheta_ext, 1, function(x) stringr::str_c(as.character(which(x > 0)), collapse = ""))
face_counts <- table(face)

A_hat <- max_linear_A_emp(X, k, alpha = 1)
A_hat_star <- (d / k) *  t(piTheta_ext)


## C3.1 ------------------------------------------------------------------------------

beta1 <- c(1, 2, 3)
u1 <- exp(6)
p1_hat <- max_linear_p_hat(A_hat_star, beta1, u1, alpha = 1)

## C3.2 ------------------------------------------------------------------------------

beta2 <- c(1, 2)
u2 <- exp(7)
p2_hat <- max_linear_p_hat(A_hat_star, beta2, u2, alpha = 1)

# Sensitivity to choice of k  --------------------------------------------------------

k_vals <- seq(from = 50, to = 1500, by = 25)
pk <- array(NA, dim = c(length(k_vals), 3))
colnames(pk) <- c("k", "p1", "p2")
pk[, "k"] <- k_vals

for (i in seq_along(k_vals)) {
  R_kplus1_k <- Rfast::nth(R, k = k_vals[i] + 1, descending = TRUE)
  ext_ind_k <- which(R > R_kplus1_k)
  piTheta_ext_k <- (X[ext_ind_k, ] / R_kplus1_k) %>% apply(1, euc_proj) %>% t()
  Ak <- (d / k_vals[i]) *  t(piTheta_ext_k)
  pk[i, "p1"] <- max_linear_p_hat(Ak, beta1, u1, alpha = 1)
  pk[i, "p2"] <- max_linear_p_hat(Ak, beta2, u2, alpha = 1)
}
pk <- as.data.frame(pk)


# Make plots ------------------------------------------------------------------------

## Ternary plot ----------------------------------------------------------------------

# library(ggtern)
# 
# p1 <- A_hat %>%
#   t() %>%
#   as.data.frame() %>%
#   ggtern(aes(X1, X2, X3)) +
#   geom_point(colour = "black", size = 0.8) +
#   xlab("") +
#   ylab("") +
#   Llab(expression(X[1])) +
#   Tlab(expression(X[2])) +
#   Rlab(expression(X[3])) +
#   theme_hidegrid() +
#   theme_hidelabels()
# 
# p2 <- A_hat_star %>%
#   t() %>%
#   as.data.frame() %>%
#   mutate(face_type = case_when(
#     (X1 > 0) + (X2 > 0) + (X3 > 0) == 3 ~ "3",
#     (X1 > 0) + (X2 > 0) + (X3 > 0) == 2 ~ "2",
#     .default = "1"
#   )) %>%
#   ggtern(aes(X1, X2, X3, colour = face_type)) +
#   geom_point(size = 0.8) +
#   xlab("") +
#   ylab("") +
#   Llab(expression(X[1])) +
#   Tlab(expression(X[2])) +
#   Rlab(expression(X[3])) +
#   labs(colour = expression("Subspace," ~ C[beta])) +
#   scale_colour_manual(labels = c(expression(bgroup("|",beta,"|") == 1),
#                                  expression(bgroup("|",beta,"|") == 2),
#                                  expression(bgroup("|",beta,"|") == 3)), values = c("blue", "red", "black")) +
#   theme_hidegrid() +
#   theme_hidelabels()
# 
# ggarrange(print(p2), print(p1), ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
# 
# detach(package:ggtern, unload = TRUE)
# 
# TernaryPlot(atip = "{1}", btip = "{2}", ctip = "{3}",
#             alab = "{1,3}", blab = "{1,2}", clab = "{2,3}", lab.font = 2, lab.offset = 0.1,
#             lab.cex = 1.4, tip.cex = 1.4,
#             grid.lines = 0)
# TernaryPoints((k / d) * t(A_hat), col = "black", pch = 4)
# TernaryText(c(1, 1, 1), labels = c("{1,2,3}"), col = "black", font = 2, cex = 1.4)
