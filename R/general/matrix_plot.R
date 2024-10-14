library(colorspace)
library(abind)
library(lattice)

levelplotmatrix <- function(M, colpal, lower = -max(abs(M)) - 0.001, upper = max(abs(M)) + 0.001, x_labels = TRUE, y_labels = TRUE,...) {
  
  if (is.list(M)) M <- abind(M, along = 3)
  
  if (length(dim(M)) == 2) Mperm <- t(M[nrow(M):1, ])
  if (length(dim(M)) == 3) Mperm <- aperm(M[nrow(M):1, , ], c(2, 1, 3))
  
  scales <- list(x = list(labels = parse(text = rownames(Mperm))),
                 y = list(labels = parse(text = colnames(Mperm))))
    
  if (!x_labels) scales$x <- list(draw = FALSE)
  if (!y_labels) scales$y <- list(draw = FALSE)
  
  lattice::levelplot(Mperm, 
                     col.regions = colpal, 
                     at = lattice::do.breaks(c(lower, upper), length(colpal)), 
                     as.table = TRUE, 
                     xlab = "", ylab = "", scales = scales,
                     ...)
}


tpdm_colpal <- function(n) sequential_hcl(n = n, "Viridis")

div_colpal <- function(n) diverging_hcl(n = n, "Blue-Red")

plot_tpdm <- function(Sigma, n_col = 13, ...) levelplotmatrix(M = Sigma, colpal = tpdm_colpal(n_col), lower = 0, ...)

plot_tpdm_eigen <- function(U, n_col = 13, ...) levelplotmatrix(M = U, colpal = div_colpal(n_col), ...)


# to customise the strip formatting, e.g.
# strip = lattice::strip.custom(bg = strip.bg, par.strip.text = list(cex = 0.7))