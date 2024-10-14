levelplotmatrix <- function(M, colpal, lower = -max(abs(M)) - 0.001, upper = max(abs(M)) + 0.001,...) {
  
  if (is.list(M)) M <- abind(M, along = 3)
  
  if (length(dim(M)) == 2) Mperm <- t(M[nrow(M):1, ])
  if (length(dim(M)) == 3) Mperm <- aperm(M[nrow(M):1, , ], c(2, 1, 3))
  
  lattice::levelplot(Mperm, 
                     col.regions = colpal, 
                     at = lattice::do.breaks(c(lower, upper), length(colpal)), 
                     as.table = TRUE, 
                     xlab = "", ylab = "",
                     ...)
}


tpdm_colpal <- function(n) colorRampPalette(c('white', 'black'))(n)

div_colpal <- function(n) colorRampPalette(c('red', 'blue'))(n)

plot_tpdm <- function(Sigma, n_col = 30, ...) levelplotmatrix(M = Sigma, colpal = tpdm_colpal(n_col), lower = 0, ...)

plot_tpdm_eigen <- function(U, n_col = 15, ...) levelplotmatrix(M = U, colpal = div_colpal(n_col), ...)


# to customise the strip formatting, e.g.
# strip = lattice::strip.custom(bg = strip.bg, par.strip.text = list(cex = 0.7))