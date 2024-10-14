library(lattice)
library(patchwork)
library(RColorBrewer)
library(colorspace)
library(plot3D)

levelplotmatrix <- function(M, colpal, lim.low = -max(abs(M)), lim.up = max(abs(M)) + 0.001, draw.axes = FALSE, ...) {
  n.colpal <- length(colpal)
  breaks <- do.breaks(c(lim.low, lim.up), n.colpal)
  strip.bg <- rgb(190, 190, 190, max = 255, alpha = 160)
  if (is.matrix(M)) {
    Mperm <- t(M[nrow(M):1, ])
  } else if (is.array(M)) {
    Mperm <- aperm(M[nrow(M):1, , ], c(2, 1, 3))
  } else stop("M must be a matrix or an array.")
  lattice::levelplot(Mperm,
                     col.regions = colpal, at = breaks,
                     colorkey = list(space = "right"), 
                     scales = list(draw = draw.axes),
                     as.table = TRUE, 
                     xlab = "", ylab = "",
                     strip = strip.custom(bg = strip.bg, par.strip.text = list(cex = 0.7)),
                     ...)
}

#spatial plots
# plot zvals at coordinates [loc$x, loc$y]
# "real" -> divergent colour scheme ranging from -max(abs(zvals)) to max(abs(zvals))
# otherwise -> sequential colour scheme ranging from 0 to max(abs(zvals))
plot.spatial.vals <- function(zvals, loc, type = "real", main = "", zmax = max(abs(zvals)) + 0.001) {
  n.colpal <- 9
  if (type == "real") {
    col.pal <- c(brewer.pal(n.colpal,"Blues")[n.colpal:2], brewer.pal(n.colpal,"Reds")[2:n.colpal])
    zmin <- -zmax
  } else {
    col.pal <- brewer.pal(n.colpal, "Reds")
    zmin <- 0
  }
  temp <- seq(zmin - 0.001, zmax + 0.001, length.out = length(col.pal) + 1)
  # map("worldHires", "uk",
  #     col = alpha("olivedrab", 0.5), fill = TRUE, bg = FALSE, boundary = TRUE,
  #     xlim = c(-3.4, -1.4), ylim = c(52.8, 55.7),
  #     mar = c(1, 1, 1, 1))
  # map.axes(xaxt = 'n', yaxt = 'n', ann = TRUE)
  range.x <- max(loc$x) - min(loc$x)
  xlim <- c(min(loc$x) - 0.05 * range.x, max(loc$x) + 0.05 * range.x)
  range.y <- max(loc$y) - min(loc$y)
  ylim <- c(min(loc$y) - 0.05 * range.y, max(loc$y) + 0.05 * range.y)
  plot(NULL, xaxt = 'n', yaxt = 'n', axes = TRUE, ann = TRUE, xlim = xlim, ylim = ylim, main = main, xlab = "", ylab = "")
  grid(lwd = 1)
  for (j in seq_along(zvals)) {
    temp.j <- sum(zvals[j] > temp)
    points(loc$x[j], loc$y[j], pch = 20, cex = 2, col = col.pal[temp.j])
  }
  colkey(side = 4, col = col.pal, clim = c(min(temp), max(temp)), labels = F, breaks = round(temp, 2), add = T, cex.axis = 0.8, dist = 0, length = 0.95)
}