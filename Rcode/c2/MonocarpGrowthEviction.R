#### Check for eviction in monocarp growth kernel


## Assumes R's working directory is set to the location of this file,
## and Monocarp Demog Funs is in the same location
source("Monocarp Demog Funs.R")
source("../utilities/Standard Graphical Pars.R")

set_graph_pars("panel4")

zvals <- seq(-2.65, 4, length = 200)
tiny <- 0.001
g1 <- G_z1z(zvals, min(zvals), m.par.true)
g1[g1 < tiny] <- NA
g2 <- G_z1z(zvals, mean(zvals), m.par.true)
g2[g2 < tiny] <- NA
g3 <- G_z1z(zvals, max(zvals), m.par.true)
g3[g3 < tiny] <- NA

matplot(zvals, cbind(g1, g2, g3), type = "l", lty = 1, lwd = c(1, 2, 3), 
        col = "black", xlab = "Final size z'", ylab = "Probabity density", 
        ylim = c(0, 1))

add_panel_label("a")

## Compute and plot fraction going to the wrong place
WrongPlace <- function(z, U) {
  fac1 <- s_z(z, m.par.true) * (1 - p_bz(z, m.par.true))
  fac2 <- integrate(function(x) G_z1z(x, z, m.par.true), U, Inf)$value
  return(fac1 * fac2)
}
zvals <- Wvals <- seq(-2.65, 4.5, length = 200)
for (j in seq_along(zvals)) Wvals[j] <- WrongPlace(zvals[j], 4)
plot(zvals, Wvals, type = "l", lty = 1, lwd = 2, col = "black", xlab = "Initial size z", 
     ylab = "Fraction wrongfully evicted", ylim = c(0, 1.1 * max(Wvals)))

add_panel_label("b")

## Recompute fraction for U=4.5
for (j in seq_along(zvals)) Wvals[j] <- WrongPlace(zvals[j], 4.5)
points(zvals, Wvals, type = "l", lty = 2, lwd = 2)

## Plot seed production b_z as a function of size
par(xaxs = "r")
zvals <- seq(3, 4.5, length = 25)
Nvals <- b_z(zvals, m.par.true)
plot(zvals, Nvals, type = "l", xlab = "Parent size z", ylab = "Seed production", 
     ylim = c(0, max(Nvals)), lwd = 2)
points(c(4, 6), rep(b_z(4, m.par.true), 2), type = "l", lty = 2, lwd = 2)

add_panel_label("c")

## Compute dominant eigenvalue as a function of the upper limit
Uvals <- lambdas <- lambda.bd <- seq(4, 5.5, length = 20)
for (j in seq_along(Uvals)) {
  IPM <- mk_K_ceiling(150, m.par.true, L = -2.65, U = Uvals[j], U1 = 10)
  lambdas[j] <- Re(eigen(IPM$K)$values[1])
}
par(yaxs = "i", xaxs = "r")
plot(Uvals, lambdas, type = "l", xlab = "Upper limit U", ylab = "Dominant eigenvalue", 
     ylim = c(1.02, 1.08), lwd = 2)

add_panel_label("d")


for (j in seq_along(Uvals)) {
  IPM.bd <- mk_K_ceiling(150, m.par.true, L = -2.65, U = Uvals[j], U1 = 4)
  lambda.bd[j] <- Re(eigen(IPM.bd$K)$values[1])
}
points(Uvals, lambda.bd, type = "l", lty = 2, lwd = 2)

dev.copy2eps(file = "../../figures/c2/MonocarpGrowthEviction.eps")
