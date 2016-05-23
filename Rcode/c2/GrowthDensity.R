###### Plot to illustrate growth densities
Gyx <- function(y, x) {
  dnorm(y, mean = 1.8 + 0.7 * x, sd = 0.1 + (x/100)^0.3)
}

yvals <- seq(1, 10, length = 1000)

tiny <- 0.001
g1 <- Gyx(yvals, 2)
g1[g1 < tiny] <- NA
g2 <- Gyx(yvals, 5)
g2[g2 < tiny] <- NA
g3 <- Gyx(yvals, 8)
g3[g3 < tiny] <- NA

#### Open the graphics window
dev.new(w = 7, h = 4)

par(cex.axis = 1.4, bty = "l", cex.lab = 1.4, yaxs = "i")
par(mgp = c(2.5, 1, 0))
matplot(yvals, cbind(g1, g2, g3), type = "l", lty = 1, lwd = c(1, 2, 3), 
        col = "black", xlab = "Final size z'", ylab = "Probabity density", 
        ylim = c(0, 1))
text(2.25, 0.7, "G(z',2)", cex = 1.3)
text(4.35, 0.6, "G(z',5)", cex = 1.3)
text(8.6, 0.5, "G(z',8)", cex = 1.3)

#### Save the plot as EPS and PDF
# dev.copy2eps(file = "../../figures/c2/GrowthDensity.eps")
# dev.copy2pdf(file = "../../figures/c2/GrowthDensity.eps")
