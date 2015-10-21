rm(list = ls(all = TRUE))

set.seed(53241986)

source("Ungulate Demog Funs.R")
source("../utilities/Standard Graphical Pars.R")
source("../utilities/MatrixImage.R")

graphics.off()
par(mfrow = c(2, 2), mar = c(0, 1, 2, 0))
nBigMatrix <- 80
IPM.true <- mk_K(nBigMatrix, m.par.true, 1.6, 3.6)
K <- IPM.true$K
meshpts <- IPM.true$meshpts
persp(meshpts, meshpts, t(K), theta = 30, xlab = "Size at time t", ylab = "Size at time t+1", 
      zlab = "K(z',z)", col = "grey80", d = 2, expand = 0.85)
add_panel_label("a")

Bnorm <- function(z1, z, mean, sd) {
  dnorm(z1, mean[1], sd) * dnorm(z, mean[2], sd)
}

C <- outer(meshpts, meshpts, Bnorm, mean = c(2, 2), sd = 0.1)
fac <- 0.1 * max(K)/max(C)

C <- C * fac
persp(meshpts, meshpts, t(K + C), theta = 30, xlab = "Size at time t", 
      ylab = "Size at time t+1", zlab = "K(z',z)", col = "grey80", d = 2, 
      expand = 0.85)
add_panel_label("b")

C2 <- outer(meshpts, meshpts, Bnorm, mean = c(2, 2), sd = 0.05)
C2 <- C2 * fac
persp(meshpts, meshpts, t(K + C2), theta = 30, xlab = "Size at time t", 
      ylab = "Size at time t+1", zlab = "K(z',z)", col = "grey80", d = 2, 
      expand = 0.85)
add_panel_label("c")

C3 <- outer(meshpts, meshpts, Bnorm, mean = c(2, 2), sd = 0.02)
C3 <- C3 * fac
persp(meshpts, meshpts, t(K + C3), theta = 30, xlab = "Size at time t", 
      ylab = "Size at time t+1", zlab = "K(z',z)", col = "grey80", d = 2, 
      expand = 0.85)
add_panel_label("d")

dev.copy2eps(file = "../../figures/c2/PerturbUngulate.eps")

Bnorm <- function(z1, z, mean, sd) {
  dnorm(z1, mean[1], sd) * dnorm(z, mean[2], sd)
}

C <- outer(meshpts, meshpts, Bnorm, mean = c(2, 2), sd = 0.1)
fac <- 0.1 * max(K)/max(C)

C <- C * fac
persp(meshpts, meshpts, t(K + C), theta = 30, xlab = "Size at time t", 
      ylab = "Size at time t+1", zlab = "K(z',z)", col = "grey80", d = 2, 
      expand = 0.85)
add_panel_label("b")

C2 <- outer(meshpts, meshpts, Bnorm, mean = c(2, 2), sd = 0.05)
C2 <- C2 * fac
persp(meshpts, meshpts, t(K + C2), theta = 30, xlab = "Size at time t", 
      ylab = "Size at time t+1", zlab = "K(z',z)", col = "grey80", d = 2, 
      expand = 0.85)
add_panel_label("c")

C3 <- outer(meshpts, meshpts, Bnorm, mean = c(2, 2), sd = 0.02)
C3 <- C3 * fac
persp(meshpts, meshpts, t(K + C3), theta = 30, xlab = "Size at time t", 
      ylab = "Size at time t+1", zlab = "K(z',z)", col = "grey80", d = 2, 
      expand = 0.85)
add_panel_label("d")

dev.copy2eps(file = "../../figures/c2/PerturbUngulate.eps")

