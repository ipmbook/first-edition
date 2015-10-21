## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Ungulate IBM to illustrate the construction of an IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

## working directory must be set here, so the source()'s below run
setwd("~/Repos/ipm_book/")
## run the utility functions
source("./Rcode/utilities/Standard Graphical Pars.R")
## run the ungulate IBM, fit demographic models & build the IPM
source("./Rcode/c2/Ungulate Calculations.R")
## set the working directory to our figures location
setwd("~/Repos/ipm_book/c4/figures")

## copy model parameters 
m.par <- m.par.true

## build the IPM associated with the  
IPM.sys <- mk_K(nBigMatrix, m.par, min.size, max.size)

## extract the mesh points and the kernel
meshpts <- IPM.sys$meshpts
h <- diff(meshpts[1:2])
K <- IPM.sys$K
P <- IPM.sys$P
F <- IPM.sys$F

## compute the stable size distribution
IPM.eig.sys <- eigen(K)
w.z <- abs(Re(IPM.eig.sys$vectors[,1]))
w.z <- w.z / sum(w.z)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## 1. mean and variance of lifespan
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## define the e vector which we will use for summing down columns
e <- matrix(1, nrow = 1, ncol = nBigMatrix)
## compute the fundamental operator
N <- solve(diag(1, nBigMatrix) - P)
## mean lifespan, given initial size
eta.bar <- (e %*% N)[,,drop=TRUE]

## sensitivty of eta
g.z1z0 <- outer(meshpts, meshpts, g_z1z, m.par) * h
eta.sens <- N * (eta.bar %*% g.z1z0)[,,drop=TRUE]

## fourth term (sensitivity of eta^2)
eta2.sens <- t(eta.bar * t(eta.sens))
plot(apply(eta2.sens, 2, sum) * h, type="l") 

## second term
eta.bar.N <- eta.sens %*% N

## first term
eta.bar.H.N <- matrix(NA, nBigMatrix, nBigMatrix)
for (i.z0 in seq_along(meshpts)) {
  # compute the kernel H
  H <- outer(g.z1z0[,i.z0], N[i.z0,], function(g.z1, N.z) g.z1 * N.z)
  # now compute the double integral associated with first term
  eta.bar.H.N[i.z0,] <- eta.bar %*% N %*% H
}

eta.var.sens <- 2 * eta.bar.H.N + 2 * eta.bar.N - eta.sens - 2 * eta2.sens
plot(apply(eta.var.sens, 2, sum) * h, type="l")