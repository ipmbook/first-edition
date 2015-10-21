## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Set everything up
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(igraph)

## working directory must be set here, so the source()'s below run
path <- ifelse(.Platform$OS.type=="windows","c:/Repos/ipm_book/", "~/repos/ipm_book/")
setwd(path)

## source the utility functions
source("./Rcode/utilities/Standard Graphical Pars.R")

## source the vital rate functions for the model
source("./Rcode/c6/Ungulate Age Demog Funs.R")

## number of simulations
nsim <- 250

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Repeatedly simulate from the IBM 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## set up the ungulate IBM
set.seed(270875)
init.pop.size <- 500
n.yrs <- 400
m.par <- m.par.true

## run the ungulate IBM
sim.data.list <- list()
for (i in seq_len(nsim)) {
  source("./Rcode/c6/Ungulate Age Simulate IBM.R")
  rows.to.keep <- sample.int(nrow(sim.data), 3000)
  sim.data.list[[i]] <- sim.data[rows.to.keep,]  
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 3 - Fit the vital rate functions for IPMs
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m.par.szag.est <- m.par.sz.est <- list()

for (i in seq_len(nsim)) {

  sim.data <- sim.data.list[[i]]
  
  # 1. Size-age model
    
  m.par <- m.par.true
  
  m.par[c("surv.int","surv.z","surv.a")] <- 
    coef(glm(Surv ~ z + a, family = binomial, data = sim.data))
  
  m.par[c("grow.int","grow.z","grow.a")] <- 
    coef(mod.grow <- lm(z1 ~ z + a, data = sim.data))
  
  m.par[c("repr.int","repr.z","repr.a")] <- 
    coef(glm(Repr ~ z + a, family = binomial, data = subset(sim.data, a > 0))) # ! remove age 0 indiv. !
  
  m.par[c("recr.int","recr.a")] <- 
    coef(glm(Recr ~ a, family=binomial, data=sim.data))
  
  m.par[c("rcsz.int","rcsz.z")] <- 
    coef(mod.rcsz <- lm(Rcsz ~ z, data=sim.data))
  
  m.par["grow.sd"] <- summary(mod.grow)$sigma
  m.par["rcsz.sd"] <- summary(mod.rcsz)$sigma
  
  m.par.szag.est[[i]] <- m.par
  
  ## 2. Size-only model
  
  m.par <- m.par.true
  
  m.par[c("surv.int","surv.z","surv.a")] <- 
    c(coef(glm(Surv ~ z, family = binomial, data = sim.data)), 0)
  
  m.par[c("grow.int","grow.z","grow.a")] <- 
    c(coef(mod.grow <- lm(z1 ~ z, data = sim.data)), 0)
  
  m.par[c("repr.int","repr.z","repr.a")] <- 
    c(coef(glm(Repr ~ z, family = binomial, data = sim.data)), 0)  
  
  m.par[c("recr.int","recr.a")] <- 
    c(coef(glm(Recr ~ 1, family=binomial, data=sim.data)), 0)
  
  m.par[c("rcsz.int","rcsz.z")] <- 
    coef(mod.rcsz <- lm(Rcsz ~ z, data=sim.data))
  
  m.par["grow.sd"] <- summary(mod.grow)$sigma
  m.par["rcsz.sd"] <- summary(mod.rcsz)$sigma
      
  m.par.sz.est[[i]] <- m.par  
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 4 - Construct Leslie matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fa.list <- pa.list <- list()
for (i in seq_len(nsim)) {
  
  sim.data <- sim.data.list[[i]]  
  ## calc pa and fa 
  pa <- ma <- numeric(26)
  for(j in 1:12) {
    Xj <-  sim.data[sim.data$a==(j-1),] 
    pa[j] <- mean(Xj$Surv)
    ma[j] <- mean(Xj$Repr, na.rm=TRUE) * mean(Xj$Recr, na.rm=TRUE) 
  }
  Xj <- sim.data[sim.data$a>=12,] 
  pa[13:26] <- mean(Xj$Surv)
  ma[13:26] <- mean(Xj$Repr, na.rm=TRUE) * mean(Xj$Recr,na.rm=TRUE)
  ma[is.na(ma)] <- 0 
  fa <- 0.5 * pa * ma
  ## store 
  fa.list[[i]] <- fa 
  pa.list[[i]] <- pa 
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 5 - Function definitions: implementing the age-size kernel and iteration
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 
## Functions to build IPM kernels P(a), F(a),
## 

## Define the survival kernel
P_z1z <- function (z1, z, a, m.par) {
  return( s_z(z, a, m.par) * g_z1z(z1, z, a, m.par) )
}

## Define the reproduction kernel
F_z1z <- function (z1, z, a, m.par) {
  return( s_z(z, a, m.par) * pb_z(z, a, m.par) * (1/2) * pr_z(a, m.par) * c_z1z(z1, z, m.par) )
}

## 
## Functions to implement the IPM
## 

## Calculate the mesh points, mesh width and store with upper/lower bounds and max age
mk_intpar <- function(m, L, U, M) {
  h <- (U - L) / m
  meshpts  <-  L + ((1:m) - 1/2) * h
  na <- M + 2
  return( list(meshpts = meshpts, M = M, na = na, h = h, m = m) )
}

## Build the list of age/process specific kernels + store the integration parameters in the same list
mk_age_IPM <- function(i.par, m.par) {
  within(i.par, {
    F <- P <- list()
    for (ia in seq_len(na)) {
      F[[ia]] <- outer(meshpts, meshpts, F_z1z, a = ia-1, m.par = m.par) * h
      P[[ia]] <- outer(meshpts, meshpts, P_z1z, a = ia-1, m.par = m.par) * h
    }
    rm(ia)
  })
}

## 
## Functions to iterate the age-size IPM
##

## iterate 'forward' one time step to project dynamics. 'x' is the list of age-specific size dists
r_iter <- function(x, na, F, P) {
  xnew <- list(0)
  for (ia in seq_len(na)) {
    xnew[[1]] <- (xnew[[1]] + F[[ia]] %*% x[[ia]])[,,drop=TRUE]
  }
  for (ia in seq_len(na-1)) {
    xnew[[ia+1]] <- (P[[ia]] %*% x[[ia]])[,,drop=TRUE]
  }
  xnew[[na]] <- xnew[[na]] + (P[[na]] %*% x[[na]])[,,drop=TRUE]
  return(xnew)
}

##  list to vector function for arpack 
matmul.r <- function(x.vec, IPM.sys) {
  with(IPM.sys, {
    x.list <- list()
    for (ia in seq_len(na)) x.list[[ia]] <- x.vec[seq.int(ia*m-m+1, ia*m)]
    unlist(r_iter(x.list, na, F, P))
  })
}

## 
## Functions to iterate the age-size IPM next gen kernel
##

## function to multiply by the next generation matrix, R
matmulR <- function(x, IPM.sys, n.iter=100) {
  with(IPM.sys, {
    xtot <- rep(0, length(x))
    for (ia in seq_len(n.iter)) {
      if (ia < na) {Pnow <- P[[ia]]; Fnow <- F[[ia]]} 
      else         {Pnow <- P[[na]]; Fnow <- F[[na]]}
      xtot <- Fnow %*% x + xtot; x <- Pnow %*% x 
    }
    xtot
  })
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 6 - Compare predictions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# integration parameters
i.par <- mk_intpar(m = 100, M = 20, L = 1.6, U = 3.7)

## 1. age-size model

lam.szag <- R0.szag <- GT.szag <- numeric()
for (i in seq_len(nsim)) {
  # make an initial state vector
  nt0 <- with(i.par, lapply(seq_len(na), function(ia) rep(0, m)))
  nt0[[1]] <- with(i.par, rep(1 / m, m))
  # 
  IPM.sys <- mk_age_IPM(i.par, m.par.szag.est[[i]])
  # estimate lambda 
  lam.szag[i] <- Re(arpack(matmul.r, extra=IPM.sys, complex=FALSE, 
                           options=list(n=with(IPM.sys, m*na), nev=1))$val[1])  
  ## estimate R0
  R0.szag[i] <- arpack(matmulR, extra=IPM.sys, complex=FALSE, 
                       options=list(n=IPM.sys$m, nev=1))$val[1]
  ## generation time
  GT.szag[i] <- log(R0.szag[i]) / log(lam.szag[i])
}


## 2. size-only model

lam.sz <- R0.sz <- GT.sz <- numeric()
for (i in seq_len(nsim)) {
  # just use the existing function for the age-size model to ket the required kernels
  IPM.sys <- mk_age_IPM(i.par, m.par.sz.est[[i]])
  IPM.sys$F[[1]] <- IPM.sys$F[[2]] # need to update because F_0 == 0
  F.sz <- IPM.sys$F[[1]]; P.sz <- IPM.sys$P[[1]]; K.sz <- F.sz + P.sz
  # pop growth rate and stable size distribution
  lam.sz[i] <- Re(eigen(K.sz)$val[1])
  # R0 via the eigen function
  N.sz <- solve(diag(dim(K.sz)[1]) - P.sz)
  R0.sz[i] <- Re(eigen(F.sz %*% N.sz, only.values=TRUE) $ values[1])
  # generation time
  GT.sz[i] <- log(R0.sz[i]) / log(lam.sz[i])
}

## 3. age-only model

lam.ag <- R0.ag <- GT.ag <- numeric()
for (i in seq_len(nsim)) {
  fa <- fa.list[[i]]
  pa <- pa.list[[i]]
  ##  
  L <- matrix(0, 26, 26) 
  L[1,] <- fa
  for(j in 1:25) L[j+1, j] <- pa[j]
  ## lambda, R0, generation time
  lam.ag[i] <- abs(eigen(L)$values[1])
  la <- c(1, cumprod(pa))
  R0.ag[i]  <- sum(la * c(fa, 0))
  GT.ag[i]  <- log(R0.ag[i]) / log(lam.ag[i])
}


hist(lam.szag)
hist(R0.szag)
hist(GT.szag)

hist(lam.sz)
hist(R0.sz)
hist(GT.sz)

hist(lam.ag)
hist(R0.ag)
hist(GT.ag)


plot(lam.szag, lam.sz)
abline(0,1)
plot(lam.szag, lam.ag)
abline(0,1)

plot(R0.szag, R0.sz)
abline(0,1)
plot(R0.szag, R0.ag)
abline(0,1)

plot(GT.szag, GT.sz)
abline(0,1)
plot(GT.szag, GT.ag)
abline(0,1)

postscript("c6/figures/SzAg-vs-Sz-vs-Ag.eps", 
           width=8, height=4, horizontal=FALSE, paper="special")
set_graph_pars("panel2")
# plot(lam.szag, lam.sz, xlim = c(0.99, 1.05), ylim = c(0.99, 1.05), 
#      pch = 20, 
#      xlab = "Growth rate, Age-size model", ylab = "Growth rate, Size-only model")
# abline(0, 1, lty = 2)
# add_panel_label("a")
# plot(lam.szag, lam.ag, xlim = c(0.99, 1.05), ylim = c(0.99, 1.05),
#      pch = 20, 
#      xlab = "Growth rate, Age-size model", ylab = "Growth rate, Age-only model")
# abline(0, 1, lty = 2)
# add_panel_label("b")
plot(GT.szag, GT.ag, xlim = c(5.85, 6.5), ylim = c(5.45, 6.1),
     pch = 20, 
     xlab = "Generation time, Age-size model", ylab = "Generation time, Age-only model")

abline(0, 1, lty = 2)
add_panel_label("c")
plot(GT.szag, GT.sz, xlim = c(5.85, 6.5), ylim = c(9.5, 12.5),
     pch = 20, 
     xlab = "Generation time, Age-size model", ylab = "Generation time, Size-only model")

abline(0, 1, lty = 2)
add_panel_label("d")
dev.off()

