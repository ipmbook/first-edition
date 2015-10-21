## Analysis of the age-structured ungulate example.

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Set everything up
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(igraph)

## working directory must be set here, so the source()'s below run
path=ifelse(.Platform$OS.type=="windows","c:/Repos/ipm_book/","~/repos/ipm_book/"); 
setwd(path); 

## source the utility functions
source("./Rcode/utilities/Standard Graphical Pars.R")

## source the vital rate functions for the model
source("./Rcode/c6/Ungulate Age Demog Funs.R")

## run the ungulate IBM
# set.seed(270875)
init.pop.size <- 500
n.yrs <- 400
m.par <- m.par.true
source("./Rcode/c6/Ungulate Age Simulate IBM.R")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Estimation
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## take a random sample of 3000 observations 
rows.to.keep <- sample.int(nrow(sim.data), 3000)
sim.data <- sim.data[rows.to.keep,]

## sort by size and print a sample to the screen
sim.data <- sim.data[order(sim.data$z),]
round(sim.data[sample(1:nrow(sim.data),12),][-c(3,4,6,7),],2)

##
## Fit the vital rate functions
##

mod.surv.glm <- glm(Surv ~ z + a, family = binomial, data = sim.data)
summary(mod.surv.glm)

grow.data <- subset(sim.data, !is.na(z1))
mod.grow <- lm(z1 ~ z + a, data = grow.data)
summary(mod.grow)

repr.data <- subset(sim.data, Surv==1 & a>0)
mod.repr <- glm(Repr ~ z + a, family = binomial, data = repr.data)
summary(mod.repr)

recr.data <- subset(sim.data, Surv==1 & Repr==1)
mod.recr <- glm(Recr ~ a, family=binomial, data=recr.data)
summary(mod.recr)

rcsz.data <- subset(sim.data, !is.na(Rcsz))
mod.rcsz <- lm(Rcsz ~ z, data=rcsz.data)
summary(mod.rcsz)

##  
## Store the estimated parameters
##

m.par.est <- c(
  ## survival
  surv      = coef(mod.surv.glm),
  ## growth 
  grow      =  coef(mod.grow),
  grow.sd   =  summary(mod.grow)$sigma,
  ## reproduce or not
  repr      =  coef(mod.repr),
  ## recruit or not
  recr      =  coef(mod.recr),
  ## recruit size
  rcsz      =  coef(mod.rcsz),
  rcsz.sd   =  summary(mod.rcsz)$sigma)

names(m.par.est) <- names(m.par.true)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 3 - Function definitions: implementing the kernel and iteration
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

## 'left' iteration one time step. 'x' is the list of age-specific size dists
l_iter <- function(x, na, F, P) {
  xnew <- list(0)
  for (ia in seq_len(na-1)) {
    xnew[[ia]] <- (x[[ia+1]] %*% P[[ia]]  + x[[1]] %*% F[[ia]])[,,drop=TRUE]
  }
  xnew[[na]] <- (x[[na]] %*% P[[na]] + x[[1]] %*% F[[na]])[,,drop=TRUE]
  return(xnew)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 4 - Compare truth and estimated model predictions of lambda, and calculate w(z) and v(z), using
## iteration for 100 steps
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# integration parameters
i.par <- mk_intpar(m = 100, M = 20, L = 1.6, U = 3.7)

# make an initial state vector
nt0 <- with(i.par, lapply(seq_len(na), function(ia) rep(0, m)))
nt0[[1]] <- with(i.par, rep(1 / m, m))

# estimate lambda with the true parameters...
IPM.sys <- mk_age_IPM(i.par, m.par.true)
IPM.sim <- with(IPM.sys, {
  x <- nt0
  for (i in seq_len(100)) {
    x1 <- r_iter(x, na, F, P)
    lam <- sum(unlist(x1))
    x <- lapply(x1, function(x) x / lam)
  }
  list(lambda = lam, x = x)
})
IPM.sim$lambda


# ... and using those estimated from the IBM 'data'
IPM.sys <- mk_age_IPM(i.par, m.par.est)
IPM.sim <- with(IPM.sys, {
  x <- nt0
  for (i in seq_len(100)) {
    x1 <- r_iter(x, na, F, P)
    lam <- sum(unlist(x1))
    x <- lapply(x1, function(x) x / lam)
  }
  list(lambda = lam, x = x)
})
IPM.sim$lambda

# store the stable age*size distribution and renormalise it to be a genuine density functions
wz <- IPM.sim$x
wz <- lapply(wz, function(wa) wa/i.par$h)

# reproductive value calculations
IPM.sim <- with(IPM.sys, {
  x <- nt0
  for (i in seq_len(100)) {
    x1 <- l_iter(x, na, F, P)
    lam <- sum(unlist(x1))
    x <- lapply(x1, function(x) x / lam)
  }
  list(lambda = lam, x = x)
})
vz <- IPM.sim$x
vz <- lapply(vz, function(va) va/i.par$h)

# plot the 
postscript(file="./c6/figures/SoaySizeAgeAndRV.eps", w=8.5, h=4.5,
           horizontal = FALSE, onefile = FALSE, paper = "special")
set_graph_pars(ptype = "panel2")
# stable distribution
z.dist <- wz[[1]]
plot(i.par$meshpts, z.dist, type="n", cex.lab=1.2,
     xlab="Log mass, z", ylab=expression(italic(w[a](z))), xlim=c(1.8, 3.7))
for (ia in seq_len(i.par$na)) lines(i.par$meshpts, wz[[ia]])
add_panel_label("a")
# reproductive value
z.dist <- vz[[1]]
plot(i.par$meshpts, z.dist, type="n", cex.lab=1.2,
     xlab="Log mass, z", ylab=expression(italic(v[a](z))), xlim=c(1.8, 3.7))
for (ia in seq_len(i.par$na)) lines(i.par$meshpts, vz[[ia]])
add_panel_label("b")
# clean up
dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 5 - Population growth rate calculated in different ways
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##
## 1. Just use iteration and a stopping rule based on the error - pretty fast
##

## calc lambda and stable size dist by iteration using 'right' multiplication
domEig_r <- function(IPM.sys, tol) {
  with(IPM.sys, {
    # initial state
    wt <- lapply(seq_len(na), function(ia) rep(0, m))
    wt[[1]] <- rep(1 / m, m)
    # iterate until converges in lambda
    lam <- 1; qmax <- 10*tol
    while (qmax > tol) {
      wt1 <- r_iter(wt, na, F, P)
      qmax <- sum(abs(unlist(wt1) - lam * unlist(wt)))
      lam <- sum(unlist(wt1)); wt <- lapply(wt1, function(x) x/lam)
    }
    list(lambda = lam, wz = wt)
  })
}

## calc lambda and reproductive value dist by iteration using 'left' multiplication
domEig_l <- function(IPM.sys, tol) {
  with(IPM.sys, {
    # initial state
    vt <- lapply(seq_len(na), function(ia) rep(0, m))
    vt[[na]] <- rep(1 / m, m)
    # iterate until converges in lambda
    lam <- 1; qmax <- 10*tol
    while (qmax > tol) {
      vt1 <- l_iter(vt, na, F, P)
      qmax <- sum(abs(unlist(vt1) - lam * unlist(vt)))
      lam <- sum(unlist(vt1))
      vt <- lapply(vt1, function(x) x/lam)
    }
    list(lambda = lam, vz = vt)
  })
}

IPM.sim.r <- domEig_r(IPM.sys, 1e-6)
IPM.sim.l <- domEig_l(IPM.sys, 1e-6)
IPM.sim.l$lambda; IPM.sim.r$lambda # these must match

##
## 2. using arpack - almost as fast
##

matmul.r <- function(x.vec, IPM.sys) {
  with(IPM.sys, {
    x.list <- list()
    for (ia in seq_len(na)) x.list[[ia]] <- x.vec[seq.int(ia*m-m+1, ia*m)]
    unlist(r_iter(x.list, na, F, P))
  })
}
eigs.r <- arpack(matmul.r, extra=IPM.sys, complex=FALSE, options=list(n=with(IPM.sys, m*na), nev=1))
Re(eigs.r$val[1])

##
## 3. using the big 'Goodman' matrix with eigen - very slow
##

## take the list of IPM kernels and use these to construct a single 'Goodman' iteration matrix
IPM_list_to_mat <- function(IPM.sys) {
  with(IPM.sys, {
    big.mat <- matrix(0, na*m, na*m)
    isize  <- list()
    for (ia in seq.int(1, na)) isize[[ia]] <- seq.int(ia*m-m+1, ia*m)
    for (ia in seq.int(1, na)) big.mat[ isize[[1]], isize[[ia]] ] <- F[[ia]]
    for (ia in seq.int(1, na-1)) big.mat[ isize[[ia+1]], isize[[ia]] ] <- P[[ia]]
    big.mat[ isize[[na]], isize[[na]]] <- P[[na]]
    big.mat
  })
}

if(FALSE){
 IPM.big.mat <- IPM_list_to_mat(IPM.sys)
 Re(eigen(IPM.big.mat)$val[1])

## 4. use the 'Goodman' matrix with arpack - faster than eigen, but still much slower than using list form

matmul.G <- function(x, IPM.big.mat) {
  return(IPM.big.mat %*% x)
}
arpack(matmul.G, extra=IPM.big.mat, complex=FALSE, options=list(n=with(IPM.sys, m*na), nev=1))$val[1]
}

## 
## 5. comparison of run times (takes quite a long time to run...) 
## 

library(microbenchmark)
n <- with(IPM.sys, m*na)
if(FALSE)
microbenchmark(times = 10,
               Re(eigen(IPM.big.mat)$val[1]),
               arpack(matmul.G, extra=IPM.big.mat, complex=FALSE, options=list(n=n, nev=1))$val[1],
               arpack(matmul.r, extra=IPM.sys, complex=FALSE, options=list(n=n, nev=1))$val[1],
               domEig_r(IPM.sys, 1e-9)$lam)
}
# Output...
# Unit: milliseconds
# min        lq          median      uq          max         neval
# 81145.5690 81252.19094 81366.23157 81474.01956 81502.0172  10
# 2706.4029  2759.95091  2870.75050  2943.77356  3006.9845   10
# 183.7572   194.30038   202.96817   214.05179   223.1214    10
# 84.0223    85.47455    95.13607    97.28803    121.5931    10

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 6 - Store anything we'll use later
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract and store lambda for use later
lam <- IPM.sim.r$lambda
# extract the stable age-size dist and reproductive value vector (+ normalise and convert to list form)
wz <- IPM.sim.r$wz
wz <- lapply(wz, function(wz) wz / (sum(wz) * i.par$h))
vz <- IPM.sim.l$vz
vz <- lapply(vz, function(vz) vz / (sum(vz) * i.par$h))
# initial size distribution of new recuits
c.init <- wz[[1]]/sum(wz[[1]])

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 7 - R0 using the kernel components
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

## use arpack to compute the leading eigenvalue of the next gen matrix
R0.arpack <- arpack(matmulR, extra=IPM.sys, complex=FALSE, options=list(n=IPM.sys$m, nev=1))$val[1]
R0.arpack

# Does this match the R0 we get by direct iteration with the next generation matrix...
sum(matmulR(c.init, IPM.sys))
# ... similar, but not the same, because lambda != 1

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 8 - Calculate the age-specific rates
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m.par <- m.par.est
IPM.sys <- mk_age_IPM(i.par, m.par)
IPM.sim <- domEig_r(IPM.sys, 1e-6)

##
## Age-specific vital rates 
##

## function to take a set of age-specific fecundities and survival probs and return a 
## Leslie matrix (not quite a Leslie matrix as we include the 'greybeard' class as well)
mk_leslie <- function(fa, pa) {
  mdim <- length(fa)
  P <- F <- matrix(0, mdim, mdim)
  F[1,] <- fa
  P[seq(1, mdim^2-mdim, by=mdim+1)+1] <- pa
  P[mdim, mdim] <- pa[length(pa)]
  return(list(P = P, F = F))
}

# stable distribution of new recruits
c.init <- IPM.sim$wz[[1]]
c.init <- c.init/sum(c.init)
# 'e' vector
e <- matrix(1, nrow=1, ncol=length(c.init))

# age-specific vital rates
na.mat <- 100
la <- fa <- numeric()
Pa <- diag(length(c.init))
for (ia in seq.int(1, na.mat)) {
  if (ia < i.par$na) {
    Pnow <- IPM.sys$P[[ia]]
    Fnow <- IPM.sys$F[[ia]]
  } else {
    Pnow <- IPM.sys$P[[i.par$na]]
    Fnow <- IPM.sys$F[[i.par$na]]    
  }
  la[ia] <- sum((e %*% Pa) * c.init)
  fa[ia] <- sum((e %*% Fnow %*% Pa) * c.init)
  Pa <- Pnow %*% Pa
}
fa <- fa / la
pa <- la[seq.int(2, na.mat)] / la[seq.int(1, na.mat-1)]
round(la, 3); round(pa, 3); round(fa, 3)

# build the age-structured transition matrix
lmat <- mk_leslie(fa, pa)
K.lmat <- lmat$F + lmat$P

# lambda
lam.lmat <- Re(eigen(K.lmat)$values[1])
lam.lmat

# fundamental matrix
N.lmat <- solve(diag(na.mat) - lmat$P)

# mean lifespan
sum(la)

# R0 three equivalent ways
eigen(lmat$F %*% N.lmat)$values[1]
(lmat$F %*% N.lmat)[1,1]
sum(la * fa)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 9 - Generation time
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# generation time using the R0 and lam from the model
log(R0.arpack) / log(lam)

# generation time using the mean age-specific vital rates
log(sum(la * fa)) / log(lam.lmat) # log(R0)/log(lam)
sum(seq_along(fa) * fa * la) / sum(fa * la) # mean age of mothers at offspring production

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 10 - Fit the wrong model -> size structure only
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##
## Fit and store the parameters for vital rate functions ignoring age
##

m.par.noage <- m.par.true

m.par.noage[c("surv.int","surv.z","surv.a")] <- 
  c(coef(glm(Surv ~ z, family = binomial, data = sim.data)), 0)

m.par.noage[c("grow.int","grow.z","grow.a")] <- 
  c(coef(lm(z1 ~ z, data = sim.data)), 0)

m.par.noage[c("repr.int","repr.z","repr.a")] <- 
  c(coef(glm(Repr ~ z, family = binomial, data = sim.data)), 0)

m.par.noage[c("recr.int","recr.a")] <- 
  c(coef(glm(Recr ~ 1, family=binomial, data=sim.data)), 0)

##
## Build the corresponding IPM, compute lambda and some basic life cycle calculations
##

# just use the existing function for the age-size model to ket the required kernels
IPM.sys.noage <- mk_age_IPM(i.par, m.par.noage)
IPM.sys.noage$F[[1]] <- IPM.sys.noage$F[[2]] # need to update because F_0 == 0
F.noage <- IPM.sys.noage$F[[1]]
P.noage <- IPM.sys.noage$P[[1]]  
K.noage <- F.noage + P.noage

# pop growth rate and stable size distribution
lam.noage <- Re(eigen(K.noage)$val[1])
wt.noage <- abs(Re(eigen(K.noage)$vec[,1]))
wt.noage <- wt.noage / (sum(wt.noage) * i.par$h)

# R0 via the eigen function
N.noage <- solve(diag(length(wt.noage)) - P.noage)
R0.noage <- Re(eigen(F.noage %*% N.noage, only.values=TRUE) $ values[1])
R0.noage

# generation time
log(R0.noage) / log(lam.noage)

##
## Age-specific vital rates + associated life cycle calculations
##

e <- matrix(1, nrow=1, ncol=length(wt.noage))
c.init.noage <- (F.noage %*% wt.noage)[,,drop=TRUE]
c.init.noage <- c.init.noage / sum(c.init.noage)

# age-spcific vital rates
na.mat <- 100
la.noage <- fa.noage <- numeric()
Pa <- diag(length(wt.noage))
for (ia in seq.int(1, na.mat)) {
  la.noage[ia] <- sum((e %*% Pa) * c.init.noage)
  fa.noage[ia] <- sum((e %*% F.noage %*% Pa) * c.init.noage)
  Pa <- P.noage %*% Pa
}
fa.noage <- fa.noage / la.noage
pa.noage <- la.noage[seq.int(2, na.mat)] / la.noage[seq.int(1, na.mat-1)]
round(la.noage, 3)
round(pa.noage, 3)
round(fa.noage, 3)

# sanity check -> do we get the same lambda if we work with the mean age-specific rates?
lmat.noage <- mk_leslie(fa.noage, pa.noage)
N.lmat.noage <- solve(diag(na.mat) - lmat.noage$P)
Re(eigen(lmat.noage$F + lmat.noage$P)$values[1])
# ...yes

# mean lifespan
sum(la.noage)

# R0
eigen(lmat.noage$F %*% N.lmat.noage)$values[1]
sum(la.noage * fa.noage)

# generation time
log(sum(la.noage * fa.noage)) / log(lam.noage)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 11 - Plot to summarise the difference between right and wrong models
##   with regard to age-specific survival/fecundity, and short-term predictions 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot(0:25, la.noage[1:26], pch=20, cex=0.8, col="red", ylim = c(0.015, 1),
     xlab="Age", ylab="Survivorship / Fecundity",type="b")
points(0:25, fa.noage[1:26], pch=4 , cex=0.8, col="red",type="b")
points(0:25, la[1:26],       pch=20, cex=0.8, col="blue",type="b")
points(0:25, fa[1:26],       pch=4 , cex=0.8, col="blue",type="b")
legend("topright", bty= "n", cex = 0.8,
       legend=c("Survivorship (size-only IPM)","Fecundity (size-only IPM)",
                "Survivorship (age-size IPM)","Fecundity (age-size IPM)"),
       col=rep(c("red","blue"), each= 2),
       pch=rep(c(20,4), times = 2)
)


##### Plot for the book, on effects of ignoring age 

pa = (la[2:26]/la[1:25])
pa.noage = (la.noage[2:26]/la.noage[1:25])pch=c(16,1))

graphics.off();  
dev.new(width=8,height=4);
set_graph_pars("panel2"); par(mgp=c(2.5,1,0)) 
matplot(0:20,cbind(pa.noage[1:21],pa[1:21]), type="o",lty=c(1,2),pch=c(16,1),
xlab="Age",col="black",ylab="Survival p(a) to age (a+1)");
legend("bottomleft",c("Size-only model","Age-size model"),lty=c(1,2),pch=c(16,1),bty="n")
add_panel_label("a")


### Iterate starting with 100 4 year old females, and 100 8 year old females
### in each case with the stable size-distribution for their age 
wInit4=wInit8=wz; 
for(i in 1:length(wz)) wInit4[[i]] = wInit8[[i]]=0*wz[[i]]	
wInit4[[4]]=100*wz[[4]]/sum(wz[[4]]); 
wInit8[[8]]=100*wz[[8]]/sum(wz[[8]])

N4 = N8 = rep(100,11); na=length(wz);
out<- with(IPM.sys, {
  x4 <- wInit4; x8 <- wInit8;
  for (i in 2:11) {
    x4 <- r_iter(x4, na, F, P);
    x8 <- r_iter(x8,na,F,P);
    N4[i] <- sum(unlist(x4));
    N8[i] <- sum(unlist(x8))
  }
  list(N4=N4,N8=N8)
})

matplot(0:10,cbind(out$N4,out$N8),ylim=c(50,200),xlab="Time (years)", ylab="Total number of females",
pch=16,type="o",lty=c(1,2),col="black",bty="n")
legend("topright",c("Age 4","Age 8"),lty=c(1,2),bty="n")
add_panel_label("b")
dev.copy2eps(file="c6/figures/UngulateSizeOnly.eps")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 12 - Perturbation analysis
## Bonus! Not in the book...
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 
## Kernel sensitivity and elasticity
## 

calc_K_sens <- function(wz, vz) {
  na <- IPM.sys$na; h <- IPM.sys$h
  denom <- h * sum(unlist(vz) * unlist(wz))
  within(list(), {
    F <- P <- list()
    for (ia in seq_len(na)) {
      F[[ia]] <- outer(vz[[1]], wz[[ia]], "*") / denom
      if (ia < na) {
        P[[ia]] <- outer(vz[[ia+1]], wz[[ia]], "*") / denom
      } else {
        P[[na]] <- outer(vz[[na]], wz[[na]], "*") / denom
      }
    }
    rm(ia)
  })
}

calc_K_elas <- function(K.sens, IPM.sys, lambda) {
  na <- IPM.sys$na; h <- IPM.sys$h
  within(list(), {
    F <- P <- list()
    for (ia in seq_len(na)) {
      F[[ia]] <- K.sens$F[[ia]] * (IPM.sys$F[[ia]] / h) / lambda
      if (ia < na) {
        P[[ia]] <- K.sens$P[[ia]] * (IPM.sys$P[[ia]] / h) / lambda
      } else {
        P[[na]] <- K.sens$P[[na]] * (IPM.sys$P[[na]] / h) / lambda
      }
    }
    rm(ia)
  })
}

K.sens <- calc_K_sens(wz, vz)
K.elas <- calc_K_elas(K.sens, IPM.sys, lam)
# get the elasticities integrated over size 
K.elas.age.F <- sapply(K.elas$F, sum) * i.par$h^2
K.elas.age.P <- sapply(K.elas$P, sum) * i.par$h^2
# get the elasticities summed over age 
K.elas.size <- Reduce("+", K.elas$F) + Reduce("+", K.elas$P)

## 
## Function sensitivity and elasticity example
## 

# 1. survival function sensitivity / elasticity

dK_by_ds_F_z1z <- function(a, m.par, meshpts) {
  pert_kern <- function (z1, z, a, m.par) {
    pb_z(z, a, m.par) * (1/2) * pr_z(a, m.par) * c_z1z(z1, z, m.par)
  }
  outer(meshpts, meshpts, pert_kern, a = a, m.par = m.par)
}

dK_by_ds_P_z1z <- function(a, m.par, meshpts) {
  pert_kern <- function (z1, z, a, m.par) {
    g_z1z(z1, z, a, m.par)
  }
  outer(meshpts, meshpts, pert_kern, a = a, m.par = m.par)
}

integrand.F <- 
  lapply(seq_len(i.par$na), 
         function (ia) K.sens$F[[ia]] * dK_by_ds_F_z1z(ia-1, m.par, i.par$meshpts))
integrand.P <- 
  lapply(seq_len(i.par$na),
         function (ia) K.sens$P[[ia]] * dK_by_ds_P_z1z(ia-1, m.par, i.par$meshpts))

s.sens.F <- lapply(integrand.F, function(x) apply(x, 2, sum) * i.par$h)
s.sens.P <- lapply(integrand.P, function(x) apply(x, 2, sum) * i.par$h)

elas.term <- lapply(seq_len(i.par$na), 
                    function (ia) s_z(i.par$meshpts, ia-1, m.par) / lam)

s.elas.F <- mapply("*", elas.term, s.sens.F, SIMPLIFY=FALSE)
s.elas.P <- mapply("*", elas.term, s.sens.P, SIMPLIFY=FALSE) 

sum(unlist(s.elas.F)) * i.par$h + sum(unlist(s.elas.P)) * i.par$h

plot(i.par$meshpts, s.elas.F[[1]], type="n", ylim=c(0,0.10))
for (ia in seq_len(9)+1) lines(i.par$meshpts, s.elas.F[[ia]], lwd=2-0.15*ia)
plot(i.par$meshpts, s.elas.P[[1]], type="n", ylim=c(0,0.55))
for (ia in seq_len(10)) lines(i.par$meshpts, s.elas.P[[ia]], lwd=2-0.15*ia)

# 2. growth function sensitivity / elasticity

dK_by_dg_z1z <- function(a, m.par, meshpts) {
  pert_kern <- function (z1, z, a, m.par) s_z(z, a, m.par)
  outer(meshpts, meshpts, pert_kern, a = a, m.par = m.par)
}

g.sens <-
  lapply(seq_len(i.par$na),
         function (ia) K.sens$P[[ia]] * dK_by_dg_z1z(ia-1, m.par, i.par$meshpts))

elas.term <- lapply(seq_len(i.par$na), 
                    function (ia) outer(i.par$meshpts, i.par$meshpts, g_z1z, a = ia-1, m.par = m.par) / lam)

g.elas <- mapply("*", elas.term, g.sens, SIMPLIFY=FALSE)




