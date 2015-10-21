rm(list=ls(all=TRUE))

library(igraph)

## working directory must be set here, so the source()'s below run
path=ifelse(.Platform$OS.type=="windows","c:/Repos/ipm_book/","~/repos/ipm_book/"); 
setwd(path); 

## source the utility functions
source("./Rcode/utilities/Standard Graphical Pars.R")

## source the vital rate functions for the model
source("./Rcode/c6/Ungulate Age Demog Funs.R")

m.par <- m.par.true

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - IPM functions, exactly the same as in Ungulate Age Calculations.R
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
P_z1z <- function (z1, z, a, m.par) {s_z(z, a, m.par) * g_z1z(z1, z, a, m.par)}
F_z1z <- function (z1, z, a, m.par) {
  return( s_z(z, a, m.par) * pb_z(z, a, m.par) * (1/2) * pr_z(a, m.par) * c_z1z(z1, z, m.par)
}

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

## Iterate 'forward' one time step to project dynamics. 'x' is the list of age-specific size dists
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

# integration parameters
i.par <- mk_intpar(m = 100, M = 20, L = 1.6, U = 3.7)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - IPM functions, exactly the same as in Ungulate Age Calculations.R
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# make an initial state vector
nt0 <- with(i.par, lapply(seq_len(na), function(ia) rep(0, m)))
nt0[[1]] <- with(i.par, rep(1 / m, m))

# estimate lambda with the true parameters...
IPM.sys <- mk_age_IPM(i.par, m.par.true)
IPM.sim <- with(IPM.sys, {
  x <- nt0
  for (i in seq_len(500)) {
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

############################
#  Real work starts here
###########################
graphics.off();  dev.new(width=8,height=4);
set_graph_pars("panel2"); 

w4=wz[[4]]/sum(wz[[4]]); w8= wz[[8]]/sum(wz[[8]]);

matplot(i.par$meshpts,cbind(w4,w8),col=c("black"),lty=c(1,2),type="l",xlim=c(2.7,3.7),
xlab="Size z", ylab="Frequency"); 
legend("topright",c("Age 4","Age 8"),lty=c(1,2))

wInit4=wInit8=wz; 
for(i in 1:length(wz)) wInit4[[i]] = wInit8[[i]]=0*wz[[i]]	
wInit4[[4]]=100*w4; wInit8[[8]]=100*w8

N4 = N8 = rep(100,11); na=length(wz)
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

 matplot(0:10,cbind(out$N4,out$N8),ylim=c(50,200),xlab="Years", ylab="Total number of females",
 pch=1,type="o",lty=c(1,2),col="black",bty="n")

