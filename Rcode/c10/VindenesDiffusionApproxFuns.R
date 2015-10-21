################################################################################
# Code copied from Ecological Archives E092-093-S1 by Vindenes et al. (2011) 
################################################################################
#============================================================================
# DIFFUSION APPROXIMATION 
#============================================================================
# The following functions perform simulations of a diffusion process with 
# density-independent  dynamics, using the expected growth rate (lambda), 
# the environmental variance (envar) and the demographic variance (demvar). 
# Letting y denote population size on log scale, the infinitesimal mean is 
# log(lambda)-1/(2*lambda^2)*(envvar+demvar/exp(y)) and the infinitesimal 
# variance is 1/(2*lambda^2)(envvar+demvar/exp(y)). The parameter "y0" is 
# the initial population size on log scale. Other parameters are "t.max" 
# (the maximum number of time steps), "delta.t" (the size of each time interval, 
# the process is stored only at integer times),  "n.sim" (the number of 
# desired realizations), and "b" (an upper limit to the process (on log scale)).  
# Each realization is stored in a separate row in the returned matrix. 
# The process is defined as extinct if y<=0  (population size exp(y)<=1).

mu.nu <- function(y, envar, demvar, lambda){
    mu <- log(lambda)-(1/(2*lambda^2))*(envar+demvar/exp(y))
    nu <- (1/lambda^2)*(envar+demvar/exp(y))
    mu[y<=0] <- 0
    nu[y<=0] <- 0
    return(list("mu"=mu, "nu"=nu))
    }
    
mu.nu.Ito <- function(y, envar, demvar, lambda){
    mu <- lambda-1-0.5*(envar+demvar/exp(y)); 
    nu <- (envar+demvar/exp(y))
    mu[y<=0] <- 0
    nu[y<=0] <- 0
    return(list("mu"=mu, "nu"=nu))
    }
    
diffusion <- function(y0, t.max=1000, delta.t=.1, n.sim=1000, b=Inf,...) {
  yy <- matrix(NA,ncol=(t.max+1),nrow=n.sim)
  tt <- seq(0,length=(t.max+1))
  yy[,1] <- y <- y0
  tt[1] <- t <- 0
  j <- 1
  for (i in 1:(t.max/delta.t)) {
    mn  <- mu.nu(y=y,...)
    y <- y + mn$mu*delta.t + rnorm(n.sim)*sqrt(mn$nu*delta.t)
    y[y>=b] <- b
    y[y<=0] <- 0
    if (i%%(1/delta.t)==0) {
    	j <- j + 1
      	yy[,j] <- y
    	}
  	}
  	return(yy) 
	}

diffusion.Ito <- function(y0, t.max=1000, delta.t=.1, n.sim=1000, b=Inf,...) {
  yy <- matrix(NA,ncol=(t.max+1),nrow=n.sim)
  tt <- seq(0,length=(t.max+1))
  yy[,1] <- y <- y0
  tt[1] <- t <- 0
  j <- 1
  for (i in 1:(t.max/delta.t)) {
    mn  <- mu.nu.Ito(y=y,...)
    y <- y + mn$mu*delta.t + rnorm(n.sim)*sqrt(mn$nu*delta.t)
    y[y>=b] <- b
    y[y<=0] <- 0
    if (i%%(1/delta.t)==0) {
    	j <- j + 1
      	yy[,j] <- y
    	}
  	}
  	return(yy) 
	}

