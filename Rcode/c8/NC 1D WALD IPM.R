## Neubert-Caswell IPM in 1D with WALD kernel for dispersal distance.
## A 'hybrid' species: Oenothera demography and C. nutans dispersal,
## plus some totally made-up height dependence of dispersal 
setwd("~/Repos/ipm_book/Rcode/c8"); 

################################################### 
#  1. Load the Oenothera IPM for local demography
#  and add some made-up height dependence in dispersal     
###################################################
source("../c2/Monocarp Demog Funs.R")
L <- (-2.65); U <- 4; # size range for IPM 

# Check: we get the right eigenvalue 
out <- mk_K(100,m.par.true,L,U);
eigen(out$K)$values[1]; 


# Effect of parent height on dispersal (made-up)
# In WALD, tau is proportional to seed release
# height H, and mu is proportional to H^2. 
# We imagine that H is proportional to z, and give
# a z=3 individual WALD parameter values close
# to those estimated by Skarpaas & Shea for C. nutans

tau_z <- function(z) {h <- pmax(z/3,0.1); 7*h*h}
mu_z <- function(z)  {h <- pmax(z/3,0.1); 3*h}; 

############################################################ 
## 2. Functions to compute the mfg of the marginalized
## WALD distribution as a function of its parameters. 
## We pass these as entries in a parameter vector w.par 
########################################################### 

## Function to compute the WALD mgf
WALDmgf <- function(s,w.par) {
    mu <- w.par["r.mu"];
    tau <- w.par["r.tau"];
    t1 <- (tau/mu); 
    t2 <- 2*(mu^2)*s/tau; 
    mgf <- exp(t1*(1-sqrt(1-t2)));
    return(mgf)
}    
    
## Function to compute the marginalize WALD mgf
margWALDmgf <- function(s,w.par) {
  (1/pi)*integrate(function(q) WALDmgf(s*cos(q),w.par),0,pi)$value;
}
# Reality check: mgf and marginalized mgf should have value 1 at s=0
w.par <- c(3,7); names(w.par) <- c("r.mu","r.tau")
WALDmgf(0,w.par); margWALDmgf(0,w.par)

# Maximum value of s for which the WALD mgf is finite.  
s.max <- function(w.par) {
  mu <- w.par["r.mu"];
  tau <- w.par["r.tau"];
  return(tau/(2*mu*mu)); 
}  

################################################### 
# 3. Function to compute the transformed kernel 
###################################################
Hs <- function(s,m,m.par,L,U) {
	out<-mk_K(m,m.par,L,U);
	P <- out$P; Fs <- out$F; zvals <- out$meshpts;
	for(j in 1:m) {
		mu <- mu_z(zvals[j]);
		tau <- tau_z(zvals[j])
		w.par <- c(mu,tau); 
		names(w.par) <- c("r.mu","r.tau")
	    Fs[,j] <- Fs[,j]*margWALDmgf(s,w.par)
	}
	return(P+Fs)     
}

################################################ 
# 4. Function to compute wave speeds c(s)
################################################
cs <- function(s,m,m.par,L,U) {
	M <- Hs(s,m,m.par,L,U) 
	L1 = abs(eigen(M)$values[1]); 
    return((1/s)*log(L1)) 
   }
cs = Vectorize(cs,"s"); 
plot(function(s) cs(s,100,m.par.true,L,U),0.05,0.2);

## Find the asymptotic wave speed c*(s) 
out=optimize(cs,lower=0.05,upper=0.16,m=100,m.par=m.par.true,L=L,U=U); 
cat("Wave speed cstar =",out$objective,"\n"); 
cat("Wave shape parameter s =",out$minimum,"\n");


