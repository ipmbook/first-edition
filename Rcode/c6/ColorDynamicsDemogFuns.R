G_z1z <- function(z1, z) dnorm(z1,mean=z,sd=0.025) 

s_z <- function(z) 0.85*exp(-35*(z+0.35)^2) + 0.9* exp(-35*(z-0.35)^2)

b_z <- function(z) 0.5*sqrt(s_z(z))*exp(-z/5); 

c_z1z <- function(z1,z) dnorm(z1,mean=z,sd=0.1); 

plot(s_z,-1,1); 


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
P_z1z <- function (z1, z) s_z(z) * G_z1z(z1, z)

## Define the fecundity kernel
F_z1z <- function (z1, z, m.par) b_z(z) * c_z1z(z1,z)

mk_K <- function(m,L, U) {
	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, P_z1z))
	F <- h * (outer(meshpts, meshpts, F_z1z))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}

out=mk_K(200,L=-1,U=1); K=out$K; 
n0=m0=rep(0,200); n0[200]=1; m0[1]=1; 
for(j in 1:5000) {
    n0=K%*%n0; m0=K%*%m0; 
    n0=n0/max(n0);   m0=m0/max(m0); 
}
matplot(out$meshpts,log(cbind(n0,m0)),type="l",col=c("blue","red")); 

