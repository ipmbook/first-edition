########## Dirt-simple kernel 
G_z1z <- function(z1, z) dnorm(z1,mean=1+0.4*z,sd=0.3) 
s_z <- function(z) 0.25 + sqrt(1+z)/4; 
b_z <- function(z) sqrt(1+z); 
c_z1z <- function(z1,z) dnorm(z1,mean=0.3,sd=0.1); 

######### Functions to build the iteration matrices 
P_z1z <- function (z1, z) s_z(z) * G_z1z(z1, z)
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

## Get lambda, w, and R0 from eigen by computing R = F(I-P)^{-1}
out=mk_K(500,-0.5,3.5); K=out$K; P=out$P; F=out$F; meshpts=out$meshpts;
eK=eigen(K); 
lambda=Re(eK$values[1]); w=eK$vectors[,1]; w=abs(w)/sum(abs(w)); 
plot(meshpts,w); 

I=diag(rep(1),nrow(P)); 
R=F%*%solve(I-P); R0=Re(eigen(R)$values[1]); R0; 

## get R0 without inverting I-P, using arpack 
PF=list(P=P,F=F); 
matmul <- function(x, PF) {
    P=PF$P; F=PF$F; 
    xtot=x; newx=x; 
    for(j in 1:50) {newx=P%*%newx; xtot=xtot+newx}
    return(F%*%xtot); 
}    
out=arpack(matmul, extra=PF, options=list(n=nrow(P), nev=1))
R0.arpack=out$values[1]; R0.arpack; 


