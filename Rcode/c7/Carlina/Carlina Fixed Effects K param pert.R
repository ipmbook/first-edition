## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fit and run fixed and mixed effects effects Carlina stochastic IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
library(lme4)
library(MCMCglmm)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c2",sep="")); 

source("../utilities/Standard Graphical Pars.R");

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c7/Carlina",sep="")); 

#load yearly parameter vector 
load("Yearly parameters.Rdata")

source("Carlina Demog Funs DI.R") 

 
#####################################################################
#Stochastic perturbation analysis
#####################################################################

nBigMatrix <- 100
n.est <- 20000
n.runin <- 500
minsize <- 1.5
maxsize <- 5
n.years <-20

stoc_pert_analysis<-function(params,n.est,n.runin,C.t,C.t.mean){
	
	year.i <- sample(1:n.years,n.est+1,replace=TRUE)

	K.year.i <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.K<-mk_K(nBigMatrix,params[,i],minsize,maxsize)
		K.year.i[i,,] <- year.K$K
	}

	h <- year.K$h; 
	meshpts <- year.K$meshpts
	
#Calculate mean kernel, v and w

	mean.kernel <- apply(K.year.i,2:3,mean)

	w <- Re(eigen(mean.kernel)$vectors[,1]); 
	v <- Re(eigen(t(mean.kernel))$vectors[,1]);

	# scale eigenvectors <v,w>=1 
	w <- abs(w)/sum(h*abs(w))
	v <- abs(v)
	v <- v/(h*sum(v*w))
    cat(h*sum(v*w)," should = 1","\n")

#Esimate Lambda s
#initialize variables	

	nt<-rep(1/nBigMatrix,nBigMatrix)
	rt.V <- rt.N <- rep(NA,n.est)
	
#Iterate model

	for (year.t in 1:n.est){
		if(year.t%%10000==0) cat("iterate: ", year.t,"\n");

			
		#iterate model with year-specific kernel
		nt1<-K.year.i[year.i[year.t],,] %*% nt
	
		sum.nt1<-sum(nt1)
		
		#Calculate log growth rates  
		
		rt.V[year.t] <- log(sum(nt1*v)/sum(nt*v))
		rt.N[year.t] <- log(sum(nt1)/sum(nt))
		nt <- nt1 / sum.nt1  
	
	}

Ls <- exp(mean(rt.V))
	
	
### Get wt and Rt time series ###
	wt<-matrix(1/nBigMatrix, nrow=n.est+1, ncol=nBigMatrix);
	
	for (i in 1:n.est) {
		
		K             <- K.year.i[year.i[i],,]
		wt[i+1,]  <-K %*% wt[i,]
		wt[i+1,]  <-wt[i+1,]/sum(wt[i+1,]);
		if(i%%10000==0) cat("wt ",i,"\n")

	}

	
### Get vt time series ###
	vt<-matrix(1/nBigMatrix, nrow=n.est+1, ncol=nBigMatrix);
	for (i in (n.est+1):2) {
		
		K           <- K.year.i[year.i[i],,]
		vt[i-1,]  <- vt[i,] %*% K
		vt[i-1,]  <- vt[i-1,]/sum(vt[i-1,]);
		if(i%%10000==0) cat("vt  ",i,"\n")

	}

elas.s <- matrix(0,nBigMatrix,nBigMatrix)
elas.s.mean <- matrix(0,nBigMatrix,nBigMatrix)

for (year.t in n.runin:(n.est-n.runin)) {

		#standard calculations needed for the various formulae

	vt1.wt                       <- outer(vt[year.t+1,],wt[year.t,],FUN="*")
	vt1.C.wt                    <- vt1.wt * C.t[year.i[year.t],,]  
	        
	 vt1.C.wt.mean        <- vt1.wt * C.t.mean[year.i[year.t],,] 
	 	
	K           <- K.year.i[year.i[year.t],,]
	
	vt1.K.wt    <- sum(vt[year.t+1,] * (K %*% wt[year.t,]))

		#calculation of the standard elasticities
	
		elas.s            <-elas.s + (vt1.C.wt) / vt1.K.wt;
		elas.s.mean <-elas.s.mean + (vt1.C.wt.mean) / vt1.K.wt;
		
}

 elas.s             <- elas.s/(n.est-2*n.runin+1)
 elas.s.mean  <- elas.s.mean/(n.est-2*n.runin+1)
 
 
 return(list(meshpts=year.K$meshpts, h=h, elas.s=elas.s, elas.s.mean=elas.s.mean, mean.kernel=mean.kernel, Ls=Ls))

}

########################################################################
#Let's do the intercept of the probability of flowering function
########################################################################

#Select the parameters to use
params.to.use <- m.par.est

#Calculate meshpts and h for evaluating the perturbation kernels
year.K      <-  mk_K(nBigMatrix,params.to.use[,1],minsize,maxsize)
meshpts <-  year.K$meshpts
h              <-  year.K$h

#First calculate the mean and sd of beta0 and perturbation kernels

beta.0.mean <- mean(params.to.use["flow.int",])
beta.0.sd       <- sd(params.to.use["flow.int",])

Ct_z1z <- function(z1,z,m.par){
		return(s_z(z, m.par) * ( m.par["p.r"] * b_z(z, m.par) * c_0z1(z1, m.par) - G_z1z(z1, z, m.par)) * 
		(1/(1+exp(m.par["flow.int"]+m.par["flow.z"]*z))) * p_bz(z, m.par) *
		 m.par["flow.int"])
	 }
	
C.pert <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.C <- h * (outer(meshpts, meshpts, Ct_z1z, m.par = params.to.use[,i]))
		C.pert[i,,] <- year.C
	}

Ct_z1z_mean <- function(z1,z,m.par){
		return(s_z(z, m.par) * ( m.par["p.r"] * b_z(z, m.par) * c_0z1(z1, m.par) - G_z1z(z1, z, m.par)) *
		(1/(1+exp(m.par["flow.int"]+m.par["flow.z"]*z)))  * p_bz(z, m.par) *
		 beta.0.mean)
	}

C.pert.mean <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.C <- h * (outer(meshpts, meshpts, Ct_z1z_mean, m.par = params.to.use[,i]))
		C.pert.mean[i,,] <- year.C
	}

pert.K <- stoc_pert_analysis(params.to.use, n.est, n.runin, C.pert, C.pert.mean)

elas.s            <- sum(pert.K$elas.s)
elas.s.mean <- sum(pert.K$elas.s.mean)
elas.s.sd       <- elas.s-elas.s.mean
cat("Stochastic elasticity ",elas.s,"\n")
cat("Stochastic elasticity mean ",elas.s.mean,"\n")
cat("Stochastic elasticity sd ",elas.s.sd,"\n")
cat("Stochastic sensitivity  mean ",pert.K$Ls*elas.s.mean/beta.0.mean,"\n")
cat("Stochastic sensitivity  sd ",pert.K$Ls*elas.s.sd/beta.0.sd,"\n")

##################################################################
#Code for time invariant parameters

stoc_pert_analysis_invariant<-function(params,n.est,n.runin,K.t.deriv,K.t.deriv2){
	
	year.i <- sample(1:n.years,n.est+1,replace=TRUE)

	K.year.i <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.K<-mk_K(nBigMatrix,params[,i],minsize,maxsize)
		K.year.i[i,,] <- year.K$K
	}

	h <- year.K$h; 
	meshpts <- year.K$meshpts
	
#Calculate mean kernel, v and w

	mean.kernel <- apply(K.year.i,2:3,mean)

	w <- Re(eigen(mean.kernel)$vectors[,1]); 
	v <- Re(eigen(t(mean.kernel))$vectors[,1]);

	# scale eigenvectors <v,w>=1 
	w <- abs(w)/sum(h*abs(w))
	v <- abs(v)
	v <- v/(h*sum(v*w))
    cat(h*sum(v*w)," should = 1","\n")

#Esimate Lambda s
#initialize variables	

	nt<-rep(1/nBigMatrix,nBigMatrix)
	rt.V <- rt.N <- rep(NA,n.est)
	
#Iterate model

	for (year.t in 1:n.est){
		if(year.t%%10000==0) cat("iterate: ", year.t,"\n");

			
		#iterate model with year-specific kernel
		nt1<-K.year.i[year.i[year.t],,] %*% nt
	
		sum.nt1<-sum(nt1)
		
		#Calculate log growth rates  
		
		rt.V[year.t] <- log(sum(nt1*v)/sum(nt*v))
		rt.N[year.t] <- log(sum(nt1)/sum(nt))
		nt <- nt1 / sum.nt1  
	
	}

Ls <- exp(mean(rt.V))
	
	
### Get wt and Rt time series ###
	wt<-matrix(1/nBigMatrix, nrow=n.est+1, ncol=nBigMatrix);
	
	for (i in 1:n.est) {
		
		K             <- K.year.i[year.i[i],,]
		wt[i+1,]  <-K %*% wt[i,]
		wt[i+1,]  <-wt[i+1,]/sum(wt[i+1,]);
		if(i%%10000==0) cat("wt ",i,"\n")

	}

	
### Get vt time series ###
	vt<-matrix(1/nBigMatrix, nrow=n.est+1, ncol=nBigMatrix);
	for (i in (n.est+1):2) {
		
		K           <- K.year.i[year.i[i],,]
		vt[i-1,]  <- vt[i,] %*% K
		vt[i-1,]  <- vt[i-1,]/sum(vt[i-1,]);
		if(i%%10000==0) cat("vt  ",i,"\n")

	}

Exp.1 <- 0
Exp.2 <- 0

for (year.t in n.runin:(n.est-n.runin)) {

		#standard calculations needed for the various formulae

	vt1.K.t.deriv.wt       <- sum(vt[year.t+1,] * (K.t.deriv[year.i[year.t],,] %*% wt[year.t,]))
	vt1.K.t.deriv2.wt    <- sum(vt[year.t+1,] * (K.t.deriv2[year.i[year.t],,] %*% wt[year.t,]))
	K                                <- K.year.i[year.i[year.t],,]
	vt1.K.wt                   <- sum(vt[year.t+1,] * (K %*% wt[year.t,]))

		#calculation the expectations 
	
		Exp.1            <- Exp.1 + vt1.K.t.deriv2.wt / vt1.K.wt;
		Exp.2            <- Exp.2 + (vt1.K.t.deriv.wt / vt1.K.wt)^2;		
}

 Exp.1              <- Exp.1 /(n.est-2*n.runin+1)
 Exp.2              <- Exp.2 /(n.est-2*n.runin+1)
 
 sens.s             <- (Ls/2)*(Exp.1-Exp.2)
 return(list(meshpts=year.K$meshpts, h=h, sens.s=sens.s, Ls=Ls))

}

#################################################
#Calculate K derivative kernels
#################################################


Kt_partial_A <- function(z1,z,m.par){
		return(s_z(z, m.par) * p_bz(z, m.par) * m.par["p.r"] * c_0z1(z1, m.par) * b_z(z, m.par))	
		 }
	
Kt.partial.A <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.C <- h * (outer(meshpts, meshpts, Kt_partial_A, m.par = params.to.use[,i]))
		Kt.partial.A[i,,] <- year.C
	}

#Note the derivative kernels are the same here but we're written the code as if they
#were different so it's easier to alter for other parameters.

Kt_partial_A2 <- function(z1,z,m.par){
		return(s_z(z, m.par) * p_bz(z, m.par) * m.par["p.r"] * c_0z1(z1, m.par) * b_z(z, m.par))
	}

Kt.partial.A2 <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.C <- h * (outer(meshpts, meshpts, Kt_partial_A2, m.par = params.to.use[,i]))
		Kt.partial.A2[i,,] <- year.C
	}

pert.K <- stoc_pert_analysis_invariant(params.to.use, n.est, n.runin, Kt.partial.A, Kt.partial.A2)

pert.K$sens.s





#old code for checking...

stoc._pert_analysis_tmp<-function(params,n.est,n.runin,C.t,C.t.mean){
	
	year.i <- sample(1:n.years,n.est+1,replace=TRUE)

	K.year.i <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.K<-mk_K(nBigMatrix,params[,i],minsize,maxsize)
		K.year.i[i,,] <- year.K$K
	}

	h <- year.K$h; 
	meshpts <- year.K$meshpts
	
#Calculate mean kernel, v and w

	mean.kernel <- apply(K.year.i,2:3,mean)

	w <- Re(eigen(mean.kernel)$vectors[,1]); 
	v <- Re(eigen(t(mean.kernel))$vectors[,1]);

	# scale eigenvectors <v,w>=1 
	w <- abs(w)/sum(h*abs(w))
	v <- abs(v)
	v <- v/(h*sum(v*w))
    cat(h*sum(v*w)," should = 1","\n")

#Esimate Lambda s
#initialize variables	

	nt<-rep(1/nBigMatrix,nBigMatrix)
	rt.V <- rt.N <- rep(NA,n.est)
	
#Iterate model

	for (year.t in 1:n.est){
		if(year.t%%10000==0) cat("iterate: ", year.t,"\n");

			
		#iterate model with year-specific kernel
		nt1<-K.year.i[year.i[year.t],,] %*% nt
	
		sum.nt1<-sum(nt1)
		
		#Calculate log growth rates  
		
		rt.V[year.t] <- log(sum(nt1*v)/sum(nt*v))
		rt.N[year.t] <- log(sum(nt1)/sum(nt))
		nt <- nt1 / sum.nt1  
	
	}

Ls <- mean(rt.V)
	
	
### Get wt and Rt time series ###
	wt<-matrix(1/nBigMatrix, nrow=n.est+1, ncol=nBigMatrix);
	
	for (i in 1:n.est) {
		
		K             <- K.year.i[year.i[i],,]
		wt[i+1,]  <-K %*% wt[i,]
		wt[i+1,]  <-wt[i+1,]/sum(wt[i+1,]);
		if(i%%10000==0) cat("wt ",i,"\n")

	}

	
### Get vt time series ###
	vt<-matrix(1/nBigMatrix, nrow=n.est+1, ncol=nBigMatrix);
	for (i in (n.est+1):2) {
		
		K           <- K.year.i[year.i[i],,]
		vt[i-1,]  <- vt[i,] %*% K
		vt[i-1,]  <- vt[i-1,]/sum(vt[i-1,]);
		if(i%%10000==0) cat("vt  ",i,"\n")

	}

elas.s <- 0
elas.s.mean <- 0

elas.s.tmp <- matrix(0,nBigMatrix,nBigMatrix)

for (year.t in n.runin:(n.est-n.runin)) {

		#standard calculations needed for the various formulae

	                     
	vt1.C.wt                    <- sum(vt[year.t+1,] * C.t[year.i[year.t],,] %*% wt[year.t,])
	        
	vt1.C.wt.mean        <- sum(vt[year.t+1,] * C.t.mean[year.i[year.t],,] %*% wt[year.t,])
	 	
	K           <- K.year.i[year.i[year.t],,]
	
	vt1.K.wt    <- sum(vt[year.t+1,] * (K %*% wt[year.t,]))
	
	vt1.wt                       <- outer(vt[year.t+1,],wt[year.t,],FUN="*")
	vt1.C.wt.tmp           <- vt1.wt * C.t[year.i[year.t],,]  
	        



		#calculation of the standard elasticities
	
		elas.s            <-elas.s + (vt1.C.wt) / vt1.K.wt;
		elas.s.tmp    <-elas.s.tmp + (vt1.C.wt.tmp) / vt1.K.wt;
		elas.s.mean <-elas.s.mean + (vt1.C.wt.mean) / vt1.K.wt;
		
}

 elas.s             <- elas.s/(n.est-2*n.runin+1)
 elas.s.tmp     <- elas.s.tmp/(n.est-2*n.runin+1)
 elas.s.mean  <- elas.s.mean/(n.est-2*n.runin+1)
 
 
 return(list(meshpts=year.K$meshpts, h=h, elas.s=elas.s, elas.s.mean=elas.s.mean, mean.kernel=mean.kernel, Ls=Ls, elas.s.tmp=elas.s.tmp))

}






