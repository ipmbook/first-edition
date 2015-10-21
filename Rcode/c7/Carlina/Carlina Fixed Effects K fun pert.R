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
#Let's do the probability of flowering function
########################################################################

#Select the parameters to use
params.to.use <- m.par.est

#Calculate meshpts and h for evaluating the perturbation kernels
year.K      <-  mk_K(nBigMatrix,params.to.use[,1],minsize,maxsize)
meshpts <-  year.K$meshpts
h              <-  year.K$h

#First calculate the mean function and perturbation kernels

p_bz.mean <- 0

for(i in 1:n.years){
	p_bz.mean <- p_bz.mean + p_bz(meshpts, params.to.use[,i])
}

p_bz.mean <- p_bz.mean/n.years

Ct_z1z <- function(z1,z,m.par){
		return( p_bz(z, m.par) * s_z(z, m.par) *
				 ( m.par["p.r"] * b_z(z, m.par) * c_0z1(z1, m.par) - G_z1z(z1, z, m.par)) )
	 }
	
C.pert <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.C <-h * (outer(meshpts, meshpts, Ct_z1z, m.par = params.to.use[,i]))
		C.pert[i,,] <- year.C
	}

Ct_z1z_mean <- function(z1,z,m.par){
		return( p_bz.mean * s_z(z, m.par) *
				  ( m.par["p.r"] * b_z(z, m.par) * c_0z1(z1, m.par) - G_z1z(z1, z, m.par)) )
	}

C.pert.mean <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.C <-h * (outer(meshpts, meshpts, Ct_z1z_mean, m.par = params.to.use[,i]))
		C.pert.mean[i,,] <- year.C
	}

pert.K <- stoc_pert_analysis(params.to.use, n.est, n.runin, C.pert, C.pert.mean)



meshpts <- pert.K$meshpts
elas.s <- apply(pert.K$elas.s,2,sum)
elas.s.mean <- apply(pert.K$elas.s.mean,2,sum)
elas.s.sd <- elas.s - elas.s.mean
sens.mean <- elas.s.mean * pert.K$Ls / p_bz.mean
set_graph_pars("panel4")
plot(meshpts,elas.s,type="l",xlab="Size (t), z",ylab=expression(e[S] ^p[b]))
add_panel_label("a")
plot(meshpts,elas.s.mean,type="l",xlab="Size (t), z",ylab=expression(e[S] ^{p[b]*","*mu}))
add_panel_label("b")
plot(meshpts,elas.s.sd,type="l",xlab="Size (t), z",ylab=expression(e[S] ^{p[b]*","*sigma}))
add_panel_label("c")
plot(meshpts,sens.mean,type="l",xlab="Size (t), z",ylab=expression(s[S] ^{p[b]*","*mu}))
add_panel_label("d")
dev.copy2eps(file="~/Repos/ipm_book/c7/figures/CarlinapbElasSens.eps")


########################################################################
#Let's do the survival function
########################################################################

#First calculate the mean function and perturbation kernels

s.zmean <- 0

for(i in 1:n.years){
	s.zmean <- s.zmean + s_z(meshpts, params.to.use[,i])
}

s.zmean <- s.zmean/n.years

Ct_z1z <- function(z1,z,m.par){
		return( s_z(z, m.par) * (1- p_bz(z, m.par) )*G_z1z(z1, z, m.par) +
		              s_z(z, m.par) * p_bz(z, m.par) * m.par["p.r"] * b_z(z, m.par) * c_0z1(z1, m.par) ) 
	}
	
C.pert <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.C <-h * (outer(meshpts, meshpts, Ct_z1z, m.par = params.to.use[,i]))
		C.pert[i,,] <- year.C
	}

Ct_z1z_mean <- function(z1,z,m.par){
		return( s.zmean * (1- p_bz(z, m.par) )*G_z1z(z1, z, m.par) +
		              s.zmean * p_bz(z, m.par) * m.par["p.r"] * b_z(z, m.par) * c_0z1(z1, m.par) ) 
	}


C.pert.mean <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.C <-h * (outer(meshpts, meshpts, Ct_z1z_mean, m.par = params.to.use[,i]))
		C.pert.mean[i,,] <- year.C
	}

pert.K <- stoc_pert_analysis(params.to.use, n.est, n.runin, C.pert, C.pert.mean)

#max(pert.K$elas.tmp-apply(pert.K$elas.s,2,sum))

meshpts <- pert.K$meshpts
elas.s <- apply(pert.K$elas.s,2,sum)
elas.s.mean <- apply(pert.K$elas.s.mean,2,sum)
elas.s.sd <- elas.s - elas.s.mean
sens.mean <- elas.s.mean * pert.K$Ls / s.zmean

#Check sum elasticities is 1

cat(sum(pert.K$elas.s)," should be 1","\n")

set_graph_pars("panel4")
plot(meshpts,elas.s,type="l",xlab="Size (t), z",ylab=expression(e[S] ^s(z)))
add_panel_label("a")
plot(meshpts,elas.s.mean,type="l",xlab="Size (t), z",ylab=expression(e[S] ^{s(z)*","*mu}))
add_panel_label("b")
plot(meshpts,elas.s.sd,type="l",xlab="Size (t), z",ylab=expression(e[S] ^{s(z)*","*sigma}))
add_panel_label("c")
plot(meshpts,sens.mean,type="l",xlab="Size (t), z",ylab=expression(s[S] ^{s(z)*","*mu}))
add_panel_label("d")

dev.copy2eps(file="~/Repos/ipm_book/c7/figures/CarlinasElasSens.eps")	
	
########################################################################
#Let's do the growth function
########################################################################

#First calculate the mean function and perturbation kernels

G.mean <- 0

for(i in 1:n.years){
	G.mean <- G.mean + outer(meshpts, meshpts, G_z1z, m.par = params.to.use[,i])
}

G.mean <- G.mean/n.years

Ct_z1z <- function(z1,z,m.par){
		return( s_z(z, m.par) * (1- p_bz(z, m.par) )*G_z1z(z1, z, m.par)  ) 
	}
	
C.pert <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.C <-h * (outer(meshpts, meshpts, Ct_z1z, m.par = params.to.use[,i]))
		C.pert[i,,] <- year.C
	}

Ct_z1z_mean <- function(z1,z,m.par){
		return( s_z(z, m.par)  * (1- p_bz(z, m.par) ) * G.mean) 
	}


C.pert.mean <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.C <-h * (outer(meshpts, meshpts, Ct_z1z_mean, m.par = params.to.use[,i]))
		C.pert.mean[i,,] <- year.C
	}
	
pert.K <- stoc_pert_analysis(params.to.use, n.est, n.runin, C.pert, C.pert.mean)

meshpts <- pert.K$meshpts
elas.s.sd <- pert.K$elas.s - pert.K$elas.s.mean
sens.mean <- pert.K$elas.s.mean * pert.K$Ls / G.mean

## set up the plots
ikeep <- ikeep <- which(meshpts>1.5 & meshpts<5) # use to extract a region to plot
set_graph_pars("panel4")
## plot the growth sensitivity and elasticity surfaces
image(meshpts[ikeep], meshpts[ikeep], t(pert.K$elas.s[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Size (t), z", ylab="Size (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(pert.K$elas.s[ikeep,ikeep]), 
        add=TRUE)
add_panel_label("a")
## plot the offspring size kernel sensitivity and elasticity surfaces
image(meshpts[ikeep], meshpts[ikeep], t(pert.K$elas.s.mean[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Size (t), z", ylab="Size (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(pert.K$elas.s.mean[ikeep,ikeep]), add=TRUE)
add_panel_label("b")
image(meshpts[ikeep], meshpts[ikeep], t(elas.s.sd [ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Size (t), z", ylab="Size (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(elas.s.sd [ikeep,ikeep]), add=TRUE)
add_panel_label("c")
image(meshpts[ikeep], meshpts[ikeep], t(sens.mean[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Size (t), z", ylab="Size (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(sens.mean[ikeep,ikeep]), add=TRUE)
add_panel_label("d")

dev.copy2eps(file="~/Repos/ipm_book/c7/figures/CarlinaGElasSens.eps")	





########################################################
#Old code sums before averaging and explicit loop...
stoc.pert.analysis.old<-function(params,n.est,n.runin,C.t,C.t.mean){
	
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


elas.s <- rep(0,nBigMatrix)
elas.s.mean <- rep(0,nBigMatrix)

for (year.t in n.runin:(n.est-n.runin)) {

		#standard calculations needed for the various formulae

	
	vt1.C.wt                    <- sapply(1:nBigMatrix,function(z0) sum(vt[year.t+1,] *    
	                                           (C.t[year.i[year.t],,z0]  * wt[year.t,z0])))
	        
	 vt1.C.wt.mean        <- sapply(1:nBigMatrix,function(z0) sum(vt[year.t+1,] *    
	                                           (C.t.mean[year.i[year.t],,z0]  * wt[year.t,z0])))
	
	K           <- K.year.i[year.i[year.t],,]
	
	vt1.K.wt    <- sum(vt[year.t+1,] * (K %*% wt[year.t,]))
	
# # 	vt1.above <- rep(NA,nBigMatrix)
		
# for(z0 in 1:nBigMatrix){
#          pr.pb.c0.G           <- params["p.r",year.i[year.t]] * b_z(meshpts[z0],params[,year.i[year.t]]) *       
#                                                        c_0z1(meshpts,params[,year.i[year.t]]) - 
#                                                       G_z1z(meshpts,meshpts[z0],params[,year.i[year.t]])
#        pb.s.wt                  <-  p_bz(meshpts[z0],params[,year.i[year.t]]) * 
#                                                       s_z(meshpts[z0],params[,year.i[year.t]]) 
#         vt1.above[z0]     <- sum(vt[year.t+1,] * pb.s.wt * (pr.pb.c0.G * wt[year.t,z0])) *h
# }
	



		#calculation of the standard elasticities
	
		elas.s            <-elas.s + (vt1.C.wt) / vt1.K.wt;
		elas.s.mean <-elas.s.mean + (vt1.C.wt.mean) / vt1.K.wt;
		
}

 elas.s            <- elas.s/(n.est-2*n.runin+1)
 elas.s.mean <- elas.s.mean/(n.est-2*n.runin+1)
 
 
 return(list(meshpts=year.K$meshpts, h=h, elas.s=elas.s, elas.s.mean=elas.s.mean, mean.kernel=mean.kernel, Ls=Ls))

}




