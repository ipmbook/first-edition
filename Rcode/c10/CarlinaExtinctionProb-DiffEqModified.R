## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use Carlina IBM/IPM to illustrate Vindenes et al. (2010) approximations 
## 
## The code here, incorporating contributions from Robin Snyder, does one more
## diffusion approximation than appears in the book: the Lande (Oikos, 1998)
## transformed diffusion on a scale such that the noise is isotropic  
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos");
setwd(paste(root,"/first-edition/Rcode/c10",sep="")); 

source("../c7/Carlina/Carlina Demog Funs DI.R");
source("../utilities/Standard Graphical Pars.R")
source("Monocarp IBM One Step.R");
source("VindenesDiffusionApproxFunsModified.R"); 

load("LowNoiseCarlinaDiffusionPars.Rdata");

constantSeeds=FALSE; 

###############################################################################
#  Function to simulate the approximate dynamics of V(t)  
#  NOTE: we do this here without passing to the diffusion approximation in 
#        continuous time, as Vindenes et al. (2010) do 
############################################################################### 
sigmaD=sqrt(sigma2D); sigmaE=sqrt(sigma2E); 

diffEQ <- function(V0,nreps,tmax,mean.loglambda,sd.loglambda,sigmaD) {            
    Vt=matrix(NA,tmax,nreps);
    Vt[1,]=V0; 
    for(j in 2:tmax) {
       ext <- (Vt[j-1,]<0); 
       Vt[j-1,ext] <- 0; 
       lam.t = rlnorm(nreps,mean.loglambda,sd.loglambda)
       Vt[j,]= Vt[j-1,]*lam.t + rnorm(nreps,mean=0,sd=sigmaD*sqrt(Vt[j-1,])); 
     }    
    return(Vt)
}    

graphics.off(); set_graph_pars("panel4"); layout(mat=matrix(c(1,2,3,3),2,2,byrow=TRUE));
par(pty="m"); 

##############################################################################
# Simulate DiffEQ for V with Carlina parameters, plot results about extinction 
##############################################################################
nreps=10000; tmax=101; V0=25; 

### Some sample trajectories 
system.time(Vt <- diffEQ(V0,nreps,tmax,mean.loglambda,sd.loglambda,sigmaD)); 
matplot(1:100,log10(1+Vt[1:100,1:100]),type="l",col="black",xlab="Time t",ylab="Log10(1+ V(t))",lty=1);         
add_panel_label("a"); 

### fraction extinct as a function of time, for sample trajectories 
extinct=apply(Vt,1,function(x) mean(x==0)); 
plot(1:tmax,extinct,type="l",col="black",xlab="Time t",ylab="Fraction extinct");  
add_panel_label("b"); 

### fraction extinct after 50 time steps, as function of initial V 
L0=ceiling(10^seq(0,2,length=8)); E25 <- numeric(length(L0)); 
j=1; 
for(L in L0) {
    Vt <- diffEQ(L,nreps=10000,tmax=50,mean.loglambda,sd.loglambda,sigmaD); 
    E25[j] <- mean(Vt[25,]==0); 
    j <- j+1; cat(j,"\n"); 
}    
par(yaxs="i"); 
plot(L0,E25,type="o",lwd=1,xlab="Initial V", ylab="Fraction extinct at t=50",log="x",
    ylim=c(0,1),pch=16,cex=1.5,lty=2); 

#############################################################################
### Next, use the Vindenes et al diffusion approximation for log V, and again  
### find fraction extinct after 50 time steps as function of initial V
#############################################################################
j=1; 
for(L in L0) {
    Vt <- diffusion(y0=L,t.max=50,delta.t=.1,n.sim=10000,b=Inf,envar=sigma2E, demvar=sigma2D, lambda=lambda); 
    Vt <- t(Vt); 
    E25[j] <- mean(Vt[nrow(Vt),]==0); 
    j <- j+1; cat(j,"\n"); 
}    
points(L0,E25,type="o",pch=1,lty=2,cex=1.5); 

#############################################################################
### Next, use the Lande 1998 diffusion approximation (Oikos 83:
### 353- 358, eq. 7b, 8), and again find fraction extinct after 50
### time steps as function of initial V
#############################################################################
j=1;
ensd = sqrt(sigma2E)
bconst = 2*log(1 + sqrt(1 + sigma2D/sigma2E))
for(L in L0) {
    Vt <- diffusion.Lande(y0=L,t.max=50,delta.t=.1,n.sim=10000,b=Inf,envar=sigma2E, demvar=sigma2D, lambda=lambda, ensd=ensd, bconst=bconst); 
    Vt <- t(Vt); 
    E25[j] <- mean(Vt[nrow(Vt),]==0); 
    j <- j+1; cat(j,"\n"); 
}    
points(L0,E25,type="o",pch=3,lty=2,cex=1.5); 

##############################################################################
### Next, use log transformation of the diffusion approximation for V   
### using Ito's Lemma, and do the same extinction calculations
##############################################################################
j=1; 
for(L in L0) {
    Vt <- diffusion.Ito(y0=L,t.max=50,delta.t=.1,n.sim=10000,b=Inf,envar=sigma2E, demvar=sigma2D, lambda=lambda); 
    Vt <- t(Vt); 
    E25[j] <- mean(Vt[nrow(Vt),]==0); 
    j <- j+1; cat(j,"\n"); 
}    
points(L0,E25,type="o",pch=2,lty=2,cex=1.5); 

add_panel_label("c"); 

#####################################################################################
# Simulate the Carlina IPM with different starting populations to look at extinction 
#####################################################################################
IBMsim <- TRUE; 
if(IBMsim) {
  Nstart=ceiling(10^seq(0,2,length=8)); 
  pext <- numeric(length(Nstart));
  for(k in 1:length(Nstart)){
    zfinal <- numeric(500); 	
    for(j in 1:500) {
      zt <- sample(initial.pop,Nstart[k],replace=TRUE); v0 <- sum(vfun(zt));
      for(T in 1:50) {
	m.par.year=m.par.true + qnorm(runif(12,0.001,0.999))*m.par.sd.true
	m.par.year[c("grow.int","rcsz.int")] <- m.par.year[c("grow.int","rcsz.int")] + 
          matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 
        zt <- oneIBMstep(zt,m.par.year,constantSeeds=FALSE); 
        if(length(zt)>25000) break; 
      }  
      if(j%%50==0) cat(k,j,v0, length(zt),"\n"); 
      zfinal[j] <- length(zt)
    }
    pext[k] <- mean(zfinal<1)
  }

  points(Nstart,pext,type="o",lty=1,lwd=2,pch=17,cex=1.5)
} 

legend("topright",legend=c("IBM simulations","Difference Eqn Approx","Diffusion Approx 3","Diffusion Approx 2", "Lande diffusion approx"),
    lty=c(1,2,2,2,2),lwd=c(2,1,1,1,1),pch=c(17,16,1,2,3),bty="n",cex=1.05); 



