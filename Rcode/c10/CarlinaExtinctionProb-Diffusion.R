## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use Carlina IBM/IPM to illustrate Vindenes et al. (2010) approximations 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 

source("../c7/Carlina/Carlina Demog Funs DI.R");
source("../utilities/Standard Graphical Pars.R");
source("Monocarp IBM One Step.R");
source("VindenesDiffusionApproxFuns.R"); 

load("LowNoiseCarlinaDiffusionPars.Rdata"); 

constantSeeds <- FALSE; 

#########################################################################################
# Run diffusion approximation with Carlina parameters, plot results about extinction 
#########################################################################################
nreps <- 10000; t.max <- 100; V0 <- 5; 

Vt <- diffusion(y0=V0,t.max=t.max,delta.t=.1,n.sim=100,b=Inf,envar=sigma2E, demvar=sigma2D, lambda=lambda)
Vt <- t(Vt); # now each column is a simulation 

graphics.off(); set_graph_pars("panel4"); 
matplot(1:101,log10(1+Vt),type="l",col="black",xlab="Time t",ylab="Log10(1+ V)",lty=1);         
add_panel_label("a"); 

extinct=apply(Vt,1,function(x) mean(x==0)); 
plot(0:t.max,extinct,type="l",col="black",xlab="Time t",ylab="Fraction extinct");  
add_panel_label("b"); 

L0=seq(0,3,length=10); E100 <- numeric(length(L0)); 
j=1; 
for(L in L0) {
    Vt <- diffusion(y0=10^L,t.max=25,delta.t=.1,n.sim=5000,b=Inf,envar=sigma2E, demvar=sigma2D, lambda=lambda); 
    Vt <- t(Vt); 
    E100[j] <- mean(Vt[nrow(Vt),]==0); 
    j <- j+1; cat(j,"\n"); 
}    
par(yaxs="i"); 
plot(10^L0,E100,type="o",xlab="Initial V", ylab="Fraction extinct at t=25",log="x",ylim=c(0,1.05)); 
add_panel_label("c"); 

########################################################################################
# Simulate the Carlina IBM with different starting populations to look at extinction 
########################################################################################
Nstart=10^L0; pext <- numeric(length(Nstart));
for(k in 1:length(Nstart)){
zfinal <- numeric(500); 	
for(j in 1:500) {
  zt <- sample(initial.pop,Nstart[k],replace=TRUE); v0 <- sum(vfun(zt));
  for(T in 2:26) {
	m.par.year=m.par.true + qnorm(runif(12,0.001,0.999))*m.par.sd.true
	m.par.year[c("grow.int","rcsz.int")] <- m.par.year[c("grow.int","rcsz.int")]+matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 
    zt <- oneIBMstep(zt,m.par.year,constantSeeds=FALSE); 
    if(length(zt)>25000) break; 
  }  
  cat(k,j,v0, length(zt),"\n"); 
  zfinal[j] <- length(zt)
}
  pext[k] <- mean(zfinal<1)
}
plot(Nstart,pext,type="o",xlab="Initial V", ylab="Fraction extinct at t=25",log="x",ylim=c(0,1.05))
add_panel_label("d"); 

#dev.copy2eps(file="../../c10/figures/CarlinaDiffusion.eps")

