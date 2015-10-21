## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use Carlina IBM/IPM to illustrate Vindenes et al. (2010) approximations 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 

source("../c7/Carlina/Carlina Demog Funs DI.R");
source("Monocarp IBM One Step.R");

constantSeeds=FALSE; 

############## Change parameters so that lambda is near 1, noise is "small" 
m.par.save <- m.par.true; 
m.par.sd.save <- m.par.sd.true; 
VarCovar.grow.rec.save <- VarCovar.grow.rec; 
chol.VarCovar.grow.rec.save <- chol.VarCovar.grow.rec;  

m.par.true["p.r"] <- 0.0004; 
m.par.sd.true <- m.par.sd.save/4; 
VarCovar.grow.rec <- VarCovar.grow.rec.save/16; 
chol.VarCovar.grow.rec <-  chol.VarCovar.grow.rec.save/4 

############################################################################
#  Get v and w of mean environment kernel; make reproductive value function 
############################################################################
# make kernel for mean environment parameters, in m.par.true; 
U <- 6.5; L <- 0; m=200;  
IPM <- mk_K(m, m.par.true, L, U) 
K <- IPM$K; meshpts<-IPM$meshpts; h<- IPM$h; 

eK <-eigen(K); 
lambda <- Re(eK$values[1]); w <- eK$vectors[,1]; w=abs(w)/sum(h*abs(w)); 
v <- abs(eigen(t(K))$vectors[,1]); v <- v/(h*sum(v*w)); 

vfun <- approxfun(meshpts,v); 

###############################################################################
#  Estimate demographic stochasticity variance parameter sigma2D using the IBM
############################################################################### 

#######
# First approach: evaluate the integral s_d^2(z) w(z) dz by midpoint rule 
# with s_d^2(z) evaluated by repeated simulations of the IBM
#######

nreps=2500; vVar <- numeric(m); 
for(j in 1:m) {
	V1<-numeric(nreps);
	for(k in 1:nreps) {
	  z <- oneIBMstep(meshpts[j],m.par.true,constantSeeds);
	  V1[k] <- sum(vfun(z))
	}
	vVar[j]=var(V1)
	cat(j,log(vVar[j]),"\n")
}
sigma2D <- sum(h*w*vVar)

##################################################################################
#  Second approach: simulate variance of V(t+1) for population in steady state
##################################################################################

### Simulate IBM with mean parameters to get initial pop. in stable state ####
init.pop.size <- 250; 
z <- rnorm(init.pop.size, mean = m.par.true["rcsz.int"], sd = m.par.true["rcsz.sd"])
for(j in 1:500) {
  z <- oneIBMstep(z,m.par.true,constantSeeds); zmin=min(z); zmax=max(z); 
  if(length(z)>25000) z<- sample(z,25000,replace=FALSE)  
  if(j%%10==0) cat(j,length(z),zmin,zmax,"\n"); 
} 

# We'll do this 5 times for different draws from steady state distribution 
sigma2D.2 <- numeric(5);
for(k in 1:5){
initial.pop <- sample(z,5000,replace=FALSE); 
V.init <- sum(vfun(initial.pop));  
###### Do 2500 replicates of taking one time step forward #####
V1 <-  matrix(NA,2500); 
for(j in 1:2500) {
   q <- oneIBMstep(initial.pop,m.par.true,constantSeeds=FALSE) ;
   V1[j] <- sum(vfun(q));
   if(j%%100==0) cat(j,"\n") 
}

### Resulting estimate of demographic stochasticity variance parameter  #######
sigma2D.2[k] <- var(V1)/V.init; 
}
sigma2D.2;


###############################################################################
#  Estimate the environmental stochasticity variance parameter, sigma2E
#  This is equivalent to Appendix A in Vindenes et al 2011 but uses the
#  second line of the expression for lambda(Z) on p.3, instead of the third. 
############################################################################### 
V.init <- h*sum(v*w); lam.t <- matrix(NA,2500); 
for(j in 1:2500){ 
    m.par.year=m.par.true + qnorm(runif(12,0.001,0.999))*m.par.sd.true
    m.par.year[c("grow.int","rcsz.int")] <- m.par.year[c("grow.int","rcsz.int")] + 
        matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 
    K.year <- mk_K(m, m.par.year, L, U)$K; 
    lam.t[j] <- h*sum(v*(K.year%*%w))/V.init; 
    if(j%%100==0) cat(j,"\n") 
}

sigma2E <- var(lam.t);  # for diffusion approximation 
mean.loglambda <- mean(log(lam.t)); # for diffEQ approximation
sd.loglambda  <- sd(log(lam.t));    # for diffEQ approximation

save.image(file="LowNoiseCarlinaDiffusionPars.Rdata"); 
