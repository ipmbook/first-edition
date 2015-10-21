## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IBM to illustrate the approximation to demographic stochasticity 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 

source("../c2/Monocarp Demog Funs.R");
source("../utilities/Standard Graphical Pars.R")
source("Monocarp IBM One Step.R"); 

constantSeeds=TRUE; 

##########################################################################
# Simulate IBM to get an initial population
##########################################################################
init.pop.size <- 250; n.yrs <-100
source("../c2/Monocarp Simulate IBM.R") 
cat(pop.size.t,"\n")

## trim an initial transient off the simulation
sim.data <- sim.data[sim.data$yr > 10,]
sim.data$yr <- sim.data$yr-10

# Extract 1000 observations to use as the starting population
sample.index <- sample(1:nrow(sim.data), size=1000, replace=FALSE)
sim.data <- sim.data[sample.index,]
initial.pop <- sim.data$z; 

#############################################################################
# Do 5000 replicates of taking 10 time steps forward from initial.pop in IBM 
#############################################################################
Nt10 <- numeric(5000); 

for(j in 1:5000) {
  nt <- initial.pop; 
  for(k in 1:20) {
     nt <- oneIBMstep(nt,m.par.true,constantSeeds);
  }
  Nt10[j] <- length(nt);
  if(j%%50==0) cat(j,range(nt),"\n") 
}

#############################################################################
# Do 5000 replicates of taking 10 time steps forward from initial.pop in IPM
# with demographic stochasticity  
#############################################################################
 L <- (-3); U= 5; m=125;
 IPM <- mk_K(m,m.par.true,L,U); 
 K <- IPM$K; meshpts=IPM$meshpts; 
 
 h <- (U-L)/m; 
 breakpts <- c(-100,meshpts[1]-h,meshpts+h,100); 
 initial.n <- hist(initial.pop,breaks=breakpts,plot=FALSE)$counts
 initial.n <- initial.n[2:(m+1)]; 
 
 Pt10 <- numeric(5000); 
 for(j in 1:5000) {
  dnt <- nt <- initial.n; 
  for(k in 1:20) {
     nt <- rpois(m,lambda=K%*%nt); # with demog stochasticity
     dnt <- K%*%dnt;               # without demog stochasticity
  }
  Pt10[j] <- sum(nt);
  if(j%%50==0) cat(j,"\n") 
}

graphics.off(); dev.new(); 
set_graph_pars("panel1");  par(yaxs="i",cex.lab=1.35); 
bwN=bw.SJ(Nt10)+bw.SJ(Pt10); bwN=bwN/2; 

plot(density((Nt10),bw=bwN),xlim=c(0,11000),xlab="Total population, t=20",ylab="Frequency",main="",cex=1.2); 
points(density((Pt10),bw=bwN),col="blue",lty=2,type="l",lwd=2);  
abline(v=sum(dnt),lty=3,col="red",lwd=4); 

legend("topright",c("IBM simulations","Poisson IPM"),col=c("black","blue"),lty=c(1,2), lwd=c(1,2),bty="n",cex=1.3); 

# dev.copy2eps(file="../../c10/figures/IterateDemogStoch.eps"); 

 