## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IBM to illustrate the approximation to demographic stochasticity 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 

source("../c2/Monocarp Demog Funs.R");
source("../utilities/Standard Graphical Pars.R")
source("../c7/Carlina/Carlina Demog Funs DI.R");
source("Monocarp IBM One Step.R");

constantSeeds=FALSE; 

##########################################################################
# Simulate IBM to get an initial population
##########################################################################
init.pop.size <- 250; n.yrs <-80

### Simulate IBM with mean parameters to get initial pop. in stable state ####
init.pop.size <- 250; 
z <- rnorm(init.pop.size, mean = m.par.true["rcsz.int"], sd = m.par.true["rcsz.sd"])
for(j in 1:100) {
  z <- oneIBMstep(z,m.par.true,constantSeeds); zmin=min(z); 
  if(length(z)>25000) z<- sample(z,25000,replace=FALSE)  
  cat(j,length(z),zmin,"\n"); 
} 

initial.pop <- sample(z,1000,replace=FALSE); 

##########################################################################
# Do 500 replicates of taking one time step forward from initial.pop 
###########################################################################
nt1 <- matrix(NA,5000,100); binlimits=seq(2,4.25,length=101);
for(j in 1:5000) {
q <- oneIBMstep(initial.pop,m.par.true,constantSeeds) ;
q <- hist(q,breaks=c(-15,binlimits,15),plot=FALSE,warn.unused=FALSE)$counts;
q <- q[2:(length(q)-1)];
nt1[j,] <- q;
if(j%%500==0) cat(j,"\n") 
}

##########################################################################
# Monte Carlo P-value of K-S test against Poisson distribution 
# x is vector of data values, B is the number of Monte-Carlo samples 
##########################################################################
KStest <- function(x,B) {
    lam <- mean(x) 
    efun <- function(y) ppois(y,lam); 
    D0 <- ks.test(x, efun)$statistic;
    Dvals <- numeric(B);
    for(j in 1:B) {
        xj <- rpois(length(x),mean(x))
        efun <- function(y) ppois(y,mean(xj)); 
        Dvals[j] <- ks.test(xj,efun)$statistic
    }    
    return(mean(Dvals>=D0)); 
}    

## Compute P-values for one-step-ahead simulations
pvals<-numeric(100); 
for(j in 1:100) {
    pvals[j]<-KStest(nt1[,j],B=250); 
    cat(j,"\n");
}    

#######################################################################
# Plot results 
#######################################################################

graphics.off(); 
# dev.new(height=4,width=8)
set_graph_pars("panel4"); par(cex.lab=1.4,cex.axis=1.25)

# Plot histogram of bin 50 across the replicate simulations 
hist(nt1[,50],xlab="Number in Interval 50", ylab="Frequency",main="")
add_panel_label("a"); 

# Plot variance versus mean for the 100 bins 
h=binlimits[2]-binlimits[1]
binMeans <- apply(nt1,2,mean);
binVars<- apply(nt1,2,var)

plot(binMeans,binVars, xlab="Bin mean",ylab="Bin Variance");
abline(0,1); add_panel_label("b")

# Plot P-values (KS test against Poisson) versus bin mean
plot(binMeans,pvals,type="p",xlab="Bin mean",ylab="KS test P-value");
abline(h=0.05,lty=2); 
add_panel_label("c"); 

# Plot covariance versus distance for all pairs of bins 
out <- cor(nt1); diag(out) <- NA;
nbins <- length(binMeans);
dmat <- h*outer(1:nbins,1:nbins,FUN=function(x,y) abs(x-y));
diag(dmat) <- NA; 
plot(dmat,out,xlab="Distance between bins",ylab="Correlation between bins",pch=16,cex=0.35);
add_panel_label("d")

#dev.copy2eps(file="../../c10/figures/MonocarpDemogStoch.eps")



