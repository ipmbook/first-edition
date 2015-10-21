## Working directory must be set here, so the source()'s below run
setwd("~/Repos/ipm_book/Rcode/c5")
source("Platte Demog Funs.R");
source("../utilities/Standard Graphical Pars.R");

## Create the P iteration and IPM structures 
L <- (-1.5); U <- 4.5; m<-200; 
out <- mk_P(m,m.par,L,U);
P <- out$P; meshpts <- out$meshpts; h <- meshpts[2]-meshpts[1];

####################################################################
# Function to compute average number of weevils per
# flowering plant 
#####################################################################
weevilLoad <- function(nt,meshpts,m.par) {
	nFlower <- h*sum(nt*p_bz(meshpts));    # total flowering plants
	eps <- exp(m.par["weevil.int"] + 1.71*meshpts); 
	nWeevil <- h*sum(nt*p_bz(meshpts)*eps)  # total weevils on those plants
	return(nWeevil/nFlower)
}

############################################################################
# Trace effects of weevil infestation by varying the intercept parameter
# e_0 in the size-dependent weevil burden function. Values from Rose et al. 
# Table 2 are annual estimates for 11 years starting with 1990. Here we 
# use average values over 2 or 3 year time periods.  
############################################################################
int <- Ntot <- wbar <- numeric(4)
int[1] <- mean (c(-17.3,-17,-16.7))
int[2] <- mean (c(-3.74,-3.39,-2.81))
int[3] <- mean(c(-1.64,-0.96,-0.75)) 
int[4] <- mean(c(-1.24,-0.52))
int[5] <- 1.5; # beyond the range of the data 

graphics.off(); 
dev.new(height=4.5,width=9)
set_graph_pars("panel2")
nt<-matrix(1,m,200);

# Iterate to find stable population for 1st value of intercept parameter 
m.par["weevil.int"] <- int[1]
for(j in 2:200) {
	nt[,j] <- Iterate(nt[,j-1],meshpts,P,m.par)
}
matplot(meshpts,nt[,200],type="l",col="black",xlab="Log root crown diameter, z",
ylab="Size distribution n(z)",lty=1,lwd=2)
legend("topright",c("1990-92","1993-95","1996-98","1999-2000"),lty=1:4,lwd=2,bty="n")
Ntot[1] <- h*sum(nt[,200])
wbar[1]<-weevilLoad(nt[,200],meshpts,m.par);

# Stable population for other values of intercept parameter 
for(k in 2:5) {
  m.par["weevil.int"] <- int[k]
  for(j in 2:200) {
	  nt[,j] <- Iterate(nt[,j-1],meshpts,P,m.par)
  }	
  if(k<5) matpoints(meshpts,nt[,200],type="l",col="black",lty=k,lwd=2)
  Ntot[k] <- h*sum(nt[,200])	
  wbar[k]<-weevilLoad(nt[,200],meshpts,m.par);
  }
add_panel_label("a")

########################################################################
## Population size as a function of mean weevil eggs per flowering plant  
########################################################################
plot(wbar[1:4],Ntot[1:4],type="p",pch=1,cex=1.2,xlab="Mean weevil eggs per flowering plant",
ylab="Equilibrium population size",ylim=c(0,max(Ntot)),xlim=c(0,max(wbar)) );
add_panel_label("b")

## Extend the range of weevil infestation 
int=wbar=Ntot = seq(min(int),max(int),length=50);
for(k in 1:50) {
  m.par["weevil.int"] <- int[k]
  for(j in 2:200) {
	  nt[,j] <- Iterate(nt[,j-1],meshpts,P,m.par)
  }	
  Ntot[k] <- h*sum(nt[,200])	
  wbar[k]<-weevilLoad(nt[,200],meshpts,m.par);
  }
points(wbar,Ntot,type="l",lwd=2)

#####################################################################
## Change the model so that all plants have the mean weevil load
## given their size, and recompute effects of weevil infestation  
#####################################################################
b_z <- function(z, m.par) {
    eps <- exp(m.par["weevil.int"] + 1.71*z); 
    N <- exp(-0.55 + 2.02*z - 0.02*eps)
    return(N)
}
P <- mk_P(m,m.par,L,U)$P

## Recompute effects of weevil infestation 
int=wbar=Ntot = seq(min(int),max(int),length=50);
for(k in 1:50) {
  m.par["weevil.int"] <- int[k]
  for(j in 2:200) {
	  nt[,j] <- Iterate(nt[,j-1],meshpts,P,m.par)
  }	
  Ntot[k] <- h*sum(nt[,200])	
  wbar[k]<-weevilLoad(nt[,200],meshpts,m.par);
  }
points(wbar,Ntot,type="l",lty=2,lwd=2)

legend("topright",c("Negative binomial","Constant"),lty=c(1,2),lwd=2,bty="n")

dev.copy2eps(file="../../c5/figures/PlatteCalculations.eps")


#################################################################################
## Compute R0(0) as a function of intercept parameter 
#################################################################################

# go back to the original model with negative binomial weevil load 
source("Platte Demog Funs.R");

# Because P does not depend on weevil load intercept parameter, we can
# compute the fundamental operator N= (I-P)^{-1} once and for all 
P=mk_P(m,m.par,L,U)$P; Nmat <- solve(diag(m)-P);  

# We can also compute the fraction of time in each state, 
# for a cohort of newborns 
Nc0 <- Nmat%*%c_0z1(meshpts,m.par) 

## Birth rate depends on weevil load intercept parameter so we have to 
## compute the total offspring number for each parameter value, 
## for Seeds=0 and hence p_e=1.  

e0vals <- seq(-20,12,length=100); 
R0vals <- numeric(length(e0vals));

m.par.temp <- m.par; 
for(j in 1:length(e0vals)) {
	m.par.temp["weevil.int"] <- e0vals[j];
	R0vals[j]=h*sum(Nc0*B_z(meshpts,m.par.temp)); 
}	

# mean weevil burden on a typical-size flowering plant 
# as a function of intercept parameter 
meanEggs <- function(e0) exp(e0+1.71*2.5); 

# Find the e0 value for extinction assuming maximum p_e
# is either 1 or 0.25
Efun<-splinefun(R0vals,e0vals);
meanEggs(Efun(1)); # max p_e = 1
meanEggs(Efun(4)) # max p_e = 0.25

graphics.off(); dev.new(height=5,width=7)
set_graph_pars("panel1"); par(pty="m"); 
par(cex.axis=1.4,cex.lab=1.4,yaxs="i",mgp=c(2.25,1,0)); 
plot(meanEggs(e0vals),R0vals,xlab="Expected weevils/plant at z=2.5",ylab=expression(R[0](0)),log="x",
	type="l",lwd=2,ylim=c(0,R0vals[1])+0.25); 

abline(h=1,lty=2,lwd=2); 
arrows(x0=meanEggs(-17.3),y0=10,x1=meanEggs(-0.52),y1=10,code=3,col="blue",lwd=3); 	
dev.copy2eps(file="../../c5/figures/PlatteR0.eps")	
dev.copy2pdf(file="../../c5/figures/PlatteR0.pdf")	

