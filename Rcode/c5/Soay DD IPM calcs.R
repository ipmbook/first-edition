setwd("~/Repos/ipm_book/Rcode/c5")
source("../utilities/Standard Graphical Pars.R")
source("Soay DD Demog Funs.R")
##########################################################
# Plot density-dependent vital rates 
##########################################################
set_graph_pars("panel2"); 
par(mfrow=c(2,1),cex.lab=1.25)
par(pty="m",mar=c(4,4,1,4))

zMum <- 2.9;  # typical size of adult female 
Nvals <- seq(146,413,length=100); # range of population sizes in data 

svals <- s_z(zMum,Nvals,m.par.true)
prvals <- pr_z(Nvals,m.par.true);
z0vals <- m.par.true["rcsz.int"] + m.par.true["rcsz.Nt"] * Nvals + m.par.true["rcsz.z"] *zMum
z0vals <- exp(z0vals)*exp(-0.5*m.par.true["rcsz.sd"]^2);

matplot(Nvals,cbind(svals,prvals),type="l",lty=c(1,2),col="black",lwd=2,xlim=c(146,430),
     xlab="", ylab="Probability");
legend("bottomleft",legend=c("Adult survival","Recruitment probability","Offspring mass"),lty=c(1,2,3),
    col=c("black","black","blue"),bty="n",lwd=c(2,2,3))
     
par(new=TRUE)
plot(Nvals, z0vals,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(146,430),
     ylim=c(9,13),lwd=3,lty=3)
axis(side=4,at=c(8:14))
mtext("Offspring mean mass(kg)",side=4,line=-1.5,cex=1.1)     

add_panel_label("a")

######################################################### 
# Set the size range to preclude eviction. 
# Computing is fast for a simple model like this, so 
# we set broad limits that are more than sufficient
#########################################################

# Smallest recruits occur when Mum is small & population is large
# How small is Mum likely to be? 
pb_z(1,m.par.true);# <0.01, so very few parents of size z<1

# How small can new recruits be? 
# Use recruit size distribution for Mum's size z=1, and female 
# population 25% higher than any in the data 
rc <- function(z1) c_z1z(z1,1,520,m.par.true) 
integrate(rc,-Inf,0.5)  # <0.005, so L=0.5 is safe 

## Largest adults: compute asymptotic mean size from mean of the 
## growth distribution, and add 3 standard deviations 
z.max <- m.par.true["grow.int"]/(1-m.par.true["grow.z"])+3*m.par.true["grow.sd"]
# Slightly over 3.4; so extend by 20% by using U=3.65  

# Is U=3.65 big enough? See if anyone is evicted, using
# survivorship function at N=0 to maximize the result 
z1Dist <- function(z1) s_z(z1,0,m.par.true)*g_z1z(z1,3.65,m.par.true) 
integrate(z1Dist,3.65,Inf) # value=0.005, so not big enough

# Try again with U=3.8 
z1Dist <- function(z1) s_z(z1,0,m.par.true)*g_z1z(z1,3.8,m.par.true) 
integrate(z1Dist,3.8,Inf) # value=0.0003, fine. 

L <- 0.5; U <- 3.8; 
m <- 165; # gives h=0.02, so all indivs within 1.5% of interval midpoint 

#################################################################
# Find the equilibrium population density
# by seeing what density gives lambda=1 
#################################################################
### Plot lambda(N) as a function of N
### Where lambda(N)=1 is an equilibrium density 
lamFun <- function(Nt) {
    IPM.dd <- mk_K(Nt, m, m.par.true, L=L, U=U)
    return(Re(eigen(IPM.dd$K)$val[1]))
}

densVals <- seq(146, 416, length=100)
lamVals <- sapply(densVals,lamFun) 
plot(densVals, lamVals, xlab="Number of Females N", 
ylab=expression(paste("Population growth rate ",lambda)),type="l")
abline(h=1,lty=2)

# find the equilibrium 
Nbar <- uniroot(function(N) lamFun(N)-1, lower=200,upper=400)$root
points(Nbar,1,type="p",pch=16,cex=1.5)
## Model predicts equilibrium of 323 females 
## Seems OK, 146-413 females observed in the data, and pop is growing

add_panel_label("b")
dev.copy2eps(file="../../c5/figures/SoayDDFuns.eps")

#########################################################################
# Local stability analysis of equilibrium using Jacobian kernel
#########################################################################

## Iterate density-dependent IPM to find equilibrium,
## start with 50 individuals at stable distribution for N=50 
IPM.dd <- mk_K(Nt=50, m=m, m.par=m.par.true, L=L,U=U)
meshpts <- IPM.dd$meshpts; h<-(U-L)/m;

n0 <- Re(eigen(IPM.dd$K)$vectors[,1]); 
n0 = 50*n0/sum(h*n0); 
Ntvals <- numeric(21); Ntvals[1] <- sum(h*n0);
ntmat <- matrix(0,m,21); ntmat[,1] <- n0; 
for(k in 2:21) {
	IPM.dd <- mk_K(Nt=Ntvals[k-1], m=m, m.par=m.par.true, L=L,U=U)
	ntmat[,k] <- IPM.dd$K%*%ntmat[,k-1]
	Ntvals[k] <- h*sum(ntmat[,k])
}
### Plot size distributions after 5,10,15,20 years 
set_graph_pars("panel4")
matplot(meshpts,ntmat[,c(1,6,11,16,21)],type="l",xlab="Size z", ylab="Size distribution n(z,t)",
col=grey(1-1:6/6),lty=1)
legend("topleft",legend=as.character(c(1,6,11,16,21)-1),lty=1,col=grey(1-1:6/6),bty="n") 
add_panel_label("a")

### iterate a few more times to converge onto the equilibrium 
nbar <- ntmat[,21]
for(k in 1:25) {
	IPM.dd <- mk_K(Nt=sum(h*nbar), m=m, m.par=m.par.true, L=L,U=U)
	nbar <- IPM.dd$K%*%nbar
}

### Compute Jacobian matrix at equilibrium and find its dominant eigenvalue
J <- mk_J(nbar=nbar, m=m, m.par=m.par.true, L=L, U=U, eps=0.01);
lam_J <- eigen(J)$values[1]; cat(lam_J,"\n")

### Compare rate of convergence to nbar with the predicted rate, 
### i.e. each iteration decreases the distance by factor of abs(lam_J)
dvals <- numeric(21); 
for(j in 1:21) dvals[j] <- h*sum(abs(ntmat[,j]-nbar)) # observed rate

dhat <- abs(lam_J)^c(0:20); dhat <- dhat*dvals[21]/dhat[21]; # prediction 

dhat[1:7] <- NA; 
plot(0:20,log(dvals),type="p",xlab="Year t", ylab="ln(Distance from equilibrium)");
points(0:20,-0.35+log(dhat),type="l",lty=2,lwd=2)
add_panel_label("b")

nreps<-10; N0<-seq(50,450,length=nreps)
Ntvals<-matrix(NA,21,nreps)
for(j in 1:nreps) {
	nt <- runif(m) # random initial size distribution 
    nt <- N0[j]*nt/sum(h*nt) # scale to total population N0[j]
	Ntvals[1,j] <- sum(h*nt)
	for(k in 2:21) {
	   IPM.dd <- mk_K(Nt=sum(h*nt), m=m, m.par=m.par.true, L=L,U=U)
	   nt <- IPM.dd$K%*%nt  # iterate 
	   Ntvals[k,j] <- sum(h*nt) # store total population 
	}
}
matplot(0:20,Ntvals,type="o",pch=1,col="black",lty=1,xlab="Year t",ylab="Total number of females")
add_panel_label("c")

nreps<-10; N0<-seq(50,450,length=nreps)
Ntvals<-matrix(NA,21,nreps)
for(j in 1:nreps) {
	IPM.dd <- mk_K(Nt=N0[j], m=m, m.par=m.par.true, L=L,U=U)
	nt <- Re(eigen(IPM.dd$K)$vectors[,1]);
	nt <- N0[j]*nt/sum(h*nt); 
	Ntvals[1,j] <- sum(h*nt)
	for(k in 2:21) {
	   IPM.dd <- mk_K(Nt=sum(h*nt), m=m, m.par=m.par.true, L=L,U=U)
	   nt <- IPM.dd$K%*%nt
	   Ntvals[k,j] <- sum(h*nt)
	}
}
matplot(0:20,Ntvals,type="o",pch=1,col="black",lty=1,xlab="Year t",ylab="Total number of females")
add_panel_label("d")
dev.copy2eps(file="../../c5/figures/SoayDDStabilityAnalysis.eps")





