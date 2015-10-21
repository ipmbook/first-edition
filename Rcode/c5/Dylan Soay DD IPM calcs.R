setwd("~/Repos/ipm_book/Rcode/c5")
source("../utilities/Standard Graphical Pars.R")
source("Ungulate DD Demog Funs.R")

# Find the equilibrium population density
# by seeing what density gives lambda=1 

dens.vals <- seq(0, 500, length=100)
## lambda = f(density)
lamd.vals <- sapply(dens.vals,
                    function(Nt) {
                        IPM.dd <- mk_K(Nt, 100, m.par.true, 1.50, 3.55)
                        Re(eigen(IPM.dd$K)$val[1])
                    })
plot(dens.vals, lamd.vals*dens.vals,
     xlab="N(t)", ylab="N(t+1)",
     xlim=range(dens.vals), ylim=range(dens.vals), type="l")
abline(0, 1, lty=2, col=grey(0.5))
## seems OK, 146-413 females observed in the data time series, and pop is growing

## Iterate density-dependent IPM to see what happens
## start with 100 individuals at stable state distribution 
L <- 1.60; U <- 3.7; m=100; h<-(U-L)/m; 
IPM.dd <- mk_K(Nt=100, m=100, m.par=m.par.true, L=L,U=U)
meshpts <- IPM.dd$meshpts; 

n0 <- Re(eigen(IPM.dd$K)$vectors[,1]); 
n0 = 100*n0/sum(h*n0); sum(h*n0); # should be 100, it is. 

## Plot initial population
plot(meshpts,n0,type="l",xlab="Size z", ylab="Initial Soay Density")

Ntvals <- numeric(25); Ntvals[1] <- 100;
ntmat <- matrix(0,m,25); ntmat[,1] <- n0; 
for(k in 2:25) {
	IPM.dd <- mk_K(Nt=Ntvals[k-1], m=m, m.par=m.par.true, L=L,U=U)
	ntmat[,k] <- IPM.dd$K%*%ntmat[,k-1]
	Ntvals[k] <- h*sum(ntmat[,k])
}
set_graph_pars("panel4")
plot(1:25,Ntvals,type="o")
matplot(meshpts,ntmat,type="l")

## Can we produce a bifurcation?
## Lets change parameters so lambda crashes fast
## when N is above the equilibrium 
m.par.pert <- m.par.true
#m.par.pert["surv.int"]<- 2
#m.par.pert["surv.nt"] <- 5/200;

## That didn't work. Try brute force:
## add a Ricker that depends on total population size. 

IPM.dd <- mk_K(Nt=100, m=100, m.par=m.par.pert, L=L,U=U)
meshpts <- IPM.dd$meshpts; 
n0 <- Re(eigen(IPM.dd$K)$vectors[,1]); 
n0 = 100*n0/sum(h*n0); sum(h*n0); # should be 100, it is. 

Ntvals <- numeric(25); Ntvals[1] <- 100;
ntmat <- matrix(0,m,25); ntmat[,1] <- n0; 
for(k in 2:25) { 
	Nt=Ntvals[k-1]
	IPM.dd <- mk_K(Nt=Nt, m=m, m.par=m.par.pert, L=L,U=U)
	ntmat[,k] <- (IPM.dd$K%*%ntmat[,k-1])*6*exp(-Nt/200)
	Ntvals[k] <- h*sum(ntmat[,k])
}
plot(1:25,Ntvals,type="o")
matplot(meshpts,ntmat,type="l")



if(FALSE){
lamd.vals.pert <- sapply(dens.vals,
                         function(Nt) {
                             IPM.dd <- mk_K(Nt, 100, m.par.pert, 1.50, 3.55)
                             Re(eigen(IPM.dd$K)$val[1])
                         })
plot(dens.vals, lamd.vals.pert*dens.vals,
     xlab="N(t)", ylab="N(t+1)",
     xlim=range(dens.vals), ylim=range(dens.vals), type="l")
abline(0, 1, lty=2, col=grey(0.5))
}
