rm(list=ls(all=TRUE))
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c6",sep="")); 

require(Matrix); require(statmod); 
source("../Utilities/matrixImage.R"); 
source("../Utilities/Standard Graphical Pars.R"); 

source("domEig.R"); source("GaussQuadSub.R"); 

############################################################### 
#        Functions that make up the kernels  
###############################################################

## Growth 
g_x <- function(xp,x) { 
    mu <-  x + 325*(x^1.3)/((144+(x^2.3)/42)^2); 
    return(dnorm(xp,mean=mu,sd=0.1))
} 

mu_x <-  function(x) x + 325*(x^1.26)/((144+(x^2.26)/42)^2)  

## Survival 
s_x <- function(x) 0.98 

## Fecundity 
b_x <- function(x) {u = -3+ 0.15*x; return(ifelse(x>10, 0.83*exp(u)/(1+exp(u)), 0)) }    

# Offspring size     
c_x <- function(xp) dnorm(xp,mean=1.1,sd=0.35) 

# Kernel functions 
p_x <- function(xp,x) s_x(x)*g_x(xp,x)
f_x <- function(xp,x) b_x(x)*s_x(x)*c_x(xp) 
k_x <- function(xp,x) p_x(xp,x)+f_x(xp,x)  

## make the midpoint rule iteration matrix 
mk_K <- function(m, L, U) {
	h <- (U - L)/m; meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, p_x))
	F <- h * (outer(meshpts, meshpts, f_x))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}

## make the Gauss-Legendre iteration matrix 
require(statmod); 
mk_GLK <- function(m, L, U) {
	out<-gaussQuadInt(L,U,m) 
    meshpts=out$nodes; wts=out$weights; 
	P <- (outer(meshpts, meshpts, p_x));
	F <- (outer(meshpts, meshpts, f_x));
    for(i in 1:m) {P[,i]=P[,i]*wts[i]; F[,i]=F[,i]*wts[i]} 
	K <- P + F
	return(list(K = K, meshpts = meshpts, wts=wts, P = P, F = F))
}

##########################################################################
# Compute dominant eigenvalue by midpoint rule and GL, for different matrix sizes
##########################################################################
msizes <- seq(100,800,length=10); lam <- lamGL <- numeric(10);   
for(j in 1:10) { 
    out <- mk_K(msizes[j],0,140); 
    eig=domEig(out$K); lam[j]=eig$lambda; 
    out <- mk_GLK(msizes[j],0,140); 
    eig=domEig(out$K); lamGL[j]=eig$lambda; 
    cat(j,lam[j],lamGL[j],"\n"); 
}

# compute dominant eigenvector for 100 grid points: (1,0,0,...0)  
out=mk_K(msizes[1],0,140); eigen(out$K)$vectors[,1]; 

#########################################################################
# Compute dominant eigenvalue by bin-to-bin, using GL(7) in destination
# bin and GL(3) in source bin to compute the bin averages 
##########################################################################
lamB2B <- numeric(10); msizes2=50+(0:9)*25; 
for(j in 5:5) {
  msize=msizes2[j]; KB2B=matrix(0,msize,msize); 
  mpr=mk_K(msize,0,140); 
  meshpts=mpr$meshpts; h=meshpts[2]-meshpts[1]; 
  out1=gaussQuadInt(-h/2,h/2,7); 
  out2=gaussQuadInt(-h/2,h/2,3); 
  for(i in 1:msize) {
  for(k in 1:msize) {
    KB2B[i,k]=quad2D(FUN=k_x,wts1=out1$weights,wts2=out2$weights, 
        nodes1=out1$nodes+meshpts[i],nodes2=out2$nodes+meshpts[k])
  } }
 KB2B=KB2B/h; lamB2B[j]=domEig(KB2B)$lambda;    
 cat(j,lamB2B[j],"\n");
 if(j==5) {x.b2b=meshpts; w.b2b=domEig(KB2B,tol=1e-12)$w; w.b2b=w.b2b/mean(w.b2b) };
} 

##########################################################################
# Plot and compare eigenvalues and size distributions 
##########################################################################   
out <- mk_K(1000,0,140); 
x.mpr=out$meshpts;
w.mpr = domEig(out$K,tol=1e-12)$w; w.mpr=w.mpr/mean(w.mpr); 

dev.new(height=4,width=8); set_graph_pars("panel2"); 
matplot(cbind(msizes,msizes2), cbind(lam,lamB2B),xlab="Size of iteration matrix", 
ylab="Dominant eigenvalue", type="o",lty=1,pch=c(16,1),lwd=1,
ylim=c(0,max(lam)),col="black"); 
abline(h=lamB2B[10],lty=2); 
legend("topright",c("Bin-to-Bin","Midpoint"),pch=c(1,16),lty=1,bty="n")
add_panel_label("a")

plot(x.mpr,sqrt(w.mpr),type="l",xlab="Tree dbh (cm)",ylab="Sqrt(relative frequency)",lty=1,xlim=c(0,140)); 
points(x.b2b,sqrt(w.b2b),type="p",col="grey50",lty=1); 
legend("topright",c("Bin-to-Bin (150)","Midpoint (1000)"),col=c("grey50","black"),
lty=1,pch=c(1,NA), bty="n");   
add_panel_label("b")
dev.copy2eps(file="../../c6/figures/TreeB2B.eps")
