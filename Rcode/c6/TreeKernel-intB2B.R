# R script to implement the P. chinensis IPM using bin-to-bin
# with the cumulative kernel for final size and GL
# for initial size 

rm(list=ls(all=TRUE))
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c6",sep="")); 

require(Matrix); require(statmod); 
source("domEig.R"); source("GaussQuadSub.R"); 
source("../Utilities/matrixImage.R"); 
source("../Utilities/Standard Graphical Pars.R"); 

############################################################### 
#        Functions that make up the kernels  
###############################################################

## Growth 
g_x1x <- function(x1,x) { 
    mu <-  x + 325*(x^1.3)/((144+(x^2.3)/42)^2); 
    return(dnorm(x1,mean=mu,sd=0.1))
} 
## Integral of Growth 
intg_x1x <- function(x1,x) {
       mu <-  x + 325*(x^1.3)/((144+(x^2.3)/42)^2); 
       return(pnorm(x1,mean=mu,sd=0.1))
} 

## Survival 
s_x <- function(x) 0.98 

## Fecundity 
b_x <- function(x) {u = -3+ 0.15*x; return(ifelse(x>10, 0.83*exp(u)/(1+exp(u)), 0)) }    

# Offspring size     
c_x1 <- function(x1) dnorm(x1,mean=1.1,sd=0.35) 
# Integral of offspring size 
intc_x1 <- function(x1) pnorm(x1,mean=1.1,sd=0.35) 

# Kernel functions 
p_x <- function(xp,x) s_x(x)*g_x1x(xp,x)
f_x <- function(xp,x) b_x(x)*s_x(x)*c_x1(xp) 
k_x <- function(xp,x) p_x(xp,x)+f_x(xp,x)  

# Cumulative kernel functions 
intp_x <- function(x1,x) s_x(x)*intg_x1x(x1,x)
intf_x <- function(x1,x) b_x(x)*s_x(x)*intc_x1(x1) 
intk_x <- function(x1,x) intp_x(x1,x)+intf_x(x1,x)  

############################################################################## 
# Function to make the midpoint rule iteration matrix 
##############################################################################
mk_K <- function(m, L, U) {
	h <- (U - L)/m; meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, p_x))
	F <- h * (outer(meshpts, meshpts, f_x))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}

mk_K_int <- function(m,L,U,order) {
		h <- (U - L)/m; meshpts <- L + ((1:m) - 1/2) * h
  		out=gaussQuadInt(-h/2,h/2,order); 
  		nodes=out$nodes; weights=out$weights; 
  		F <- P <- matrix(0,m,m);
		for(i in 1:m){
			for(j in 1:m){
			      fvals1=intf_x(meshpts[i]-h/2,meshpts[j]+out$nodes);	
			      fvals2=intf_x(meshpts[i]+h/2,meshpts[j]+out$nodes);	
				  F[i,j]=sum(weights*(fvals2-fvals1))
				  pvals1=intp_x(meshpts[i]-h/2,meshpts[j]+out$nodes);	
			      pvals2=intp_x(meshpts[i]+h/2,meshpts[j]+out$nodes);	
				  P[i,j]=sum(weights*(pvals2-pvals1))
			}
		}
		P<- P/h; F<-F/h; K=P+F; 	
		return(list(K = K, meshpts = meshpts, P = P, F = F))
}

########## No-loop version of mk_K_int: much faster

mk_K_int_noloop <- function(m, L, U, order) {
  h <- (U - L) / m 
  meshpts <- L + ((1:m) - 1/2) * h 
  out <- gaussQuadInt(-h/2, h/2, order) 
  quad <- expand.grid(x1=meshpts, x=meshpts, map=seq.int(order))
  quad <- transform(quad, nodes=out$nodes[map], weights=out$weights[map])
  F <- with(quad, {
    fvals1 <- intf_x(x1 - h/2, x + nodes) 
    fvals2 <- intf_x(x1 + h/2, x + nodes) 
    weights * (fvals2 - fvals1) 
  })
  P <- with(quad, {
    fvals1 <- intp_x(x1 - h/2, x + nodes) 
    fvals2 <- intp_x(x1 + h/2, x + nodes) 
    weights * (fvals2 - fvals1) 
  })
  dim(F) <- dim(P) <- c(m, m, order)
  F <- apply(F, c(1,2), sum) / h
  P <- apply(P, c(1,2), sum) / h
  K <- F + P
  return(list(K = K, meshpts = meshpts, P = P, F = F))
}

#########################################################################
# Compute dominant eigenvalue by bin-to-bin, using cumulative kernel in 
# in destination bin and GL(3) in source bin 
##########################################################################
lam <- lamB2B <- numeric(10); msizes2=50+(0:9)*25; 
for(j in 1:10) {
  msize=msizes2[j]; 
  mpr=mk_K(msize,0,140); lam[j]=domEig(mpr$K)$lambda; 
  KB2B <- mk_K_int_noloop(msize,L=0,U=140,order=3)$K
  lamB2B[j]=domEig(KB2B)$lambda;    
 cat(j,lam[j],lamB2B[j],"\n");
 if(j==5) {x.b2b=mpr$meshpts; w.b2b=domEig(KB2B,tol=1e-12)$w; w.b2b=w.b2b/mean(w.b2b) };
} 

##########################################################################
# Plot and compare eigenvalues and size distributions 
##########################################################################   
out <- mk_K(1000,0,140); 
x.mpr=out$meshpts;
w.mpr = domEig(out$K,tol=1e-12)$w; w.mpr=w.mpr/mean(w.mpr); 

dev.new(height=4,width=8); set_graph_pars("panel2"); 
matplot(msizes2, cbind(lam,lamB2B),xlab="Size of iteration matrix", 
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
# dev.copy2eps(file="../../c6/figures/TreeB2B.eps")
