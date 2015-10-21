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
    return(dnorm(x1,mean=mu,sd=sd.g))
} 
## Integral of Growth 
intg_x1x <- function(x1,x) {
       mu <-  x + 325*(x^1.3)/((144+(x^2.3)/42)^2); 
       return(pnorm(x1,mean=mu,sd=sd.g))
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
# Functions to make iteration matrices 
##############################################################################

# midpoint rule 
mk_K <- function(m, L, U) {
	h <- (U - L)/m; meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, p_x))
	F <- h * (outer(meshpts, meshpts, f_x))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F, h = h))
}


# Gauss-Legendre 
mk_GLK <- function(m, L, U) {
	out<-gaussQuadInt(L,U,m) 
    meshpts=out$nodes; wts=out$weights; 
	P <- (outer(meshpts, meshpts, p_x));
	F <- (outer(meshpts, meshpts, f_x));
	for(i in 1:m) {P[,i]=P[,i]*wts[i]; F[,i]=F[,i]*wts[i]} 
	K <- P + F
	return(list(K = K, meshpts = meshpts, wts=wts, P = P, F = F))
}



# B2B with cumulative kernel (z') and Gauss-Legendre (z) 
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
		return(list(K = K, meshpts = meshpts, P = P, F = F,h=h))
}

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
  return(list(K = K, meshpts = meshpts, P = P, F = F,h=h))
}

sd.g=1; 
system.time(mk_K_int(300,L=0,U=140,order=3))
system.time(mk_K_int_noloop(500,L=0,U=140,order=3))


#########################################################################
# Compute dominant eigenvalue by bin-to-bin, using cumulative kernel in 
# in destination bin and GL(3) in source bin 
##########################################################################
sdvals=c(0.1,1); results=as.list(1:2); 
for(k in 1:2) {
sd.g=sdvals[k]; 
msizes=cbind(100+(0:8)*100, 100+(0:8)*50);  
lam <- lamB2B <- lamGL <- Ft <- FtGL <- FtB2B <- numeric(6); 

for(j in 1:9) {
  msize <- bsize <- msizes[j,k]; 
  
  mpr <- mk_K(msize,0,140); eig <- domEig(mpr$K) 
  lam[j] <- eig$lambda; 
  h <- mpr$h; 
  wmpr <- eig$w; wmpr <- 1000*wmpr/sum(h*wmpr); 
  Ft[j]<- h*sum(wmpr*b_x(mpr$meshpts))
 
  gl <- mk_GLK(msize,0,140); eig <- domEig(gl$K); 
  lamGL[j] <- eig$lambda;
  weights=gl$wts; 
  wgl <- eig$w; wgl <- 1000*wgl/sum(weights*wgl); 
  FtGL[j]<- sum(weights*wgl*b_x(gl$meshpts)); 
  
  b2b <- mk_K_int_noloop(bsize,0,140,3); eig <- domEig(b2b$K); 
  lamB2B[j] <- eig$lambda; 
  h <- b2b$h; 
  wb <- eig$w; wb <- 1000*wb/sum(h*wb); 
  FtB2B[j]<- h*sum(wb*b_x(b2b$meshpts)); 
  cat(j,Ft[j],FtGL[j],FtB2B[j],"\n");
  
}
results[[k]]=cbind(msizes,Ft,FtGL,FtB2B);  
}

results1=results[[1]]; results2=results[[2]]; 

graphics.off(); dev.new(height=4.25,width=8);  
set_graph_pars("panel2"); par(yaxs="i"); 
matplot(results1[,1],results1[,3:5],type="o",lty=1, col="black",lwd=1,cex=1.2,
pch=c(16,1,2),xlab="Matrix size (# meshpoints)",ylab="Annual total fecundity",ylim=c(0,120)); 
title(main="Growth StdDev = 0.1",cex.main=0.8); 
legend("topright",c("Midpoint","Gauss-Legendre","Bin-to-bin"),lty=1,lwd=1,pch=c(16,1,2),
bty="n",cex=0.9); 
add_panel_label("a"); 

matplot(results2[,2],results2[,3:5],type="o",lty=c(1,2,3), col="black", lwd=1,cex=1.2,
pch=c(16,1,2),xlab="Matrix size (# meshpoints)",ylab=" ",ylim=c(80,110)); 
title(main="Growth StdDev = 1",cex.main=0.8); 
add_panel_label("b"); 

dev.copy2eps(file="../../c6/figures/ThreeWayF.eps")
dev.copy2pdf(file="../../c6/figures/ThreeWayF.pdf")

#plot(mpr$meshpts,log(wmpr),type="l");
#points(b2b$meshpts,log(wb),col="red",type="l"); 
#points(gl$meshpts,log(wgl),col="blue",type="l"); 
# dev.copy2eps(file="../../c6/figures/TreeB2B.eps")
