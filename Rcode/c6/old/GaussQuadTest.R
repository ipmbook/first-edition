## Comparison of midpoint rule, Gauss-Legendre, and 
## subinterval-Gauss-Legendre quadrature
 
rm(list=ls(all=TRUE))
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c6",sep="")); 

require(statmod); source("GaussQuadSub.R"); 
source("../utilities/Standard Graphical Pars.R"); 
set_graph_pars("panel4"); par(yaxs="i")

L=1; U=5;  order=9; # order of GL on each subinterval 
fx <-function(x) 10*x*exp(-C*(x-a)^2) + 15*exp(-D*(x-b)^2); 

## 250 random functions with broad peaks 
gqs=gq=mpr=trp=sim=matrix(0,30,250); mvals=numeric(30); 
xvals=seq(L,U,length=250); fvals=matrix(0,250,250); 
for(k in 1:250) { 
  a<-runif(1,1,2); b<-runif(1,3,4); C=runif(1,5,25); D=runif(1,2,5); 
  fvals[,k]=fx(xvals); 
  exact=integrate(fx,L,U,subdivisions=2000)$value;  
  # loop over number of mesh points 
  for(j in 1:30) {
    intervals=j; m=order*intervals; mvals[j]=m; 
     
    out <- gaussQuadSub(L,U,order=order,intervals=intervals); 
    gqs[j,k] <- sum(out$weights*fx(out$nodes))/exact-1;
  
    out <- gaussQuadInt(L,U,order=m); 
    gq[j,k] <- sum(out$weights*fx(out$nodes))/exact-1; 
    
    out <- midpointInt(L,U,m); 
    mpr[j,k] <- sum(out$weights*fx(out$nodes))/exact-1; 
} 
cat("k=",k,"\n"); 
}
mpr=apply(abs(mpr),1,mean); 
gq=apply(abs(gq),1,mean); 
gqs=apply(abs(gqs),1,mean); 

matplot(xvals,fvals[,1:4],type="l",xlab="x",ylab="f(x)",lty=1:4); 
add_panel_label("a"); 
matplot(mvals,log10(cbind(mpr,gqs,gq)),type="l",lty=1:4,col=1:4,pch=1:4,lwd=2,
xlab="Number of mesh points",ylab="log10 Relative Error",xlim=c(0,100));
add_panel_label("b"); 
legend("right",c("Midpt","GLsub(9)","GL"),col=1:4,lty=1:4,bty="n",cex=1.1); 

## repeat, with sharper peaked functions 
gqs=gq=mpr=trp=sim=matrix(0,40,250); mvals=numeric(40); 
xvals=seq(L,U,length=250); fvals=matrix(0,250,250); 
## 250 random functions with broad peaks 
gqs=gq=mpr=trp=sim=matrix(0,30,250); mvals=numeric(30); 
xvals=seq(L,U,length=250); fvals=matrix(0,250,250); 
for(k in 1:250) { 
  a<-runif(1,1,2); b<-runif(1,3,4); C=runif(1,100,500); D=runif(1,5,10); 
  fvals[,k]=fx(xvals); 
  exact=integrate(fx,L,U,subdivisions=2000)$value;  
  # loop over number of mesh points 
  for(j in 1:30) {
    intervals=j; m=order*intervals; mvals[j]=m; 
     
    out <- gaussQuadSub(L,U,order=order,intervals=intervals); 
    gqs[j,k] <- sum(out$weights*fx(out$nodes))/exact-1;
  
    out <- gaussQuadInt(L,U,order=m); 
    gq[j,k] <- sum(out$weights*fx(out$nodes))/exact-1; 
    
    out <- midpointInt(L,U,m); 
    mpr[j,k] <- sum(out$weights*fx(out$nodes))/exact-1; 
} 
cat("k=",k,"\n"); 
}

mpr=apply(abs(mpr),1,mean); 
gq=apply(abs(gq),1,mean); 
gqs=apply(abs(gqs),1,mean);

matplot(xvals,fvals[,1:4],type="l",xlab="x",ylab="f(x)",lty=1:4); 
add_panel_label("c"); 
matplot(mvals,log10(cbind(mpr,gqs,gq)),type="l",lty=1:4,col=1:4,pch=1:4,lwd=2,
xlab="Number of mesh points",ylab="log10 Relative Error",xlim=c(0,200));
add_panel_label("d"); 

dev.copy2eps(file="../../c6/figures/GaussQuadTest.eps")

	