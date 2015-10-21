# Compare wave speed of Gaussian and Laplace

graphics.off(); 
setwd("~/Repos/ipm_book/Rcode/c8"); 

# set sigma=1;
sigma=1; b=sqrt(1/2); 

# Values of R0
R0=seq(1,100,length=250); 
r=seq(0,log(500),length=250); R0=exp(r);  


#Gaussian wave speed
cstar.G = sqrt(2*r)*sigma
#plot(R0,cstar.G,type="l"); 

# Laplace wave speed
eps=0.0001;
svals=seq(eps,1/b-eps,length=5000); 
M=1/(1-(b*svals)^2); 

nx=length(r); cstar.L=numeric(nx); 
for(j in 1:nx) {
	csvals=(1/svals)*(r[j]+log(M));
	cstar.L[j]=min(csvals);
}
dev.new(width=6,height=4)
par(cex.lab=1.3,cex.axis=1.3,bty="l",mgp=c(2.5,1,0),mar=c(4,4,1,1),xlog=TRUE,yaxs="i")
matplot(R0,cbind(cstar.G,cstar.L),type="l",lwd=2,log="x",
   xlab="Net reproductive rate R0",ylab="Wave speed c*",col="black"); 	


legend("topleft",legend=c("Gaussian","Laplace (bi-exponential)"),col="black",lwd=2,lty=c(1,2),bty="n",cex=1.2); 

dev.copy2eps(file="~/Repos/ipm_book/c8/figures/GaussLaplaceWaves.eps"); 
