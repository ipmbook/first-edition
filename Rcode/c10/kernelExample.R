require(mgcv); 

########### Set working directory and set graphing pars 
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 
source("../utilities/Standard Graphical Pars.R");

graphics.off();

set_graph_pars("panel1")
par(cex.axis=1.35,cex.lab=1.5);
N <- 150; # size of the 'data' set 

################# Generate a 'data' set 
z <- runif(N,1,50); z <- sort(z); 
z1bar <- 40*z/(15+z); z1 <- z1bar+4*rt(N,df=3);

################# fit the data with a spline
X <- data.frame(z=z,z1=z1)
gamGrow <- gam(z1~s(z,m=3),data=X);

################# extract the residuals and estimate growth variance
gamResids <- resid(gamGrow); 
sse <- sum(gamResids^2); 
sdhat <- sqrt(sse/gamGrow$df.residual)	

############## estimate bandwidth for kernel density
h <- bw.SJ(gamResids);

############# scale residuals so kernel variance = estimated growth variance
alpha <- sdhat/sqrt(sdhat^2+h^2);
hResid <- alpha*gamResids

############ define the density estimate
kfun <- function(z) {
	mean(dnorm(z,mean=hResid,sd=h))
} 
kfun <- Vectorize(kfun) 

############ Use the fit and density estimate to define the growth kernel
Gz1_z <- function(z1,z) {
	Gdata=data.frame(z=z);
	z1bar <- predict(gamGrow,newdata=Gdata,type="response")	
	return(kfun(z1-z1bar)) 
}

###################################################################
# Look at the variation over many replicate fits. 
# This uses different variable names to avoid overwriting
# the objects used in Gz1_z above 
###################################################################
nReps=500; px <- seq(-4*sdhat,4*sdhat,length=500)
yhat <- matrix(0,length(px),nReps) 
for(j in 1:nReps) {
  x <- runif(N,1,50); x <- sort(x); 
  x <- seq(1,50,length=N);
  y0 <- 40*x/(15+x); y <- y0+4*rt(N,df=5);
  X <- data.frame(x=x,y=y)
  fitj <- gam(y~s(x,m=3),data=X);
  Resids <- resid(fitj); 
  sse <- sum(Resids^2); 
  sdj <- sqrt(sse/fitj$df.residual)	
  hj <- bw.SJ(Resids);
  aj <- sdj/sqrt(sdj^2+hj^2);
  hjResid <- aj*Resids
  kj  <- function(z) mean(dnorm(z,mean=hjResid,sd=hj)) 
  kj <- Vectorize(kj) 
  yhat[,j]=kj(px)
  cat(j,sdj,"\n")
}

r90 <- function(x) quantile(x,probs=c(0.05,0.95))
yr <- apply(yhat,1,r90)


matplot(px,cbind(kfun(px),0.25*dt(px/4,df=3)),type="l", ylim=c(0,max(yr)),
col=c("red","black"),lty=c(2,1),lwd=2,xlab="Growth residual",ylab="Probability density");

polygon(c(px,rev(px)),c(yr[1,],rev(yr[2,])),col="turquoise",border="NA");

matpoints(px,cbind(kfun(px),0.25*dt(px/4,df=3)),type="l", ylim=c(0,max(yr)),
col=c("red","black"),lty=c(2,1),lwd=2,xlab="Growth residual",ylab="Probability density")

df=3; Gsd=4*sqrt(df/(df-2))
matpoints(px,dnorm(px,mean=0,sd=Gsd),type="l",lty=3,lwd=3)

dev.copy2eps(file="../../c10/figures/kernelExample.eps")
