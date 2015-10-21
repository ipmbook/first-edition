require(mgcv); 

########### Set working directory and set graphing pars 
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c2",sep="")); 
source("../utilities/Standard Graphical Pars.R");

graphics.off(); set_graph_pars("panel1")
par(cex.axis=1.35,cex.lab=1.5);

N <- 150; # size of the 'data' set 

############# Generate and plot one typical data set, and the true mean 
z <- runif(N,1,50); z <- sort(z); 

z1bar <- 40*z/(15+z); z1 <- z1bar+rnorm(N,mean=0,sd=5);

plot(z,z1,xlab="Size at time t",ylab="Size at time t+1",cex=1.2); 
points(z,z1bar,col="black",type="l",lwd=3)

#####################################################################
# fit the data with a spline
# Note: m=3 specifies 3rd derivative penalty so the fit can have
# nonzero curvature at the endpoints of the data range.
# With the default (m=2), curvature -> 0 at the endpoints.
#####################################################################
X <- data.frame(z=z,z1=z1)
gamGrow <- gam(z1~s(z,m=3),data=X);

############# plot the spline fit to the mean
newdata <- data.frame(z=seq(1,50,length=250))
z1hat <- predict(gamGrow,newdata=newdata,type="response")
points(newdata$z,z1hat,type="l",col="red",lty=2,lwd=3);
sse <- sum(gamGrow$residuals^2); sdhat<- sqrt(sse/gamGrow$df.residual)	

################# Use the fit to define the growth kernel function.
Gz1_z <- function(z1,z) {
	Gdata=data.frame(z=z);
	z1bar <- predict(gamGrow,newdata=Gdata,type="response")	
	return(dnorm(z1,mean=z1bar,sd=sdhat))
}

# Look at the variation over many replicate fits. 
# Use different variable names to avoid overwriting
# the objects used in Gz1_z above 
nReps=1000; 
yhat <- matrix(0,250,nReps) 
sighat <- numeric(nReps)
newdata <- data.frame(x=seq(1,50,length=250))
for(j in 1:nReps) {
  x <- runif(N,1,50); x <- sort(x); 
  x <- seq(1,50,length=N);
  y0 <- 40*x/(15+x); y <- y0+rnorm(N,mean=0,sd=5);
  X <- data.frame(x=x,y=y)
  fitj <- gam(y~s(x,m=3),data=X);
  yhat[,j] <- predict(fitj,newdata=newdata,type="response")
  sse <- sum(fitj$residuals^2); sighat[j]<- sqrt(sse/fitj$df.residual)	
  cat(j,sighat[j],"\n")
}

# plot pointwise 5th and 95th percentiles of the plot
r90 <- function(x) quantile(x,probs=c(0.05,0.95))
yr <- apply(yhat,1,r90)
#matpoints(newdata$x,t(yr),type="l",lty=3,col="blue",lwd=4);
polygon(c(newdata$x,rev(newdata$x)),c(yr[1,],rev(yr[2,])),col="turquoise",border="NA")

points(z,z1); 
points(z,z1bar,col="black",type="l",lwd=3)
points(newdata$x,z1hat,type="l",col="red",lty=2,lwd=3);

#dev.copy2eps(file="../../c10/figures/gamExample.eps")


