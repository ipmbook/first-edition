# Fit a dispersal kernel using the distance from feeder
# to final resting point for the 158 jay-cached acorns, f
# from Pons & Pausas 2007. Data provided by J. Pausas, 
# July 31, 2013, and used with permission. 

graphics.off();
require(car); require(stats4);
source("../utilities/Standard Graphical Pars.R")


# The data 
J <- read.csv("Pausas2Ellner.csv"); 

# Does feeder of origin have an effect? Yes, but not a large one
fit=lm(log(dist)~factor(feeder),data=J); 
summary(fit);

# Does it look like an exponential? No.  
set_graph_pars("panel4")

truehist(J$dist,nbins=15,xlab="Distance",ylab="Probability density", main="",col="grey75");  
add_panel_label("a")
xvals=seq(0,max(J$dist));
points(xvals,dexp(xvals,rate=1/mean(J$dist)),type="l",lty=2,lwd=2);

qqPlot(J$dist,distribution="exp",xlab="Exponential quantiles",ylab="Distance",main="");
add_panel_label("b")

# Does it look like a lognormal? Yes! 
truehist(log(J$dist),nbins=15,xlab="Log distance",ylab="Probability density", main="",col="grey75");  
mulog <- mean(log(J$dist)); sdlog<-sd(log(J$dist)); 
xvals <- seq(mulog-3*sdlog,mulog+3*sdlog,length=100);
points(xvals,dnorm(xvals,mulog,sdlog),type="l",lty=2,lwd=2);
add_panel_label("c")

qqPlot(log(J$dist),xlab="Gaussian quantiles",ylab="Log distance",main="");
add_panel_label("d")
dev.copy2eps(file="../../c8/figures/PonsPausas.eps")

# Do a formal test for non-normality of log(distance)
shapiro.test(log(J$dist)); 

# Compare dist~lognormal or exponential or gamma using AIC
# Let R do the work: fit by maximum likelihood and use AIC() 

# m,s are mean and std dev of log(distance)
nllLognormal <- function(m,s) {
	d<-dlnorm(J$dist,meanlog=m,sdlog=s,log=TRUE)
	return( -sum(d) )
}	

# m is the mean of distance
nllExp <- function(m) {
	d<-dexp(J$dist,rate=1/m,log=TRUE);
	return(-sum(d));
}

# a,s are shape and scale parameters 
nllGamma <- function(a,s) {
	d<-dgamma(J$dist,shape=a,scale=s);
	return( -sum(log(d)));
}

fit1 <- mle(minuslogl=nllLognormal,start=list(m=mulog,s=sdlog),
	method="Nelder-Mead",control=list(trace=4,maxit=5000));
	
meanD <- mean(J$dist);	
fit2 <- mle(minuslogl=nllExp,start=list(m=meanD), method="Brent",
	lower=meanD/2,upper=2*meanD);

# estimate gamma parameters by method of moments;
varD=var(J$dist);
s0 <- varD/meanD; a0 <- meanD/s0; 

fit3 <- mle(minuslogl=nllGamma,start=list(a=a0,s=s0),
	method="Nelder-Mead",control=list(trace=4,maxit=5000));

AIC(fit1,fit2,fit3)
