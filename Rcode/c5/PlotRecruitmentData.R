# Plot and analyze recruitment data on Oenothera and Platte thistle

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c5",sep="")); 

source("../utilities/Standard Graphical Pars.R")
require(MASS); require(stats4);

graphics.off(); 
dev.new(height=4,width=9);
set_graph_pars("panel4")
par(mfrow=c(1,3),cex.lab=1.5,cex.axis=1.5,cex.main=1.3)

# Load the data 
Oeno <- read.csv("OenotheraRecruit.csv"); 
Platte <- read.csv("PlatteThistleFig4.csv")
# average two digitisations of the Platte plotted data 
Platte$seeds <- (Platte$seeds+Platte$seeds2)/2;
Platte$recruits <- round((Platte$recruits+Platte$recruits2)/2);

# Plot Oenothera recruitment data and regression line 
plot(Oeno$seeds,Oeno$recruits,xlab="Seeds per m^2, year t",
ylab="Recruits per m^2, year t+1",main="Oenothera");
fit <- lm(recruits~seeds,data=Oeno)
abline(fit)
add_panel_label("a")

# Plot Platte recruitment data 
plot(Platte$seeds,Platte$recruits,xlab="Seed production, year t",ylab="Recruits, year t+1",
main="Platte thistle    ");
add_panel_label("b")

## Try Poison regression
PlattePois <- glm(recruits~log(seeds),family="poisson",data=Platte);
# This turns out to be highly overdispersed:
# Residual deviance:  710.87  on  9  degrees of freedom

### So try Negative Binomial regression, thanks to MASS  
PlatteNB1 <- glm.nb(recruits~log(seeds),link="log",data=Platte)
PlatteNB2 <- glm.nb(recruits~log(seeds)-1,link="log",data=Platte)
#   Intercept is nonsignificant so use the second fit
# 	Parameter estimates are slope=0.67, size=2.14

# Plot the fitted mean 
px=seq(0,1.1*max(Platte$seeds),length=100); 
py=px^0.67;
points(px,py,type="l",lty=1,lwd=1)

## add percentiles
py1 <- qnbinom(0.05,mu=px^0.67,size=2.14)
py2 <- qnbinom(0.95,mu=px^0.67,size=2.14)
matpoints(px,cbind(py1,py2),type="l",col="black",lty=2);

### Diagnostic: see if a log-log plot actually looks linear 
plot(log(Platte$seeds),log(Platte$recruits),xlab="log Seed production, year t",
ylab="log Recruits, year t+1", main="Platte thistle     ");
add_panel_label("c")

## Fit linear regression. Intercept is non-significant,
## so plot the fit without intercept  
PlatteLN1 <- lm(log(recruits)~log(seeds),data=Platte)
PlatteLN2 <- lm(log(recruits)~log(seeds)-1,data=Platte)
abline(PlatteLN2); 

## add percentiles of fitted distribution 
shat <- sqrt(mean(PlatteLN2$residuals^2))
px=seq(4,9,length=100); 
py1=PlatteLN2$coef[1]*px + qnorm(0.05)*shat;
py2=PlatteLN2$coef[1]*px + qnorm(0.95)*shat;
matpoints(px,cbind(py1,py2),type="l",lty=2,col="black");

dev.copy2eps(file="../../c5/figures/PlotRecruitmentData.eps")

## It's nice that we have glm.nb in MASS, but when the model we need hasn't 
## been coded for us somebody else, there's always maximum likelihood. 
## Fit with mle, first two parameters, then one (no intercept on log-log scale).  

nllNegbin <- function(a,b,size) {
	d<-dnbinom(Platte$recruits,mu=a*Platte$seeds^b,size=size);
	return( -sum(log(d)));
}

fitnb1 <- mle(minuslogl=nllNegbin,start=list(a=1,b=0.67,size=100),
	method="Nelder-Mead",control=list(trace=1,maxit=5000));

nllNegbin2 <- function(b,size) {
	d<-dnbinom(Platte$recruits,mu=Platte$seeds^b,size=size);
	return( -sum(log(d)));
}

fitnb2 <- mle(minuslogl=nllNegbin2,start=list(b=0.67,size=100),
	method="Nelder-Mead",control=list(trace=1,maxit=5000));

AIC(fitnb1,fitnb2) # simpler model is has lower AIC, though not by much