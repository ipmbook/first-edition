## IBM in 1 dimension with Gaussian kernel

rm(list=ls(all=TRUE))
require(SuppDists); 

# so figures will get saved to where they should 
setwd("/Users/stephenellner/Repos/ipm_book/c8/figures")

## Define population parameters 
sigma=4; 
R0 = 2

## Movement kernel: abs(distance)~Gaussian
displacement<-function(n) {rnorm(n,mean=0,sd=sigma) }
plot(density(displacement(5000)));

## Define simulation parameters; 
init.pop.size <- 250;
n.yrs <-500;
n.reps <- 10;

## vectors etc. to hold results for plotting 
pop.size.t <- rcrit.t <- matrix(NA,n.yrs,n.reps)
spread.rate <- matrix(NA,20,n.reps);
xvals=list(100); 

################################################# 
# Iterate the model: loop over K and replicates 
#################################################
for(K in 3:20) {
xcrit = K-2*(R0-1);	
for(rep in 1:n.reps) {
yr <- 1; X <- displacement(init.pop.size); 
while(yr <= n.yrs) {
    ## Calculate local population density
    out <- density(X,bw=sigma);
    denfun <- approxfun(out$x,length(X)*out$y,rule=2);
    ymax=max(length(X)*out$y);
    if(yr%%100==50) plot(denfun,1.05*min(X),1.05*max(X),ylim=c(0,ymax),yaxs="i",main=yr);
    density <- denfun(X)
    xfar <- X[density>0.1];
    rcrit.t[yr,rep]=diff(range(xfar))/2; 

    if((K==10)&(yr<=100)&(rep==1))  xvals[[yr]] <- X;
    
    ## apply density-dependent mortality 
    survivalProb <- 1/(1+0.5*pmax(density-xcrit,0)); 
    Surv <- rbinom(n=length(X),prob=survivalProb,size=1);
    X <- X[Surv==1];
  
    ## Store population size just before breeding 
    pop.size <- length(X) 
    pop.size.t[yr,rep] <- pop.size
    
    ## Generate new recruits, and put them in parent locations
    newRecruits <- rpois(pop.size,lambda=R0);
    X=rep(X,newRecruits); 
    
    ## everybody move 
    X=X+displacement(length(X)); 

    if(yr%%100==0) cat(paste(rep,yr, pop.size.t[yr,rep], "\n", sep=" "))
    yr <- yr+1
   
} # end loop over years 

# Estimate spread rate from the last 100 years of spread
t.end <- seq(n.yrs-100,n.yrs,by=1);
fit<- lm(rcrit.t[t.end,rep]~t.end)
spread.rate[K,rep] <- fit$coef[2]
} # end loop over replicates for one value of K 

## Plot the growth of occupied area over time 
matplot(1:n.yrs,rcrit.t,type="l",lty=c(1,2,3,4),col="black",xlab="Years",
ylab="Radius of area occupied"); 

# Estimate spread rate from the last 100 years of spread for each replicate. 
#for(rep in 1:n.reps) {
#    fit<- lm(rcrit.t[t.end,rep]~t.end)
#    spread.rate[K,rep] <- fit$coef[2]
#   }

cat(K,spread.rate[K,],"\n"); 
} # End loop over different values of K 
############### End iterations of the model  ######################

###################################################################
# Plot the results 
###################################################################
spread.rate=spread.rate[3:20,];

## Compute IPM wave speed
Gaussmgf <- function(s) {  exp(sigma^2*s^2/2) }
cs <- function(s) { (1/s)*log(R0*Gaussmgf(s)) }

plot(cs,0.05,2)
svals=seq(0.05,2,length=1000); 
csvals=cs(svals);
jmin=which(csvals==min(csvals));
cstar=csvals[jmin]; cat(cstar);

dev.new(width=9,height=4)
par(cex.lab=1.3,cex.axis=1.3,bty="l",mgp=c(2.5,1,0),mar=c(4,4,2,1),mfrow=c(1,2))

x=xvals[[100]]; out=density(x,bw=20); plot(out$x,length(x)*out$y,type="l",
		xlab="Location x",ylab="Population density",ylim=c(0,16.5),yaxs="i");
x=xvals[[50]]; out=density(x,bw=20); points(out$x,length(x)*out$y,type="l",lty=3,lwd=2); 
x=xvals[[25]]; out=density(x,bw=20); points(out$x,length(x)*out$y,type="l",lty=2,lwd=2);
mtext("A)",side=3,adj=0,cex=1.65)

px=(25*cstar + 2 *sigma)*c(-1,1); py=rep(15.4,2)
points(px,py,type="l",lty=2,lwd=2)

px=(50*cstar + 2*sigma)*c(-1,1); py=rep(15.6,2)
points(px,py,type="l",lty=3,lwd=2)

px=(100*cstar + 2*sigma)*c(-1,1); py=rep(15.8,2)
points(px,py,type="l",lty=1,lwd=1)

meanSpreadRate <- apply(spread.rate,1,mean);
plot(3:20,meanSpreadRate,pch=16,type="p",ylim=c(3.8,cstar),xlab="Population limit K",
     ylab="Mean simulated wave speed");
abline(h=cstar,lty=3,lwd=2)     
mtext("B)",side=3,adj=0,cex=1.65)
dev.copy2eps(file="GaussianIBMWaveSpeeds.eps");
dev.copy2pdf(file="GaussianIBMWaveSpeeds.pdf")




