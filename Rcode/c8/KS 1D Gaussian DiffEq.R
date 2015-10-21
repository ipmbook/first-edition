## IBM in 1 dimension with Gaussian kernel

rm(list=ls(all=TRUE))
require(SuppDists); 

## Define parameters 
sigma=4; 
R0 = 2

x=seq(-1200,1200,by=.5);
dij=outer(x,x,"-"); 
dij=0.5*dnorm(dij,mean=0,sd=sigma); 
sum(dij[,100]); 

n.yrs <-150
pop.size.t <- rcrit.t <- matrix(NA,n.yrs)

X=as.numeric(abs(x)<2);  
## iterate the model 
yr <- 1;
while(yr <= n.yrs) {
    ## Calculate local population density
    if(yr%%100==50) plot(x,X,yaxs="i",main=yr);
    xfar <- x[X>0.25];
    rcrit.t[yr]=diff(range(xfar))/2; 


    ## apply density-dependent mortality 
    survivalProb <- 1/(1+0.5*pmax(X-2,0)); 
    X <- R0*X*survivalProb; 
  
    ## Store population size
    pop.size <- sum(X) 
    pop.size.t[yr] <- pop.size
    
  
    ## everybody move 
    X=dij%*%X; 


    if(yr%%25==0) cat(paste(yr, pop.size.t[yr], "\n", sep=" "))
    yr <- yr+1
   
}

## Plot the growth of occupied area over time 
plot(1:n.yrs,rcrit.t,type="l",lty=c(1,2,3,4),col="black",xlab="Years",
ylab="Radius of area occupied"); 

# Estimate spread rate from the last 25 years of spread. 
t.end <- seq(n.yrs-100,n.yrs,by=1);
fit<- lm(rcrit.t[t.end]~t.end)
spread.rate <- fit$coef[2]

## moment generating function 
Gaussmgf <- function(s) {  exp(sigma^2*s^2/2) }

## wave speed c(s)
cs <- function(s) { (1/s)*log(R0*Gaussmgf(s)) }

plot(cs,0.05,2)

svals=seq(0.05,2,length=500); 
csvals=cs(svals);
jmin=which(csvals==min(csvals));
cat(csvals[jmin],"\n"); 
cat(spread.rate,"\n"); 
