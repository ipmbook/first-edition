## IBM in 1 dimension with C. nutans kernel (Inverse Gaussian) 

rm(list=ls(all=TRUE))
require(SuppDists); 

## Define parameters 
r.nu      =   4 
r.lambda =    3
sinvGauss(r.nu,r.lambda)$Mean;
sinvGauss(r.nu,r.lambda)$SD;
qinvGauss(0.95,r.nu,r.lambda)
sigma.d=sinvGauss(r.nu,r.lambda)$SD;

R0 = 2

## Movement kernel: abs(distance)~InverseGaussian
displacement<-function(n) {
  r <- rinvGauss(n,nu=r.nu,lambda=r.lambda);
  u <- runif(n,-1,1); 
  return(r*u/abs(u)) 
}
plot(density(displacement(5000)),xlim=c(-10,10),yaxs="i");

init.pop.size <- 250
n.yrs <-500
n.reps <- 3
pop.size.t <- rcrit.t <- matrix(NA,n.yrs,n.reps)

for(rep in 1:n.reps) {

X=displacement(init.pop.size); 
## iterate the model 
yr <- 1;
while(yr <= n.yrs) {
    ## Calculate local population density
    out <- density(X,bw=2*sigma.d);
    denfun <- approxfun(out$x,length(X)*out$y,rule=2); 
    if(yr%%25==0) {plot(out$x,length(X)*out$y,xlab="x",type="l",
                        ylab="Population density", yaxs="i",main=yr);}
    density <- denfun(X)
    xfar <- X[density>0.1];
    rcrit.t[yr,rep]=diff(range(xfar))/2; 


    ## apply density-dependent mortality 
    survivalProb <- 1/(1+pmax(density-5,0)); 
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


    if(yr%%100==0) cat(paste(yr, pop.size.t[yr,rep], "\n", sep=" "))
    yr <- yr+1
   
}

}

## Plot the growth of occupied area over time 
matplot(1:n.yrs,rcrit.t,type="l",lty=c(1,2,3,4),col="black",xlab="Years",
ylab="Radius of area occupied"); 

# Estimate spread rate from the last 25 years of spread. 
t.end <- seq(n.yrs-100,n.yrs,by=1);
spread.rate <- numeric(n.reps); 
for(rep in 1:n.reps) {
    fit<- lm(rcrit.t[t.end,rep]~t.end)
    spread.rate[rep] <- fit$coef[2]
}
cat(spread.rate,"\n"); 


## moment generating function 
Waldmgf <- function(s) {
  nu <- r.nu; lambda <- r.lambda;
  t1 <- (lambda/nu); 
  t2 <- 2*(nu^2)*s/lambda; 
  mgf <- exp(t1*(1-sqrt(1-t2)));
  return(mgf)
}  

BiWaldmgf <- function(s) {
  0.5*(Waldmgf(s)+Waldmgf(-s))
}  

s.max <- r.lambda/(2*r.nu*r.nu)

cs <- function(s) {
  (1/s)*log(R0*BiWaldmgf(s))
}

plot(cs,0.95*s.max,s.max)

svals=seq(0.0001,s.max,length=1000); 
csvals=cs(svals);
jmin=which(csvals==min(csvals));
cat(csvals[jmin]);

