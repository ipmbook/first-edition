## IBM in 1 dimension with C. nutans kernel (Inverse Gaussian) 

#rm(list=ls(all=TRUE))
require(SuppDists); 

## Define parameters 
r.nu      =   2.8 
r.lambda =    3.3
R0 = 1.5

## Movement kernel: abs(distance)~InverseGaussian
## Function to generate random movements 
displacement<-function(n) {
  r <- rinvGauss(n,nu=r.nu,lambda=r.lambda);
  u <- runif(n,-1,1); 
  return(r*u/abs(u)) 
}

## moment generating function 
Waldmgf <- function(s) {
  nu <- r.nu
  lambda <- r.lambda;
  t1 <- (lambda/nu); 
  t2 <- 2*(nu^2)*s/lambda; 
  mgf <- exp(t1*(1-sqrt(1-t2)));
  return(mgf)
}  

BiWaldmgf <- function(s) {
  0.5*(Waldmgf(s)+Waldmgf(-s))
}  

s.max <- function() {
  nu <- r.nu
  lambda <- r.lambda;
  return(lambda/(2*nu*nu)); 
}

cs <- function(s) {
    (1/s)*log(R0*BiWaldmgf(s))
}

plot(cs,0.95*s.max(),s.max())

svals=seq(0.05,s.max(),length=1000); 
csvals=cs(svals);
jmin=which(csvals==min(csvals));
cat(csvals[jmin]);

