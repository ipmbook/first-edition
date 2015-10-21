## IBM in 1 dimension with C. nutans kernel (Inverse Gaussian) 

#rm(list=ls(all=TRUE))
require(SuppDists); 

## Define parameters 
sigma=8; 
R0 = 2.25


## moment generating function 
Gaussmgf <- function(s) {
  mgf <- exp(sigma^2*s^2/2); 
  return(mgf)
}  

cs <- function(s) {
    (1/s)*log(R0*Gaussmgf(s))
}

plot(cs,0.05,0.25)

svals=seq(0.05,0.25,length=1000); 
csvals=cs(svals);
jmin=which(csvals==min(csvals));
cat(csvals[jmin]);