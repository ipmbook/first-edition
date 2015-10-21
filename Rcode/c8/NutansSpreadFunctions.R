## IBM to generate the data for the simple recruitment limited IPM - there's no recruitment limitation 1st
## Size is on a log scale so don't worry about the -ive ones!

## The parameters and life cycle follow Kachi and Hirose J.Ecol. 1985, see Rees and Rose PRSB 2002, but for
## simplicity we have used a logistic function for survival rather than the linear function with an upper bound
## as in the previous papers

## MORE COMMENTS STILL TO BE ADDED

rm(list=ls(all=TRUE))
require(SuppDists); 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so formulae are easier to
## read. We'll use a 'model.term' scheme to name elements of the vector, except that
## eps.pi is the product of pi and epsilon; pi and epsilon are only used in the calculation of
## seedling production, which is Poisson with mean epsilon*pi*(expected number of flowering heads) 

m.par.true <- c(surv.int  =  -2.27,
                surv.z    =   0.57,
                flow.int  =  -2.107,
                flow.z    =   0.86,
                grow.int  =   2.751,
                grow.z    =   0.407,
                grow.sd.int   = 1.5, 
                grow.sd.z     = 0.1,
                rcsz.int  =  -0.771,
                rcsz.sd   =   1.719,
                heads.int  =  6.363,
                heads.z    =  0.0056,
                eps.pi    =   7.1,
                r.nu      =   3.9, 
                r.lambda =    1.6)

## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets

## Growth function, given you are size z now returns the pdf of size z1 next time
g_z1z <- function(z1, z, m.par)
{
    mean <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year
    sd <- m.par["grow.sd.int"]/(1+m.par["grow.sd.z"]*z)       # sd about mean
    p.den.grow <- dnorm(z1, mean = mean, sd = sd)             # pdf that you are size z1 given you were size z
    return(p.den.grow)
}

## Survival function, logistic regression
s_z <- function(z, m.par)
{
    linear.p <- m.par["surv.int"] + m.par["surv.z"] * z  # linear predictor
    p <- 1/(1+exp(-linear.p))                            # logistic transformation to probability
    return(p)
}

## Probability of flowering function, logistic regression
p_bz <- function(z, m.par)
{
    linear.p <- m.par["flow.int"] + m.par["flow.z"] * z      # linear predictor
    p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
    return(p)
}

## Seedling production function
## Expected seedling production of a size-z plant that flowers
b_z <- function(z, m.par)
{
    heads <- m.par["heads.int"] + m.par["heads.z"]*exp(z);  
    return(m.par["eps.pi"]*heads);
}

## Recruit size pdf
c_z1 <- function(z1, m.par)
{
    mean <- m.par["rcsz.int"]
    sd <- m.par["rcsz.sd"]
    p.deRecr <- dnorm(z1, mean = mean, sd = sd)              # pdf of a size z1 recruit
    return(p.deRecr)
}


displacement<-function(n,m.par) {
    #r <-rlnorm(n,meanlog=m.par["logr.mean"],sdlog=m.par["logr.sd"])
    r <- rinvGauss(n,nu=m.par["r.nu"],lambda=m.par["r.lambda"]); 
    theta=runif(n,0,2*pi);
    return(r*cbind(cos(theta),sin(theta)))
}

#######################
# Section 2
#######################
Waldmgf <- function(s,m.par) {
    nu <- m.par["r.nu"];
    lambda <- m.par["r.lambda"];
    t1 <- (lambda/nu); 
    t2 <- 2*(nu^2)*s/lambda; 
    mgf <- exp(t1*(1-sqrt(1-t2)));
    return(mgf)
}    
    
margWaldmgf <- function(s,m.par) {
  out=integrate(function(theta) EWaldmgf(s*cos(theta),m.par), lower=0,upper=pi)$value;
  return(out/pi)
}
margWaldmgf <- Vectorize(margWaldmgf,vectorize.args="s"); 

s.max <- function(m.par) {
  nu <- m.par["r.nu"];
  lambda <- m.par["r.lambda"];
  return(lambda/(2*nu*nu)); 
}  

Kernel <- function(z1,z,m.par) {
      pb <- p_bz(z,m.par); 
      s_z(z,m.par)*( (1-pb)*g_z1z(z1,z,m.par) + pb*b_z(z,m.par)*c_z1(z1,m.par))
}

Fkernel <- function(z1,z,m.par) {
  pb <- p_bz(z,m.par); 
  s_z(z,m.par)*pb*b_z(z,m.par)*c_z1(z1,m.par)
}

Pkernel <- function(z1,z,m.par) {
  pb <- p_bz(z,m.par); 
  s_z(z,m.par)*(1-pb)*g_z1z(z1,z,m.par)
}


m<- 100; L <- -3; U <- 6;
h = (U - L)/m
meshpts = L + ((1:m) - 1/2) * h
Kmat <- Fmat <- Pmat <- matrix(NA,m,m); 
for(i in 1:m){
  for(j in 1:m){
    Kmat[i,j]=Kernel(meshpts[i],meshpts[j],m.par.true)
    Fmat[i,j]=Fkernel(meshpts[i],meshpts[j],m.par.true)
    Pmat[i,j]=Pkernel(meshpts[i],meshpts[j],m.par.true)
  }
}
Kmat=h*Kmat; Pmat=h*Pmat; Fmat=h*Fmat; 

cspeed=function(s,m.par) {
    Mkernel=EWaldmgf(s,m.par)*Fmat + Pmat;
    lam.s = Re(eigen(Mkernel)$values[1]);
    return( (1/s)*log(lam.s))
}
Cspeed <- Vectorize(cspeed,vectorize.args="s"); 

#plot(function(s) Cspeed(s,m.par.true),0.01,0.999*s.max(m.par.true), 
#     ylab="C",ylim=c(0,50)); 

plot(function(s) Cspeed(s,m.par.true),0.01,.15,
     ylab="C",ylim=c(0,50));

    
### Use empirical kernel
nd <- 250000; 
D <- displacement(nd,m.par.true);
xdisp <- D[,1];

EWaldmgf <- function(s,m.par) {
    (1/nd)*sum(exp(s*xdisp))
}

EWaldmgf <- Vectorize(EWaldmgf,vectorize.args="s");

plot(function(s) EWaldmgf(s,m.par.true)/s,0.01,.12,ylab="M/s");
    

