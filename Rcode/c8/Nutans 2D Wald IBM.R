## C. nutans 2D spread IBM, InverseGaussian dispersal kernel
## uses density() to estimate pdf of r


rm(list=ls(all=TRUE))
require(SuppDists); 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so formulae are easier to
## read. We'll use a 'model.term' scheme to name elements of the vector, except that
## eps.pi is the product of pi and epsilon; pi and epsilon are only used in the calculation of
## seedling production, which is Poisson with mean epsilon*pi*(expected number of flowering heads) 

m.par <- c(surv.int  =  -2.27,
                surv.z    =   0.57,
                flow.int  =  -2.107,
                flow.z    =   0.86,
                grow.int  =   2.751,
                grow.z    =   0.407,
                grow.sd.int   = 3.02, 
                grow.sd.z     = -0.228,
                rcsz.int  =  -0.771,
                rcsz.sd   =   1.719,
                heads.int  =  6.363,
                heads.z    =  0.0056,
                eps.pi    =   7.1*0.8,
                r.nu      =   3, 
                r.lambda =    1.5)

m.par.true <- m.par; 
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets

## Growth function, given you are size z now returns the pdf of size z1 next time
g_z1z <- function(z1, z, m.par)
{
    mean <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year
    sd <- m.par["grow.sd.int"]*exp(m.par["grow.sd.z"]*mean)       # sd about mean
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

sigma.d <- sinvGauss(m.par["r.nu"],m.par["r.lambda"])$Variance^0.5;

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Parameters controlling the IBM simulation
init.pop.size <- 1000
n.yrs <-150
n.reps <- 3
pop.size.t <- rcrit.t <- matrix(NA,n.yrs,n.reps)

for(rep in 1:n.reps) {
z <- rnorm(init.pop.size, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])
loc <- rep(1,init.pop.size)
X=cbind(displacement(init.pop.size,m.par),z,loc); 

## iterate the model 
yr <- 1;
while(yr <= n.yrs) {
    ## calculate population size
    pop.size <- nrow(X)

    ## generate binomial random number for surviving (or not) to the time of flowering. 
    ## This is a vector of 0's and 1's, with 1 = survive. 
    Surv <- rbinom(n=pop.size,prob=s_z(X[,3],m.par),size=1); 
    X <- X[Surv==1,];
  
    ## Store population size just before flowering time 
    pop.size <- nrow(X) 
    pop.size.t[yr,rep] <- pop.size
    
    ## Update locations of new recruits, flagged by X[,4]==0. 
    ## For these, the 1st two columns are parent location
    newRecruit <- (X[,4]==0);
    nNew <- sum(newRecruit);
    X[newRecruit,1:2] <- X[newRecruit,1:2]+displacement(nNew,m.par);
    X[newRecruit,4] <- rep(1,nNew); 
    
    ## Estimate density as a function of r 
    r=sqrt(X[,1]^2+X[,2]^2);
    out <- density(r);
    denfun <- approxfun(out$x,length(r)*out$y,rule=2);
    if(yr%%5==1) plot(function(r) denfun(r)/r,5,max(r),yaxs="i",main=yr);
    density <- denfun(r)/r; 
    rfar <- r[density>0.1];
    rcrit.t[yr,rep]=max(rfar);
 
    ## apply density-dependent mortality 
    survivalProb <- 1/(1+0.5*pmax(denfun(r)/r-2,0)); 
    Surv <- rbinom(n=length(r),prob=survivalProb,size=1);
    X <- X[Surv==1,];
    

    ## generate binomial random number for the probability of flowering given your
    ## size z. This is a vector of 0's and 1's, you get a 1 if you flower
    Repr <- rbinom(n=nrow(X), prob=p_bz(X[,3], m.par), size=1)

    ## number of plants that flowered
    num.Repr <- sum(Repr)

    ## Calculate seedling production
    ## We assume a plant make a Poisson distributed number of seedlings 
    ## with mean given by b_z
    ## Note: decrease seedlings to slow down population growth
    Seedlings <- rpois(num.Repr, b_z(X[Repr==1,3], m.par))
    
    # Parent locations of the seedlings 
    Seedlings.x <- rep(X[Repr==1,1],Seedlings); 
    Seedlings.y <- rep(X[Repr==1,2],Seedlings); 
    
    ## total number of recruits
    Recr <- sum(Seedlings); 

    ## generate new recruit sizes
    ## rnorm generated normally distributed random numbers
    Rcsz <- rnorm(Recr, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])
    
    ## Create matrix of recruit sizes and locations. Last column
    ## is 0 to flag: (x,y) are the location of my parent. 
    R <- cbind(Seedlings.x,Seedlings.y,Rcsz,rep(0,Recr)); 
    
 
    ## Kill off the flowering plants
    X <- X[Repr==0,];

    ## Let the non-flowering plants grow
    mean.z1 <- m.par["grow.int"] + m.par["grow.z"] * X[,3]           # mean size next year
    sd.z1 <- m.par["grow.sd.int"]*exp(m.par["grow.sd.z"]*mean.z1)       # sd about mean

    X[,3] <- rnorm(n = nrow(X), mean.z1,sd.z1)

    ## combine survivors and recruits
    X <- rbind(X,R)
    cat(paste(yr, pop.size.t[yr,rep], "\n", sep=" "))
    yr <- yr+1

}
}

## Plot the growth of occupied area over time 
matplot(1:n.yrs,rcrit.t,type="l",lty=c(1,2,3,4),col="black",xlab="Years",
ylab="Radius of area occupied"); 

# Estimate spread rate from the last 25 years of spread. 
t.end <- seq(n.yrs-25,n.yrs,by=1);
spread.rate <- numeric(n.reps); 
for(rep in 1:n.reps) {
    spread.rate[rep] <- lm(rcrit.t[t.end,rep]~t.end)$coef[2]
}
cat(spread.rate); 
#dmax=max(abs(X[,1:2])); 
#plot(X[,1],X[,2],pch=16,xlim=c(-dmax,dmax),ylim=c(-dmax,dmax),cex=0.5); 
