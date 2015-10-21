## IBM to generate the data for the simple recruitment limited IPM - there's no recruitment limitation 1st
## Size is on a log scale so don't worry about the -ive ones!

## The parameters and life cycle follow Kachi and Hirose J.Ecol. 1985, see Rees and Rose PRSB 2002, but for
## simplicity we have used a logistic function for survival rather than the linear function with an upper bound
## as in the previous papers

## MORE COMMENTS STILL TO BE ADDED

rm(list=ls(all=TRUE))

library(doBy)
#set.seed(53241986)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##
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
                grow.sd.int   = 1.51, 
                grow.sd.z     = -0.114,
                rcsz.int  =  -0.771,
                rcsz.sd   =   1.719,
                heads.int  =  6.363,
                heads.z    =  0.0056,
                eps.pi    =   7.1*0.6)

##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets, say the true values and
## estimated ones
##

## Growth function, given you are size z now returns the pdf of size z1 next time
g_z1z <- function(z1, z, m.par)
{
    mean <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year
    #sd <- m.par["grow.sd.int"]/(1+m.par["grow.sd.z"]*z)       # sd about mean
    sd <- m.par["grow.sd.int"]*exp(m.par["grow.sd.z"]*z)      # sd about mean
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##
## Parameters controlling the IBM simulation
##

m.par <- m.par.true

init.pop.size <- 500
n.yrs <-50
z <- rnorm(init.pop.size, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])

## Calculate initial pop size and mean size

pop.size.t <- numeric(0)
mean.z.t <- mean(z)
mean.fl.z.t <- NA

## iterate the model using the "true" parameters and store data in a data.frame
yr <- 1
while(yr != n.yrs & max(pop.size.t) < 25000) {

    ## calculate population size
    pop.size <- length(z)

    ## generate binomial random number for surviving (or not) to the time of flowering. 
    ## This is a vector of 0's and 1's, with 1 = survive. 
    Surv <- rbinom(n=pop.size,prob=s_z(z,m.par),size=1); 
    z <- z[Surv==1];    # only keep the survivors 
    pop.size <- length(z); 
    pop.size.t <- c(pop.size.t, length(z))
    
    ## generate binomial random number for the probability of flowering given your
    ## size z. This is a vector of 0's and 1's, you get a 1 if you flower
    Repr <- rbinom(n=pop.size, prob=p_bz(z, m.par), size=1)

    ## number of plants that flowered
    num.Repr <- sum(Repr)

    ## calculate seedling production
    Seedlings <- rep(0, pop.size)

    ## we'll assume a plant make a Poisson distributed number of seedlings 
    ## with mean given by b_z
    ## Note: decrease seedlings to slow down population growth
    Seedlings[Repr==1] <- rpois(num.Repr, b_z(z[Repr==1], m.par))

    ## total number of recruits
    Recr <- sum(Seedlings); 

    ## generate new recruit sizes
    ## rnorm generated normally distributed random numbers
    Rcsz <- rnorm(Recr, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])

    ## let the non-flowering plants grow
    z <- z[Repr==0]; # kill off the flowering plants

    mean.z1 <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year
    sd.z1 <- m.par["grow.sd.int"]/(1+m.par["grow.sd.z"]*z)       # sd about mean

    z1<- rnorm(n = length(z), mean.z1,sd.z1)

    ## create new population size vector
    z <- c(Rcsz, z1)
    cat(paste(yr, pop.size.t[yr], "\n", sep=" "))
    yr <- yr+1

}
plot(log(pop.size.t)); 
