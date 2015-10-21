## IBM to generate the data for the simple recruitment limited IPM - there's no recruitment limitation 1st
## Size is on a log scale so don't worry about the -ive ones!

## The parameters and life cycle follow Kachi and Hirose J.Ecol. 1985, see Rees and Rose PRSB 2002, but for
## simplicity we have used a logistic function for survival rather than the linear function with an upper bound
## as in the previous papers

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Define the demographic functions and parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so the formulae, below, are easier to
## read. We'll use 'model.term' scheme to name elements of the vector

m.par.true <- c(## survival
                surv.int  =  -0.65,
                surv.z    =   0.75,
                ## flowering
                flow.int  = -18.00,
                flow.z    =   6.9,
                ## growth
                grow.int  =   0.96,
                grow.z    =   0.59,
                grow.sd   =   0.67,
                ## recruit size
                rcsz.int  =   -.08, 
                rcsz.sd   =   0.76,
                ## seed size
                seed.int  =   1.00,
                seed.z    =   2.20,
                ## recruitment probability
                p.r       =   0.007)  


##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets, say the true values and
## estimated ones
##

## Growth function, given you are size z now returns the pdf of size z1 next time

G_z1z <- function(z1, z, m.par)
{
    mu <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year
    sig <- m.par["grow.sd"]                                    # sd about mean
    p.den.grow <- dnorm(z1, mean = mu, sd = sig)             # pdf that you are size z1 given you were size z
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

## Seed production function

b_z <- function(z, m.par)
{
    N <- exp(m.par["seed.int"] + m.par["seed.z"] * z)    # seed production of a size z plant
    return(N)
}

## Recruit size pdf

c_0z1 <- function(z1, m.par)
{
    mu <- m.par["rcsz.int"]
    sig <- m.par["rcsz.sd"]
    p.deRecr <- dnorm(z1, mean = mu, sd = sig)              # pdf of a size z1 recruit
    return(p.deRecr)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
P_z1z <- function (z1, z, m.par) {

    return((1 - p_bz(z, m.par)) * s_z(z, m.par) * G_z1z(z1, z, m.par))

}

## Define the fecundity kernel
F_z1z <- function (z1, z, m.par) {

    return( p_bz(z, m.par) * b_z(z, m.par) * m.par["p.r"] * c_0z1(z1, m.par))

}

mk_K <- function(m, m.par, L, U) {

	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
	F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}


# flexible size limits, defaults set for Oenothera model
mk_K_ceiling <- function(m, m.par, L, U, U1 = U) {
	# mesh points 
	h <- (U - L)/m;
	meshpts <- L + ((1:m) - 1/2) * h;
	P <- h * (outer(meshpts, pmin(meshpts,U1), P_z1z, m.par = m.par));
	F <- h * (outer(meshpts, pmin(meshpts,U1), F_z1z, m.par = m.par));
	K <- P + F;
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}

#Function to calculate mean and variance in reproductive output for a size z

get_mean_var_Repr <- function(init.z,n.samp) {

# initial population sizes and ages
z   <- rep(init.z,1000)
Repr.out <- NULL

repeat {

    ## calculate population size
    pop.size <- length(z)

    ## generate binomial random number for the probability of flowering, where the probability of flowering
    ## depends on your size z, this is a vector of 0's and 1's, you get a 1 if you flower
    Repr <- rbinom(n=pop.size, prob=p_bz(z, m.par.true), size=1)

    ## number of plants that flowered
    num.Repr <- sum(Repr)

    if(num.Repr>0) {
    	Seeds <- rpois(num.Repr, m.par.true["p.r"] * b_z(z[Repr==1],m.par.true))
    	Repr.out <- c(Repr.out,Seeds)
    }

    ## generate new recruit sizes
    ## rnorm generated normally distributed random numbers
    Rcsz <- rep(init.z,100)

    ## for the non-reproductive plants generate random number for survival
    Surv <- rep(NA, pop.size)
    Surv[Repr==0] <- rbinom(n = pop.size - num.Repr, prob = s_z(z[Repr==0], m.par.true), size = 1)
    num.die <- sum(Surv==0, na.rm=TRUE)
    
    if(num.die>0) Repr.out <- c(Repr.out,rep(0,num.die))

    ## index for individuals that did not flower and survived
    i.subset <- which(Repr==0 & Surv==1)

    ## let them grow
    E.z1 <- m.par.true["grow.int"]+m.par.true["grow.z"]*z[i.subset]
    z1 <- rnorm(n = pop.size - num.Repr - num.die, mean = E.z1, sd = m.par.true["grow.sd"])
    
    z <- c(Rcsz, z1)
	 
    if(length(Repr.out)>n.samp) break

}

return(c(mean(Repr.out),var(Repr.out),mean(Repr.out>0),var(Repr.out>0)))

}










