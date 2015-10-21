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

m.par <- c(surv.int  =  -0.65,
                surv.z    =   0.75,
                flow.int  = -18.00,
                flow.z    =   6.9,
                grow.int  =   0.96,
                grow.z    =   0.59,
                grow.sd   =   0.67,
                rcsz.int  =   -.08, 
                rcsz.sd   =   0.76, 
                seed.int  =   1.00,
                seed.z    =   2.20,
                p.r       =   0.005325442)  #changed p.r so R0=1


##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets, say the true values and
## estimated ones
##

## Growth function, given you are size z now returns the pdf of size z1 next time

G_z1z <- function(z1, z, m.par)
{
    mean <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year
    sd <- m.par["grow.sd"]                                    # sd about mean
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

## Seed production function

b_z <- function(z, m.par)
{
    N <- exp(m.par["seed.int"] + m.par["seed.z"] * z)    # seed production of a size z plant
    return(N)
}

## Recruit size pdf

c_0z1 <- function(z1, m.par)
{
    mean <- m.par["rcsz.int"]
    sd <- m.par["rcsz.sd"]
    p.deRecr <- dnorm(z1, mean = mean, sd = sd)              # pdf of a size z1 recruit
    return(p.deRecr)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
P_z1z <- function (z1, z, m.par) {

    return((1 - p_bz(z, m.par)) * s_z(z, m.par) * G_z1z(z1, z, m.par))

}

## Define the fecundity kernel
F_z1z <- function (z1, z, m.par) {

    return( p_bz(z, m.par) * b_z(z, m.par) * m.par["p.r"] * c_0z1(z1, m.par))

}

mk_K <- function(m, m.par) {

	# upper and lower integration limits
	L <- min.size - 2
	U <- max.size + 4

	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
	F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}













