## IBM to generate the data for the simple ungulate example
## Size is on a log scale so don't worry about the -ive ones!

## The parameters and life cycle correspond to the Soay sheep of St Kilda, estimated from 26 years of data
## (1985-2010)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Define the demographic functions and parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so the formulae, below, are easier to
## read. We'll use 'model.term' scheme to name elements of the vector

m.par.true <- c(## survival
                surv.int  = -9.65e+0,
                surv.z    =  3.77e+0,
                ## growth 
                grow.int  =  1.41e+0,
                grow.z    =  5.57e-1,
                grow.sd   =  7.99e-2,
                ## reproduce or not
                repr.int  = -7.23e+0,
                repr.z    =  2.60e+0,                
                ## recruit or not
                recr.int  =  1.93e+0,
                ## recruit size
                rcsz.int  =  3.62e-1,
                rcsz.z    =  7.09e-1,
                rcsz.sd   =  1.59e-1)

##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets, say the true values and
## estimated ones
##

## Growth function, given you are size z now returns the pdf of size z1 next time
g_z1z <- function(z1, z, m.par){
    mean <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year
    sd <- m.par["grow.sd"]                                    # sd about mean
    p.den.grow <- dnorm(z1, mean = mean, sd = sd)             # pdf that you are size z1 given you were size z
    return(p.den.grow)
}

## Survival function, logistic regression
s_z <- function(z, m.par){
    linear.p <- m.par["surv.int"] + m.par["surv.z"] * z       # linear predictor
    p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
    return(p)
}

## Reproduction function, logistic regression
pb_z <- function(z, m.par){
    linear.p <- m.par["repr.int"] + m.par["repr.z"] * z       # linear predictor
    p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
    return(p)
}

## Recruitment function (N.B - from birth in spring to first summer), logistic regression
pr_z <- function(m.par) {
    linear.p <- m.par["recr.int"]                             # linear predictor
    p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
    return(p)
}

## Recruit size function
c_z1z <- function(z1, z, m.par){
    mean <- m.par["rcsz.int"] + m.par["rcsz.z"] * z           # mean size of recuits next year
    sd <- m.par["rcsz.sd"]                                    # sd about mean
    p.den.rcsz <- dnorm(z1, mean = mean, sd = sd)             # pdf that offspring are size z1 given you were size z
    return(p.den.rcsz)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
P_z1z <- function (z1, z, m.par) {
    return( s_z(z, m.par) * g_z1z(z1, z, m.par) )
}

## Define the reproduction kernel
F_z1z <- function (z1, z, m.par) {
    return( s_z(z, m.par) * pb_z(z, m.par) * (1/2) * pr_z(m.par) * c_z1z(z1, z, m.par) )
}

## Build the discretized kernel
mk_K <- function(m, m.par, L, U) {
	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
	F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}


