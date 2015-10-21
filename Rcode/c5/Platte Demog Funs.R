## Demographic functions for the Platte Thistle model with recruitment limitation
## and effects of weevil herbivory, from Rose et al. (2005). Effects of native
## herbivores are omitted. 



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Define the demographic functions and parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## We will only vary one parameter, the intercept for size-specific
## mean weevil egg load, which increases over time to reflect
## spread of weevils into the thistle population. 
m.par <- c(## herbivory intercept
                weevil.int  = -17)  

## Growth function, pdf of size z1 given size z 
G_z1z <- function(z1, z, m.par) {
    mu <- 0.83 + 0.69*z           # mean size next year
    sig <- sqrt(0.19)             # sd about mean
    p.den.grow <- dnorm(z1, mean = mu, sd = sig)  # pdf of z1 given current size z
    return(p.den.grow)
}

## Survival function, logistic regression

s_z <- function(z, m.par) {
    linear.p <- -0.62 + 0.85*z  # linear predictor
    p <- 1/(1+exp(-linear.p))   # logistic transformation to probability
    return(p)
}

## Probability of flowering function, logistic regression
p_bz <- function(z, m.par) {
    linear.p <- -10.22 + 4.25*z  # linear predictor
    p <- 1/(1+exp(-linear.p))    # logistic transformation to probability
    return(p)
}

## Seed production function of a size z plant, taking into account
## size dependent mean weevil load and negative exponential distribution
## of weevils. See Rose et al. (2005) Appendix, Ecological Archives E086-022-A1
b_z <- function(z, m.par) {
    eps <- exp(m.par["weevil.int"] + 1.71*z); 
    N <- exp(-0.55 + 2.02*z)/((1+eps/16)^0.32)
    return(N)
}

## Recruit size distribution 
c_0z1 <- function(z1, m.par) {
    pRecr <- dnorm(z1, mean = 0.75, sd = sqrt(0.17))  # pdf of recruit size z1
    return(pRecr)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build and iterate the IPM 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel function
P_z1z <- function (z1, z, m.par) {
	return((1 - p_bz(z, m.par)) * s_z(z, m.par) * G_z1z(z1, z, m.par))
}

## Define the total recruitment function
B_z <- function (z, m.par) {
    return( p_bz(z, m.par) * b_z(z, m.par) )
}

## function to build the iteration matrix for the survival kernel
## using midpoint rule with m meshpoints
mk_P <- function(m, m.par, L, U) {
	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
	return(list(meshpts = meshpts,P = P))
}

## Function to iterate the IPM using midpoint rule
## with mesh points meshpts
Iterate <- function(nt,meshpts,P,m.par) {
	h <- meshpts[2]-meshpts[1]
	N <- h*sum(nt*B_z(meshpts,m.par))
	recruits <- c_0z1(meshpts,m.par)*(N^0.67) 
	return(recruits + P%*%nt)	
}




















