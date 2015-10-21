## Parameters and functions for the Carlina recruitment limited IPM 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Define the demographic functions and parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so the formulae, below, are easier to
## read. We'll use 'model.term' scheme to name elements of the vector

m.par.true <- c(## survival
                surv.int       =  -2.284390,
                surv.z          =   0.9040202,
                ## flowering
                flow.int       = -16.190,
                flow.z          =   3.883,
                ## growth
                grow.int         =   1.135735,
                grow.z            =   0.734684,
                grow.sd          =   0.2846671,
                ## recruit size
                rcsz.int       =   3.162432, 
                rcsz.sd   =   0.4988126,
                ## seed size
                seed.int  =   -4.00,
                seed.z    =   2.00,
                ## recruitment probability
                p.r       =   0.007)  


m.par.sd.true <- c(## survival
                surv.int.sd        =   1.161576,
                surv.z.sd           =   0.4136052,
                ## flowering
                flow.int.sd       =  1.029,
                flow.z.sd          =   0.0,
                ## growth
                grow.int.sd         =   0.1928372,
                grow.z.sd            =   0.1329790,
                grow.sd.sd          =  0.0,
                ## recruit size
                rcsz.int.sd       =   0.2706014, 
                rcsz.sd.sd       =   0.0,
                ## seed size
                seed.int.sd  =   0.0,
                seed.z.sd     =   0.0,
                ## recruitment probability
                p.r.sd           =   0.00) 
                
#actual numbers of recruits per year

nrec<-c(20,42,12,17,8,19,58,45,44,2,56,25,75,92,94,6,4,34,104)

store.nrec<-nrec;

                
#Variance-covariance matrix for random effects estimated simultaneously using WinBugs

VarCovar.grow.rec<-matrix(c(0.03718619,0.04022101,0.04022101,0.07322512),nrow=2)

chol.VarCovar.grow.rec=chol(VarCovar.grow.rec)

store.VarCovar.grow.rec=VarCovar.grow.rec

mean.grow.rec<-c(1.135735,3.162432)

##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets, say the true values and
## estimated ones
##

## Growth function, given you are size z now returns the pdf of size z1 next time

G_z1z <- function(z1, z, m.par)
{
    mu <- m.par["grow.int"] + m.par["grow.z"] * z              # mean size next year
    sig <- m.par["grow.sd"]                                    # sd about mean
    p.den.grow <- dnorm(z1, mean = mu, sd = sig)               # pdf that you are size z1 given you were size z
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

##NOTE there is no p.r in the F_z1z function as this is where density dependence acts and it needs to be 
##calculated each year for a density dependent model.

## Define the survival kernel
P_z1z <- function (z1, z, m.par) {

    return((1 - p_bz(z, m.par)) * s_z(z, m.par) * G_z1z(z1, z, m.par))

}

## Define the fecundity kernel
F_z1z <- function (z1, z, m.par) {

    return( s_z(z, m.par) * p_bz(z, m.par) * b_z(z, m.par)  * c_0z1(z1, m.par))

}

mk_K <- function(m, m.par, L, U) {
	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
	F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F, h=h))
}






