## IBM to generate the data for the simple density dependent ungulate example
## Density dependence appears in adult survival, recruitment probabiity, and 
## the size distribution of recruits. 

## The parameters and life cycle correspond to the Soay sheep of St Kilda, 
## estimated from 26 years of data (1985-2010)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Define the demographic functions and parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so the formulae, below, are easier to
## read. We'll use 'model.term' scheme to name elements of the vector

m.par.true <- c(## survival (density dependent)
                surv.int  =  1.06e+0,
                surv.z    =  2.09e+0,
                surv.Nt   = -1.80e-2,
                ## growth (NOT density dependent)
                grow.int  =  1.41e+0,
                grow.z    =  5.57e-1,
                grow.sd   =  7.99e-2,
                ## reproduce or not (NOT density dependent)
                repr.int  = -7.23e+0,
                repr.z    =  2.60e+0,                
                ## recruit or not (density dependent)
                recr.int  =  4.43e+0,
                recr.Nt   = -9.19e-3,
                ## recruit size (density dependent)
                rcsz.int  =  5.40e-1,
                rcsz.z    =  7.10e-1,
                rcsz.Nt   = -6.42e-4,
                rcsz.sd   =  1.59e-1)

##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different parameter sets
##

## Growth function, given you are size z now returns the pdf of size z1 next time
g_z1z <- function(z1, z, m.par) {
    mean <- m.par["grow.int"] + m.par["grow.z"] * z           
    sd <- m.par["grow.sd"]                                    
    p.den.grow <- dnorm(z1, mean = mean, sd = sd)
    return(p.den.grow)
}

## Survival function, logistic regression
s_z <- function(z, Nt, m.par) {
    linear.p <- m.par["surv.int"] + m.par["surv.Nt"] * Nt + m.par["surv.z"] * z 
    p <- 1/(1+exp(-linear.p))                                                 
    return(p)
}

## Reproduction function, logistic regression
pb_z <- function(z, m.par) {
    linear.p <- m.par["repr.int"] + m.par["repr.z"] * z       
    p <- 1/(1+exp(-linear.p))                                 
    return(p)
}

## Recruitment function (N.B - from birth in spring to first summer), logistic regression
pr_z <- function(Nt, m.par) {
    linear.p <- m.par["recr.int"] + m.par["recr.Nt"] * Nt
    p <- 1/(1+exp(-linear.p))
    return(p)
}

## Recruit size function
c_z1z <- function(z1, z, Nt, m.par) {
    mean <- m.par["rcsz.int"] + m.par["rcsz.Nt"] * Nt + m.par["rcsz.z"] * z
    sd <- m.par["rcsz.sd"]
    p.den.rcsz <- dnorm(z1, mean = mean, sd = sd)
    return(p.den.rcsz)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
P_z1z <- function (z1, z, Nt, m.par) {

    return( s_z(z, Nt, m.par) * g_z1z(z1, z, m.par) )

}

## Define the reproduction kernel
F_z1z <- function (z1, z, Nt, m.par) {

    return( 0.5* s_z(z, Nt, m.par) * pb_z(z, m.par) * pr_z(Nt, m.par) * c_z1z(z1, z, Nt, m.par) )

}

## Build the iteration matrix for the discretized kernel
mk_K <- function(Nt, m, m.par, L, U) {
	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, P_z1z, Nt = Nt, m.par = m.par))
	F <- h * (outer(meshpts, meshpts, F_z1z, Nt = Nt, m.par = m.par))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build the Jacobian kernel 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Iteration matrix for Jacobian at an equilibrium nbar 
mk_J <- function(nbar, m, m.par, L, U,eps=0.01) {
    h <- (U-L)/m;  Nbar <- h*sum(nbar) # compute total population  

	# first compute Q = dK/dN %*% nbar
    dK <- (mk_K(Nbar+0.5*eps, m, m.par, L, U)$K - mk_K(Nbar-0.5*eps, m, m.par, L, U)$K)/eps
    Q <- dK%*%nbar 
    
    # Next the kernel R(z',z) = Q(z')W(z), i.e. R[i,j]=Q[i]*W[j]
    # In this case W==1 so R[i,] = Q[i] 
    R <- matrix(0,m,m)
    for(i in 1:m) R[i,] <- Q[i] 
    
    # Finally, Q = R + K; note factor of h to give iteration matrix
    # for J rather than kernel values, to parallel the output of mk_K 
    J <- mk_K(Nbar, m, m.par, L, U)$K + h*R 
    return(J)
}    
    
    
    

