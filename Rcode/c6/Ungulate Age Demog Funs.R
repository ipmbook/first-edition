## IBM to generate the data for the age-structured ungulate example. Size is again on a log scale

## The parameters and life cycle correspond to the Soay sheep of St Kilda, estimated from 26 years of data
## (1985-2010)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Define the demographic functions and parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so the formulae, below, are easier to
## read. We'll use 'model.term' scheme to name elements of the vector

m.par.true <- c(## survival
  surv.int  = -1.70e+1,
  surv.z    =  6.68e+0,
  surv.a    = -3.34e-1,
  ## growth 
  grow.int  =  1.27e+0,
  grow.z    =  6.12e-1,
  grow.a    = -7.24e-3,
  grow.sd   =  7.87e-2,
  ## reproduce or not
  repr.int  = -7.88e+0,
  repr.z    =  3.11e+0,                
  repr.a    = -7.80e-2,                
  ## recruit or not
  recr.int  =  1.11e+0,
  recr.a    =  1.84e-1,
  ## recruit size
  rcsz.int  =  3.62e-1,
  rcsz.z    =  7.09e-1,
  rcsz.sd   =  1.59e-1)

##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets, say the true values and
## estimated ones. Note the addition of the 'a' argument for the new functions
##

## Growth function, given you are size z and age a now returns the pdf of size z1 next time
g_z1z <- function(z1, z, a, m.par){
  mean <- m.par["grow.int"] + m.par["grow.z"] * z + m.par["grow.a"] * a
  sd <- m.par["grow.sd"]
  p.den.grow <- dnorm(z1, mean = mean, sd = sd)
  return(p.den.grow)
}

## Survival function, logistic regression
s_z <- function(z, a, m.par){
  linear.p <- m.par["surv.int"] + m.par["surv.z"] * z + m.par["surv.a"] * a 
  p <- 1/(1+exp(-linear.p))                           
  return(p)
}

## Reproduction function, logistic regression
pb_z <- function(z, a, m.par){
  if (a==0) {
    p <- 0
  } else {
    linear.p <- m.par["repr.int"] + m.par["repr.z"] * z + m.par["repr.a"] * a
    p <- 1/(1+exp(-linear.p))
  }
  return(p)
}

## Recruitment function, logistic regression
pr_z <- function(a, m.par) {
  linear.p <- m.par["recr.int"] + m.par["recr.a"] * a
  p <- 1/(1+exp(-linear.p))
  return(p)
}

## Recruit size function
c_z1z <- function(z1, z, m.par){
  mean <- m.par["rcsz.int"] + m.par["rcsz.z"] * z
  sd <- m.par["rcsz.sd"]                          
  p.den.rcsz <- dnorm(z1, mean = mean, sd = sd)   
  return(p.den.rcsz)
}

