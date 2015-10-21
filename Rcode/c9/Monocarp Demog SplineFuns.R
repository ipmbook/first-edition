#################################################################################### 
#  Functions to create a version of the chapter 2 monocarp model where the flowering
#  strategy is a function-valued trait specified by a spline basis and coefficients 
#  This relies on Rcode/c2/Monocarp Demog Funs.R; functions defined there and used
#  without change are not repeated here.
####################################################################################

## Probability of flowering function 
## Input argument B must be a basis object of type 'bspline'
## created by create.bspline.basis() in the fda package  
p_bz_spline <- function(z, B, cj) {
    linear.p <- eval.basis(B,z)%*%cj;      # linear predictor
    p <- 1/(1+exp(-linear.p))              # logistic transformation to probability
    return(p) 
}

###################################################################################
## Functions to build IPM kernels. 
###################################################################################

## Define the fecundity kernel
F_z1z_spline <- function (z1, z, m.par, B, cj) {
    return( p_bz_spline(z, B, cj) * b_z(z, m.par) * m.par["p.r"] * c_0z1(z1, m.par))
}

## Define the survival kernel
P_z1z_spline <- function(z1, z, m.par, B, cj) {
    return((1 - p_bz_spline(z, B,cj)) * s_z(z, m.par) * G_z1z(z1, z, m.par))

}

## Define the fecundity kernel, to operate on vectors z1 and z 
## This is much faster than using outer() when building an iteration matrix 
F_z1z_spline_vec <- function (z1, z, m.par, B, cj) {
	a <- matrix(c_0z1(z1, m.par),ncol=1); 
	b <- matrix(p_bz_spline(z, B, cj) * b_z(z, m.par),nrow=1); 
	return(m.par["p.r"]*(a%*%b)); 
}

## Make the iteration matrices 
mk_K_spline <- function(m, m.par, L, U, B, cj) {
	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, P_z1z_spline, m.par = m.par, B = B, cj = cj))
	F <- h* F_z1z_spline_vec(meshpts, meshpts, m.par=m.par, B=B, cj=cj)
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}





