## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IBM with heritable variation in flowering intercept to illustrate evolutionary dynamics
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))
library(doBy)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c9",sep="")); 

source("../c2/Monocarp Demog Funs.R");
source("../utilities/Standard Graphical Pars.R");

# Set simulation parameters
Seeds_z <- function (z, m.par) {
            return( p_bz(z, m.par) * b_z(z, m.par) )
     }

R0_calc <- function (params) {
	IPM.true <- mk_K(nBigMatrix, params, L.z, U.z)
	# to keep close to the formulae in the text next we define the F and P iteration matricies
    P <- IPM.true$P;  F <- IPM.true$F;

    # Fundamental operator 
     N <- solve(diag(nBigMatrix)-P); 

    # Compute R0 as dominant eigenvalue of FN
    R <- F %*% N
    R0 <- abs(eigen(R)$values[1])
    return(R0)
}

lambda_beta <- function (beta.0,params){
	params["flow.int"] <- beta.0
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	lambda <- Re(eigen(IPM.kernel $K,only.values = TRUE)$values[1])
	return(lambda)
}

#Do some ESS calculations
nBigMatrix <- 100; 

U.z <-  6.5
L.z <-  -5.5

find_ESS <- function (params,tol) {
	repeat{
	
#calculate p.r so lambda=1 for resident - current strategy
	params["p.r"] <- 1
	equ.p.r <- 1/R0_calc(params)
	params["p.r"] <- equ.p.r

#find best invader
	opt.beta <- optimize(lambda_beta,lower=-30,upper=-15,params=params,maximum=TRUE)

#if best invader had lambda=1 to some tolerance we're done
	if(abs(opt.beta$objective-1)<tol) break()
	
#make best invader the resident
	params["flow.int"] <- opt.beta$maximum
	
	cat(opt.beta$objective,"   ",opt.beta$maximum,"\n")
	
}
return(opt.beta$maximum)
}

ESS.iterative.invasion <- find_ESS(m.par.true,0.0000001)

#Plot pips

n.betas<-250
beta.s <- seq(-30,-20,length=n.betas)
params <- m.par.true

lambdas.pip <- matrix(NA,ncol=n.betas,nrow=n.betas)

for(i in 1:n.betas){
	params["flow.int"] <- beta.s[i]

	params["p.r"] <- 1
	equ.p.r <- 1/R0_calc(params)
	params["p.r"] <- equ.p.r

    lambdas.pip[i,] <- sapply(beta.s,lambda_beta,params=params)
    cat(i,"\n"); 

}

dev.new(height=6,width=8);

set_graph_pars("panel1"); 
image(beta.s,beta.s,lambdas.pip<1, xlab=expression("Resident strategy " * italic(beta) [0]),
ylab=expression("Invader strategy " * italic(beta) [0]))

text(-25,-28,expression(italic(lambda)*" < 1"),cex=1.5)
text(-25,-22,expression(italic(lambda)*" < 1"),cex=1.5)
text(-28,-25,expression(italic(lambda)*" > 1"),cex=1.5)
text(-22,-25,expression(italic(lambda)*" > 1"),cex=1.5)

points(-28,-28,pch=19)
arrows(-28,-28,-28,-26,length=0.1)
points(-26,-26,pch=19)
arrows(-28,-26,-26.05,-26,length=0.1)
arrows(-26,-26,-26,-24.92,length=0.1)
arrows(-26,-24.92,-24.92,-24.92,length=0.1)
points(-24.92,-24.92,pch=19)

# dev.copy2eps(file="../../c9/figures/OenotheraPip-250.eps")

# IPM.true <- mk_K(nBigMatrix, m.par.true, L.z, U.z)
# 
# meshpts <- IPM.true$meshpts
# h <- diff(meshpts)[1]
# 
# lambda0 <- Re(eigen(IPM.true$K,only.values = TRUE)$values[1])

n.ints <- 200

R0s <- rep(NA,n.ints)
p.r.fl <- rep(NA,n.ints)
fl.ints <- seq(-24,-26,length=n.ints)
#fl.ints <- seq(2,10,length=100)

for(i in 1:n.ints) {
	params <- m.par.true
	params["flow.int"] <- fl.ints[i]
	
	R0s[i] <- R0_calc(params)
	
	params["p.r"] <- 1
	equ.p.r <- 1/R0_calc(params)
	p.r.fl[i] <- equ.p.r 
	
	}
	
ESS.fl.int <- fl.ints[which(R0s==max(R0s))]

ESS.fl.int

dev.new(height=4,width=8);
set_graph_pars("panel2"); 

plot(fl.ints, R0s,type="l", xlab=expression("Flowering intercept, "*beta[0]),ylab=expression(R[0]))
abline(v=ESS.iterative.invasion,col="red")
add_panel_label("a")
	
ESS.fl.int <- fl.ints[which(p.r.fl==min(p.r.fl))]

ESS.fl.int

plot(fl.ints, p.r.fl,type="l", xlab=expression("Flowering intercept, "*beta[0]),ylab=expression("Recruitment probability, "* p[r]))
abline(v=ESS.iterative.invasion,col="red")
add_panel_label("b")

# dev.copy2eps(file="~/Repos/ipm_book/c9/figures/OenotheraR0prOpt.eps")









