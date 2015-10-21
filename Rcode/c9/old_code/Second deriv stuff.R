## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IBM with individual variation in flowering intercept to illustrate evolutionary dynamics
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
setwd("~/Repos/ipm_book/Rcode/c2")
#setwd("c:/Repos/ipm_book/Rcode/c2")
#setwd("~/Repos/ipm_book/Rcode/c9/old_code")
#setwd("c:/Repos/ipm_book/Rcode/c9/old_code")

source("Monocarp Demog Funs.R");
source("../utilities/Standard Graphical Pars.R");

source("../c9/old_code/Monocarp VarDynamics Funs.R"); 

# Set simulation parameters
init.pop.size <- 10000
Recr <- 7500
n.yrs <- 5000
beta.off.sd <- 0.05
init.beta.sd <- 0.1

#  This generates the simulation data, it takes a while so we use 
#  a saved version
load("../c9/SimDataSmallVar.Rdata")



########################################################################### 
# FIRST PLOT: IBM mean and variance of flowering intercept over time, 10 reps. 
###########################################################################

graphics.off(); dev.new(); 
set_graph_pars("panel4"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",
	ylab="Mean flowering intercept",xlim=c(0,n.yrs))
add_panel_label("a")
matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",
ylab="Variance flowering intercept",xlim=c(0,n.yrs))
add_panel_label("b")
#dev.copy2eps(file="~/Repos/ipm_book/c10/figures/OenotheraIBMevol.eps")

nBigMatrix.z <- 100; 
nBigMatrix.beta <- 100; 

#######################  solve beta(0) = -20 first
U.beta <- max(max.beta[1:5]) +  1 
L.beta <-  min(min.beta[1:5]) - 1
h.beta <- (U.beta - L.beta)/nBigMatrix.beta 
meshpts.beta <- L.beta + ((1:nBigMatrix.beta) - 1/2) * h.beta

U.z <- max(max.z[1:5]) + 0.5 
L.z <- min(min.z[1:5]) - 0.5 
h.z <- (U.z - L.z)/nBigMatrix.z
meshpts.z <- L.z + ((1:nBigMatrix.z) - 1/2) * h.z

########################   solve beta(0) = -30 next
U.beta <- max(max.beta[6:10]) +  1 
L.beta <-  min(min.beta[6:10]) - 1
h.beta <- (U.beta - L.beta)/nBigMatrix.beta 
meshpts.beta <- L.beta + ((1:nBigMatrix.beta) - 1/2) * h.beta

U.z <- max(max.z[6:10]) + 0.5 
L.z <- min(min.z[6:10]) - 0.5 
h.z <- (U.z - L.z)/nBigMatrix.z
meshpts.z <- L.z + ((1:nBigMatrix.z) - 1/2) * h.z



#Something is going wrong with the 2nd differential calculation - put resident to ESS and
#test for stabilizing selection at ESS

#This all now works - just needed to use undate w when doing the calculation.


nBigMatrix <- nBigMatrix.z; 

#beta.0 to test
beta.0 <- seq(-40,-10,length=100)
params <- m.par.true

#put ESS as resident
params["flow.int"] <- -24.92
#params["flow.int"] <- -30
params["p.r"] <- 1
equ.p.r <- 1/R0_calc(params)
params["p.r"] <- equ.p.r
delta <- 0.001

IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
b0 <- params["flow.int"]; 
params["flow.int"] <- b0 + delta
Kplus <- mk_K(nBigMatrix, params, L.z, U.z)
params["flow.int"] <- b0 - delta
Kminus <- mk_K(nBigMatrix, params, L.z, U.z)
params["flow.int"] <- b0 	

#Calculate the derivative of lambda with respect to beta
lambda                <-  Re(eigen(IPM.kernel$K,only.values = TRUE)$values[1])   
lambda.up           <- Re(eigen(Kplus$K,only.values = TRUE)$values[1])
lambda.down      <- Re(eigen(Kminus$K,only.values = TRUE)$values[1])
d.lambda.d.beta <- (lambda.up - lambda.down)/(2*delta)
d2.lam.d.beta2   <- (lambda.up - 2 * lambda + lambda.down)/(delta*delta)

cat("first deriv",d.lambda.d.beta," second deriv",d2.lam.d.beta2,"\n")

eigen.sys <- eigen(IPM.kernel$K)
w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))

eigen.sys <- eigen(Kplus$K)
wplus <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))

eigen.sys <- eigen(Kminus$K)
wminus <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))


store.w <- w


#if you change the next line so w does not change then the signs of the 2nd derivs are out
#wminus <- wplus  <-w
	
#Calculate the derivative of W - survival only - with respect to beta
	
	ave.surv <- sum(IPM.kernel$P %*% w)
	ave.surv.Up <- sum(Kplus$P %*% wplus)
	ave.surv.Down <- sum(Kminus$P %*% wminus)
	d.W.surv.d.beta <- (ave.surv.Up - ave.surv.Down)/(2*delta)
	d2.W.surv.d.beta2 <- (ave.surv.Up - 2 * ave.surv + ave.surv.Down)/(delta*delta)
	
cat("first deriv",d.W.surv.d.beta," second deriv",d2.W.surv.d.beta2,"\n")

#Calculate the derivative of W - fecundity only - with respect to beta
	ave.fec <- sum(IPM.kernel$F %*% w)
	ave.fec.Up <- sum(Kplus$F %*% wplus)
	ave.fec.Down <- sum(Kminus$F %*% wminus)
	d.W.fec.d.beta <- (ave.fec.Up - ave.fec.Down)/(2*delta)
	d2.W.fec.d.beta2 <- (ave.fec.Up - 2 * ave.fec + ave.fec.Down)/(delta*delta)

cat("first deriv",d.W.fec.d.beta," second deriv",d2.W.fec.d.beta2,"\n")

#plot fitness measures against beta.0

lam.fit <- rep(NA,100)
fec.fit  <- rep(NA,100)
surv.fit <- rep(NA,100)

for(i in 1:100){
	
	params["flow.int"] <- beta.0[i]
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	lam.fit[i]     <- Re(eigen(IPM.kernel$K,only.values = TRUE)$values[1])
	
	eigen.sys <- eigen(IPM.kernel$K)
	w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))
	

	#w <- store.w
	

	surv.fit[i] <- sum(IPM.kernel$P %*% w)
	fec.fit[i]   <- sum(IPM.kernel$F %*% w)
	
}

par(mfrow=c(2,2),bty="l",pty="s")

plot(beta.0,lam.fit,type="l")
abline(v = - 24.33)
lm(lam.fit~beta.0+I(beta.0*beta.0))
#second differential of lm
2*coef(lm(lam.fit~beta.0+I(beta.0*beta.0)))[3]
#what I get from finite differencing
d2.lam.d.beta2


plot(beta.0,surv.fit,type="l")
abline(v = - 24.33)
lm(surv.fit~beta.0+I(beta.0*beta.0))
#second differential of lm
2*coef(lm(surv.fit~beta.0+I(beta.0*beta.0)))[3]
#what I get from finite differencing
d2.W.surv.d.beta2

plot(beta.0,fec.fit,type="l")
abline(v = - 24.33)
lm(fec.fit~beta.0+I(beta.0*beta.0))
#second differential of lm
2*coef(lm(fec.fit~beta.0+I(beta.0*beta.0)))[3]
#what I get from finite differencing
d2.W.fec.d.beta2

plot(beta.0,lam.fit,type="l")
points(beta.0,surv.fit+fec.fit,col="red")


round(d.lambda.d.beta,6)

round(d.W.surv.d.beta + d.W.fec.d.beta,6)

round(d2.lam.d.beta2,6)

round(d2.W.surv.d.beta2 + d2.W.fec.d.beta2,6)


