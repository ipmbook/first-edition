## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IBM with individual variation in flowering intercept to illustrate evolutionary dynamics
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
setwd("~/Repos/ipm_book/Rcode/c2")
source("Monocarp Demog Funs.R");

source("../utilities/Standard Graphical Pars.R");

# Set simulation parameters
setwd("~/Repos/ipm_book/Rcode/c10")

init.pop.size <- 10000
Recr <- 7500
n.yrs <-2000
beta.off.sd <- 0.5
init.beta.sd <- 0.1

#this generates the simulation data, it takes a while so we use 
#a saved version

# # mean.flow.ints <- matrix(NA,nrow=n.yrs,ncol=10)
# var.flow.ints <- matrix(NA,nrow=n.yrs,ncol=10)
# min.z <- rep(NA,10)
# max.z <- rep(NA,10)
# min.beta <- rep(NA,10)
# max.beta <- rep(NA,10)

# for(reps in 1:10){
	# init.mean.flow.int <- ifelse(reps<6, -20, -30)
	# source("Monocarp Simulate Evol IBM.R") 
    # mean.flow.ints[,reps] <- mean.flow.int
    # var.flow.ints[,reps] <- var.flow.int
    # min.z[reps] <- min.z.t
    # max.z[reps] <- max.z.t
    # min.beta[reps] <- min.flow.int
    # max.beta[reps] <- max.flow.int
    # cat(reps,"   ",init.mean.flow.int,"\n")
# }

#save(mean.flow.ints,var.flow.ints,min.z,max.z,min.beta,max.beta,file="SimData.Rdata")

load("SimData.Rdata")

set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Mean flowering intercept")
add_panel_label("a")
matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Variance flowering intercept")
add_panel_label("b")

dev.copy2eps(file="~/Repos/ipm_book/c10/figures/OenotheraIBMevol.eps")

#Implement a 2-dimensional IPM 

iterate_model <- function(params,n.iter,init.beta.mean,init.beta.sd) {
	
	#nt each column is the density of z for a given flowering intercept
    nt <- matrix(NA,nrow=nBigMatrix.z,ncol=nBigMatrix.beta)
    z.s <- dnorm(meshpts.z, mean = 2.9, sd = m.par.true["rcsz.sd"])
    flow.int.s <- dnorm(meshpts.beta, mean = init.beta.mean, sd = init.beta.sd)
    nt <- init.pop.size*matrix(outer(z.s,flow.int.s),nrow=nBigMatrix.z)
	
	#array of P matricies for each beta
    P.beta <- array(NA,c(nBigMatrix.beta,nBigMatrix.z,nBigMatrix.z))

	#array to store 
    Seeds.beta <- array(NA,c(nBigMatrix.beta,nBigMatrix.z))

	#Number of seeds by a size z plant
    Seeds_z <- function ( z, m.par) {
            return( p_bz(z, m.par) * b_z(z, m.par) )
     }

	for(i in 1:nBigMatrix.beta){
		params["flow.int"] <- meshpts.beta[i]
		P.beta[i,1:nBigMatrix.z,1:nBigMatrix.z] <- h.z * (outer(meshpts.z, meshpts.z, P_z1z, m.par = 			params))
		Seeds.beta[i,1:nBigMatrix.z] <- Seeds_z(meshpts.z, params)
	}

	#Inheritance kernel
	M <- matrix(NA,ncol=nBigMatrix.beta,nrow=nBigMatrix.beta)
	
	#Set up inheritance kernel
	for(i in 1:nBigMatrix.beta){
 		 M[,i] <- dnorm(meshpts.beta,mean=meshpts.beta[i],sd=beta.off.sd)
  	}

	#pdf of offspring sizes
	off.pdf <- c_0z1(meshpts.z,m.par.true)

	mean.beta <- rep(NA,n.iter)
	mean.beta[1] <- init.beta.mean
	var.beta <- rep(NA,n.iter)
	var.beta[1] <- init.beta.sd*init.beta.sd
	prop.Recr <- rep(NA,n.iter)


	index <- 1:nBigMatrix.beta

	#Iterate model

	for(gen in 2:n.iter){
		
	#calculate seeds produced by beta i  plants
	seeds.from.betai <- sapply(index,function(i) h.z*sum(Seeds.beta[i,] %*% nt[,i]))

	#redistribute seeds according to inheriatnce kernel
	seeds.with.betai <- h.beta * (M %*% seeds.from.betai)

	f.recruits.with.betai <- seeds.with.betai / sum(seeds.from.betai)

	nt1 <- sapply(index,function(i) P.beta[i,,] %*% nt[,i] + Recr * off.pdf * f.recruits.with.betai[i])
	
	n.surv <- sapply(index,function(i) P.beta[i,,] %*% nt[,i])
	
	n.surv <- sum(n.surv) * h.z * h.beta

	nt <- nt1

	mean.beta[gen] <- sum(apply(nt1,2,sum)*meshpts.beta/sum(nt1))
	var.beta[gen] <- sum(apply(nt1,2,sum)*meshpts.beta*meshpts.beta/sum(nt1)) - 	  mean.beta[gen]*mean.beta[gen]
	prop.Recr[gen] <- Recr/(Recr+n.surv)

cat(gen,"  ",mean.beta[gen],"   ",Recr/(Recr+n.surv),"\n")


}

return(list(mean.beta=mean.beta,var.beta=var.beta,prop.Recr=prop.Recr))
	
}


nBigMatrix.z <- 100; 
nBigMatrix.beta <- 100; 

#solve beta = -20 first

U.beta <- max(max.beta[1:5]) * 0.9
L.beta <-  min(min.beta[1:5]) * 1.1

h.beta <- (U.beta - L.beta)/nBigMatrix.beta 
meshpts.beta <- L.beta + ((1:nBigMatrix.beta) - 1/2) * h.beta

U.z <- max(max.z[1:5]) * 1.1 
L.z <- min(min.z[1:5]) * 1.1 

h.z <- (U.z - L.z)/nBigMatrix.z
meshpts.z <- L.z + ((1:nBigMatrix.z) - 1/2) * h.z

#nt each column is the density of z for a given flowering intercept


sol.2D.minus20 <- iterate_model(m.par.true,n.yrs,-20,init.beta.sd)

#solve beta = -30 next

U.beta <- max(max.beta[6:10]) * 0.9
L.beta <-  min(min.beta[6:10]) * 1.1

h.beta <- (U.beta - L.beta)/nBigMatrix.beta 
meshpts.beta <- L.beta + ((1:nBigMatrix.beta) - 1/2) * h.beta

U.z <- max(max.z[6:10]) * 1.1 
L.z <- min(min.z[6:10]) * 1.1 

h.z <- (U.z - L.z)/nBigMatrix.z
meshpts.z <- L.z + ((1:nBigMatrix.z) - 1/2) * h.z

#nt each column is the density of z for a given flowering intercept

sol.2D.minus30 <- iterate_model(m.par.true,n.yrs,-30,init.beta.sd)

set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Mean flowering intercept")
points(1:n.yrs,sol.2D.minus20$mean.beta,col="green",type="l")
points(1:n.yrs,sol.2D.minus30$mean.beta,col="blue",type="l")
abline(h=-24.92462,col="red")
add_panel_label("a")
matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Variance flowering intercept")
points(1:n.yrs,sol.2D.minus20$var.beta,col="green",type="l")
points(1:n.yrs,sol.2D.minus30$var.beta,col="blue",type="l")
add_panel_label("b")

###################################################################################################
#Approximate dynamics using Iwasa et al. 1991 Evolution
####################################################################################################

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


Approx_dynamics <- function(params,init.mean.flow.int,n.iter,delta,var.beta.t,prop.Recr) {


	params["flow.int"]<- init.mean.flow.int
	beta.t.mean <- rep(NA,n.iter)
	beta.t.var <- rep(NA,n.iter)
	S.t.mean <- rep(NA,n.iter) #this is mean S
	
	beta.t.mean[1] <- init.mean.flow.int
	beta.t.var[1] <- init.beta.sd*init.beta.sd
	S.t.mean[1] <- init.beta.sd*init.beta.sd + init.mean.flow.int*init.mean.flow.int
	

for(gen in 2:n.iter){

#Calculate p.r so the current mean strategy is at equilibrium lambda=1
	params["p.r"] <- 1
	equ.p.r <- 1/R0_calc(params)
	params["p.r"] <- equ.p.r

#Calculate the derivative of lambda with respect to beta
	params["flow.int"] <- params["flow.int"] + delta
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	lambda.up <- Re(eigen(IPM.kernel$K,only.values = TRUE)$values[1])
	params["flow.int"] <- params["flow.int"] - delta
	
	params["flow.int"] <- params["flow.int"] - delta
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	lambda.down <- Re(eigen(IPM.kernel$K,only.values = TRUE)$values[1])
	params["flow.int"] <- params["flow.int"] + delta

	d.lambda.d.beta <- (lambda.up - lambda.down)/(2*delta)
	
#Calculate the derivative of W - survival only - with respect to beta
	params["flow.int"] <- params["flow.int"] + delta
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	eigen.sys <- eigen(IPM.kernel$K)
	w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))
	ave.surv.Up <- sum(IPM.kernel$P %*% w)
	params["flow.int"] <- params["flow.int"] - delta
	
	params["flow.int"] <- params["flow.int"] - delta
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	eigen.sys <- eigen(IPM.kernel$K)
	w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))
	ave.surv.Down <- sum(IPM.kernel$P %*% w)
	params["flow.int"] <- params["flow.int"] + delta

	d.W.surv.d.beta <- (ave.surv.Up - ave.surv.Down)/(2*delta)
	
   IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	eigen.sys <- eigen(IPM.kernel$K)
	w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))
	W.surv <- sum(IPM.kernel$P %*% w)
	
#Calculate the derivative of W - fecundity only - with respect to beta
	params["flow.int"] <- params["flow.int"] + delta
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	eigen.sys <- eigen(IPM.kernel$K)
	w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))
	ave.fec.Up <- sum(IPM.kernel$F %*% w)
	params["flow.int"] <- params["flow.int"] - delta
	
	params["flow.int"] <- params["flow.int"] - delta
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	eigen.sys <- eigen(IPM.kernel$K)
	w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))
	ave.fec.Down <- sum(IPM.kernel$F %*% w)
	params["flow.int"] <- params["flow.int"] + delta

	d.W.fec.d.beta <- (ave.fec.Up - ave.fec.Down)/(2*delta)

	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	eigen.sys <- eigen(IPM.kernel$K)
	w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))
	W.fec <- sum(IPM.kernel$F %*% w)

#Calculate the mean intercept

	#We're cheating here as we're plugging in the var beta from the 2 D IPM var.beta.t[gen-1] 
	beta.t.mean[gen] <- beta.t.mean[gen-1] + d.lambda.d.beta * var.beta.t[gen-1]

#mean beta after survival
	beta.th.mean <- beta.t.mean[gen-1] + d.W.surv.d.beta * var.beta.t[gen-1] / W.surv
#Genetic variance of beta t-2	
	G.beta <- S.t.mean[gen-1] - beta.t.mean[gen-1] * beta.t.mean[gen-1]
#Genetic variance of beta^2
	G.S <- beta.t.mean[gen-1]^4 + 6 * beta.t.mean[gen-1] ^2 * G.beta + 3 * G.beta * G.beta - S.t.mean[gen-1] * S.t.mean[gen-1]
	
	#S.t.mean.Recr <- S.t.mean[gen-1] + G.S * d.lambda.d.beta / (2 * beta.t.mean[gen-1])
#mean S in recruits	
	S.t.mean.Recr <- S.t.mean[gen-1] + G.S * d.W.fec.d.beta / (2 * beta.t.mean[gen-1] * W.fec)
#mean S in survivors	
	S.t.mean.Surv <- S.t.mean[gen-1] + G.S * d.W.surv.d.beta / (2 * beta.t.mean[gen-1] * W.surv)
#mean S in population	
	S.t.mean[gen] <- prop.Recr[gen] * S.t.mean.Recr + (1-prop.Recr[gen]) * S.t.mean.Surv
#Variance of beta in recruits, note addition of mutational variance at end
	Var.Recr <- S.t.mean.Recr - beta.t.mean[gen]*beta.t.mean[gen] + beta.off.sd * beta.off.sd
#Variance of beta in survivors
	Var.Surv <- S.t.mean.Surv - beta.th.mean*beta.th.mean
#overall variance of beta at time t
	beta.t.var[gen] <- prop.Recr[gen] * Var.Recr + (1-prop.Recr[gen]) * Var.Surv +
	 				prop.Recr[gen]*(1-prop.Recr[gen])*(beta.t.mean[gen]-beta.th.mean)^2
		
	params["flow.int"] <- beta.t.mean[gen] 
	
	cat(gen,"  ",beta.t.mean[gen],"  ",beta.t.var[gen],"  ",var.beta.t[gen],"\n")
		
}

return(beta.t.mean)
}

nBigMatrix <- nBigMatrix.z; 

beta.minus20 <- Approx_dynamics(m.par.true,-20,n.yrs,0.1,sol.2D.minus20$var.beta,sol.2D.minus20$prop.Recr)
beta.minus30 <- Approx_dynamics(m.par.true,-30,n.yrs,0.1,sol.2D.minus30$var.beta,sol.2D.minus30$prop.Recr)


set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Mean flowering intercept")
points(1:n.yrs,sol.2D.minus20$mean.beta,col="green",type="l")
points(1:n.yrs,sol.2D.minus30$mean.beta,col="blue",type="l")
points(1:n.yrs,beta.minus20,col="red",type="l")
points(1:n.yrs,beta.minus30,col="red",type="l")

add_panel_label("a")
matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Variance flowering intercept")
points(1:n.yrs,sol.2D.minus20$var.beta,col="green",type="l")
points(1:n.yrs,sol.2D.minus30$var.beta,col="blue",type="l")
add_panel_label("b")

dev.copy2eps(file="~/Repos/ipm_book/c10/figures/OenotheraIBMevolPreds.eps")





























































