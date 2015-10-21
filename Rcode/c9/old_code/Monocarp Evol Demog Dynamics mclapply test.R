## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IBM with individual variation in flowering intercept to illustrate evolutionary dynamics
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
library(parallel)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
setwd("~/Repos/ipm_book/Rcode/c2")
source("Monocarp Demog Funs.R");

source("../utilities/Standard Graphical Pars.R");

# Set simulation parameters
setwd("~/Repos/ipm_book/Rcode/c10")

init.pop.size <- 10000
Recr <- 7500
n.yrs <-500
beta.off.sd <- 0.5
init.beta.sd <- 0.1

#this generates the simulation data, it takes a while so we use 
#a saved version

# mean.flow.ints <- matrix(NA,nrow=n.yrs,ncol=10)
# var.flow.ints <- matrix(NA,nrow=n.yrs,ncol=10)
# min.z <- rep(NA,10)
# max.z <- rep(NA,10)
# min.beta <- rep(NA,10)
# max.beta <- rep(NA,10)
# 
# for(reps in 1:10){
# 	init.mean.flow.int <- ifelse(reps<6, -20, -30)
# 	source("Monocarp Simulate Evol IBM.R") 
#     mean.flow.ints[,reps] <- mean.flow.int
#     var.flow.ints[,reps] <- var.flow.int
#     min.z[reps] <- min.z.t
#     max.z[reps] <- max.z.t
#     min.beta[reps] <- min.flow.int
#     max.beta[reps] <- max.flow.int
#     cat(reps,"   ",init.mean.flow.int,"\n")
# }
# 
# save(mean.flow.ints,var.flow.ints,min.z,max.z,min.beta,max.beta,file="SimData.Rdata")

load("SimData.Rdata")

set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Mean flowering intercept")
add_panel_label("a")
matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Variance flowering intercept")
add_panel_label("b")

#dev.copy2eps(file="~/Repos/ipm_book/c10/figures/OenotheraIBMevol.eps")

#Implement a 2-dimensional IPM 


iterate_model <- function(params,n.iter,init.beta.mean,init.beta.sd,meshpts.beta,meshpts.z) {

    nt <- matrix(NA,nrow=nBigMatrix.z,ncol=nBigMatrix.beta)
    nt1 <- nt
    n.surv.i <- nt
    z.s <- dnorm(meshpts.z, mean = 2.9, sd = m.par.true["rcsz.sd"])
    flow.int.s <- dnorm(meshpts.beta, mean = init.beta.mean, sd = init.beta.sd)
    nt <- init.pop.size*matrix(outer(z.s,flow.int.s),nrow=nBigMatrix.z)
	
	## Matrix to store distribution of flowering intercepts 
	nBeta <- length(meshpts.beta);
	betaDist <- matrix(0, nBeta,n.iter)
	betaDist[,1] = apply(nt,2,sum); betaDist[,1]=betaDist[,1]/sum(betaDist[,1])
	
	#array of P matrices for each beta
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
	for(i in 1:nBigMatrix.beta){
 		 M[,i] <- dnorm(meshpts.beta,mean=meshpts.beta[i],sd=beta.off.sd)
  	}

	#pdf of offspring sizes
	off.pdf <- c_0z1(meshpts.z,m.par.true)
	beta2 <- meshpts.beta*meshpts.beta
	mean.beta <- rep(NA,n.iter)
	mean.beta[1] <- init.beta.mean
	var.beta <- rep(NA,n.iter)
	var.beta[1] <- init.beta.sd*init.beta.sd
	prop.Recr <- rep(NA,n.iter)
	mean.beta2 <- rep(NA,n.iter)
	mean.beta2[1] <- sum(flow.int.s*beta2)/sum(flow.int.s)
	var.beta2 <- rep(NA,n.iter);
	
	########### Compute var(beta2) using formula Var=E[(X-Xbar)^2], not Var=E[X^2]-E[X]^2. 
	beta2Dev <- beta2-mean.beta2[1]; 
	var.beta2[1] <- sum(flow.int.s*beta2Dev^2)/sum(flow.int.s)
	
	index <- 1:nBigMatrix.beta
	
	#Iterate the model
	for(gen in 2:n.iter){
		
	#calculate seeds produced by beta_i plants
	seeds.from.betai <- sapply(index,function(i) h.z*sum(Seeds.beta[i,] %*% nt[,i]))
	#seeds.from.betai <- simplify2array(mclapply(index,function(i) h.z*sum(Seeds.beta[i,] %*% nt[,i])))
	#for(i in index) seeds.from.betai[i] <- h.z*sum(Seeds.beta[i,] %*% nt[,i])

	#redistribute seeds according to inheritance kernel
	seeds.with.betai <- h.beta * (M %*% seeds.from.betai)

	f.recruits.with.betai <- seeds.with.betai/sum(seeds.from.betai)

	nt1 <- sapply(index,function(i) P.beta[i,,] %*% nt[,i] + Recr * off.pdf * f.recruits.with.betai[i])
	#nt1 <- simplify2array(mclapply(index,function(i) P.beta[i,,] %*% nt[,i] + Recr * off.pdf * f.recruits.with.betai[i]))
	#for(i in index) nt1[,i] <- P.beta[i,,] %*% nt[,i] + Recr * off.pdf * f.recruits.with.betai[i]
	
	n.surv.i <- sapply(index,function(i) P.beta[i,,] %*% nt[,i])
	#n.surv.i <- simplify2array(mclapply(index,function(i) P.beta[i,,] %*% nt[,i]))
	#for(i in index) n.surv.i[,i] <- P.beta[i,,] %*% nt[,i]
	n.surv <- sum(n.surv.i) * h.z * h.beta
	nt <- nt1
	
	prop.Recr[gen] <- Recr/(Recr+n.surv)
	betaDist[,gen] = apply(nt,2,sum); betaDist[,gen] <- betaDist[,gen]/sum(betaDist[,gen]); 
	betaFreq <- betaDist[,gen] 
	
	#### Compute variances of beta and beta^2 using formula Var=E[(X-Xbar)^2] 
	mean.beta[gen] <- sum(betaFreq*meshpts.beta); 
	betaDev <- meshpts.beta - mean.beta[gen];
	var.beta[gen] <- sum(betaFreq*betaDev^2); 
	
	mean.beta2[gen] <- sum(betaFreq*beta2); #beta2 = meshpts.beta^2 
	beta2Dev <- beta2 - mean.beta2[gen] 
	var.beta2[gen] <- sum(betaFreq*beta2Dev^2); 
	
	if(gen%%10==0) cat(gen,"  ",mean.beta[gen],"   ",Recr/(Recr+n.surv),"\n")
}


return(list(mean.beta=mean.beta,var.beta=var.beta,mean.beta2=mean.beta2,var.beta2=var.beta2,
			prop.Recr=prop.Recr,nt=nt,betaDist=betaDist))
}

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

system.time(iterate_model(m.par.true,n.yrs,-20,init.beta.sd,meshpts.beta,meshpts.z))

sol.2D.minus20 <- iterate_model(m.par.true,n.yrs,-20,init.beta.sd,meshpts.beta,meshpts.z)

########################   solve beta(0) = -30 next
U.beta <- max(max.beta[6:10]) +  1 
L.beta <-  min(min.beta[6:10]) - 1
h.beta <- (U.beta - L.beta)/nBigMatrix.beta 
meshpts.beta <- L.beta + ((1:nBigMatrix.beta) - 1/2) * h.beta

U.z <- max(max.z[6:10]) + 0.5 
L.z <- min(min.z[6:10]) - 0.5 
h.z <- (U.z - L.z)/nBigMatrix.z
meshpts.z <- L.z + ((1:nBigMatrix.z) - 1/2) * h.z


sol.2D.minus30 <- iterate_model(m.par.true,n.yrs,-30,init.beta.sd,meshpts.beta,meshpts.z)



set_graph_pars("panel4"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Mean flowering intercept")
points(1:n.yrs,sol.2D.minus20$mean.beta,col="green",type="l")
points(1:n.yrs,sol.2D.minus30$mean.beta,col="blue",type="l")
abline(h=-24.92462,col="red")
add_panel_label("a")
matplot(,var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Variance flowering intercept")
points(1:n.yrs,sol.2D.minus20$var.beta,col="green",type="l")
points(1:n.yrs,sol.2D.minus30$var.beta,col="blue",type="l")
add_panel_label("b")

## plot beta distributions over time 
beta20 <- sol.2D.minus20$betaDist;
matplot(meshpts.beta,beta20[,c(1,100,300,500)], xlab="Beta", ylab="Relative Frequency",type="l"); 

beta30 <- sol.2D.minus30$betaDist;
matplot(meshpts.beta,beta30[,c(1,100,300,500)], xlab="Beta", ylab="Relative Frequency",type="l"); 


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


Approx_dynamics <- function(params,init.mean.flow.int,n.iter,delta,var.beta.t) {


	params["flow.int"]<- init.mean.flow.int
	beta.t.mean <- rep(NA,n.iter)
	
	beta.t.mean[1] <- init.mean.flow.int
	

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

#Calculate the mean intercept
	beta.t.mean[gen] <- beta.t.mean[gen-1] + d.lambda.d.beta * var.beta.t[gen-1]
		
	params["flow.int"] <- beta.t.mean[gen]
	
	if(gen%%10==0) cat(gen,"  ",beta.t.mean[gen],"  ",d.lambda.d.beta,"\n")
		
}

return(beta.t.mean)
}

nBigMatrix <- nBigMatrix.z; 

beta.minus20 <- Approx_dynamics(m.par.true,-20,n.yrs,0.1,sol.2D.minus20$var.beta)
beta.minus30 <- Approx_dynamics(m.par.true,-30,n.yrs,0.1,sol.2D.minus30$var.beta)

dev.new(); ### start plot 
set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",
	ylab="Mean beta",xlim=c(1,n.yrs)) 
points(1:n.yrs,sol.2D.minus20$mean.beta,col="blue",type="l",lty=1,lwd=4)
points(1:n.yrs,sol.2D.minus30$mean.beta,col="turquoise",type="l",lty=1,lwd=4)
points(1:n.yrs,beta.minus20,col="red",type="l",lty=2,lwd=2)
points(1:n.yrs,beta.minus30,col="red",type="l",lty=2,lwd=2)
legend("topright",legend=c("IPM -20","IPM -30","Approx"),col=c("blue","turquoise","red"),
    lty=1,bty="n",lwd=c(3,3,2));  

add_panel_label("a")

matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",
ylab="Variance flowering intercept",,xlim=c(0,n.yrs))
points(1:n.yrs,sol.2D.minus20$var.beta,col="blue",type="l",lwd=2)
points(1:n.yrs,sol.2D.minus30$var.beta,col="turquoise",type="l",lwd=2)
add_panel_label("b")


#dev.copy2eps(file="~/Repos/ipm_book/c10/figures/OenotheraIBMevolPreds.eps")





























































