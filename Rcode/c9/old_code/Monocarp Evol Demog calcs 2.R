## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IBM with individual variation in flowering intercept to illustrate evolutionary dynamics
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
setwd("~/Repos/ipm_book/Rcode/c2")
source("Monocarp Demog Funs.R");

# 
# s_z <- function(z,m.par) {
# 	 p <- 0.36 + 0.17 *z
# 	 p[p>0.7] <- 0.7
# 	 p[p<0.0] <- 0.0
# 	 return(p)
#  }

source("../utilities/Standard Graphical Pars.R");

# Set simulation parameters
setwd("~/Repos/ipm_book/Rcode/c10")

init.pop.size <- 10000
Recr <- 5000
n.yrs <-2000
beta.off.sd <- 0.3
init.beta.sd <- 0.1

# # mean.flow.ints <- matrix(NA,nrow=n.yrs,ncol=10)
# var.flow.ints <- matrix(NA,nrow=n.yrs,ncol=10)

# for(reps in 1:10){
	# init.mean.flow.int <- ifelse(reps<6, -20, -30)
	# source("Monocarp Simulate Evol IBM.R") 
    # mean.flow.ints[,reps] <- mean.flow.int
    # var.flow.ints[,reps] <- var.flow.int
    # cat(reps,"   ",init.mean.flow.int,"\n")
# }

# save(mean.flow.ints,var.flow.ints,file="SimData.Rdata")

load("SimData.Rdata")

set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Mean flowering intercept")
add_panel_label("a")
matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Variance flowering intercept")
add_panel_label("b")

dev.copy2eps(file="~/Repos/ipm_book/c10/figures/OenotheraIBMevol.eps")

#Implement a 2-dimensional IPM 

iterate_model <- function(params,n.iter,init.z.pdf,init.beta.pdf) {
	
	#nt each column is the density of z for a given flowering intercept
    nt <- matrix(NA,nrow=nBigMatrix.z,ncol=nBigMatrix.beta)
    nt <- init.pop.size*matrix(outer(z.s,flow.int.s),nrow=nBigMatrix.z)
	
	#array of P matricies for each beta
    P.beta <- array(NA,c(nBigMatrix.beta,nBigMatrix.z,nBigMatrix.z))

    Seeds.beta <- array(NA,c(nBigMatrix.beta,nBigMatrix.z))

    Seeds_z <- function ( z, m.par) {
            return( p_bz(z, m.par) * b_z(z, m.par) )
     }

	for(i in 1:nBigMatrix.beta){
		params["flow.int"] <- meshpts.beta[i]
		P.beta[i,1:nBigMatrix.z,1:nBigMatrix.z] <- h.z * (outer(meshpts.z, meshpts.z, P_z1z, m.par = 			params))
		Seeds.beta[i,1:nBigMatrix.z] <- Seeds_z(meshpts.z, params)
	}

	R.beta <- matrix(NA,ncol=nBigMatrix.beta,nrow=nBigMatrix.beta)

	for(i in 1:nBigMatrix.beta){
 		 R.beta[,i] <- dnorm(meshpts.beta,mean=meshpts.beta[i],sd=beta.off.sd)
  	}

#seeds.from.betai <- rep(NA,nBigMatrix.beta)

	off.pdf <- c_0z1(meshpts.z,m.par.true)

	mean.beta <- rep(NA,n.iter)
	var.beta <- rep(NA,n.iter)

	index <- 1:nBigMatrix.beta

#Iterate model

	for(gen in 1:n.iter){

# 
# for(i in 1:nBigMatrix){
# 
# 	seeds.from.betai[i] <- h.z*sum(Seeds.beta[i,] %*% nt[,i])
# 	
# }

	seeds.from.betai <- sapply(index,function(i) h.z*sum(Seeds.beta[i,] %*% nt[,i]))

#total.seeds <- sum(seeds.from.betai)

	seeds.with.betai <- h.beta * (R.beta %*% seeds.from.betai)

#total.seeds

#sum(seeds.with.betai)

	f.recruits.with.betai <- seeds.with.betai / sum(seeds.with.betai)
# 
# for(i in 1:nBigMatrix){
# 
# 	nt1[,i] <- P.beta[i,,] %*% nt[,i] + Recr * off.pdf * f.recruits.with.betai[i]
# 	
# }

	nt1 <- sapply(index,function(i) P.beta[i,,] %*% nt[,i] + Recr * off.pdf * f.recruits.with.betai[i])

	nt <- nt1

	mean.beta[gen] <- sum(apply(nt1,2,sum)*meshpts.beta/sum(nt1))
	var.beta[gen] <- sum(apply(nt1,2,sum)*meshpts.beta*meshpts.beta/sum(nt1)) - 	  mean.beta[gen]*mean.beta[gen]


cat(gen,"  ",mean.beta[gen],"\n")


}

return(list(mean.beta=mean.beta,var.beta=var.beta))
	
}


nBigMatrix.z <- 100; 
nBigMatrix.beta <- 150; 

U.beta <- -15 * 0.9
L.beta <- -35 * 1.1

h.beta <- (U.beta - L.beta)/nBigMatrix.beta 
meshpts.beta <- L.beta + ((1:nBigMatrix.beta) - 1/2) * h.beta

U.z <- 6.2 * 1.1 
L.z <- -4 * 1.1 

h.z <- (U.z - L.z)/nBigMatrix.z
meshpts.z <- L.z + ((1:nBigMatrix.z) - 1/2) * h.z

#nt each column is the density of z for a given flowering intercept

z.s <- dnorm(meshpts.z, mean = 2.9, sd = m.par.true["rcsz.sd"])

flow.int.s <- dnorm(meshpts.beta, mean = -20, sd = init.beta.sd)

sol.2D.minus20 <- iterate_model(m.par.true,n.yrs,z.s,flow.int.s)

flow.int.s <- dnorm(meshpts.beta, mean = -30, sd = init.beta.sd)

sol.2D.minus30 <- iterate_model(m.par.true,n.yrs,z.s,flow.int.s)

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

#dev.copy2eps(file="~/Repos/ipm_book/c10/figures/OenotheraIBMevol.eps")




# # 




# par(mfrow=c(2,2),pty="s",bty="l")

# plot(1:n.yrs,mean.flow.int.minus20,type="l",ylim=c(-30,-20))
# points(1:n.iter,mean.beta,type="l",col="green")
# abline(h=ESS.fl.int,col="red") #adding ESS from calculation below


# # plot(1:n.yrs,var.flow.int.minus20,type="l",)
# points(1:n.iter,var.beta,type="l",col="green")

# plot(meshpts.beta,apply(nt1,2,sum)/sum(nt1),type="l")

# params <- m.par.true

# params["flow.int"]<- ESS.fl.int

# plot(meshpts.z,p_bz(meshpts.z,params),type="l",col="red")

# params["flow.int"]<- -24.24019

# points(meshpts.z,p_bz(meshpts.z,params),type="l",col="green")

# params["flow.int"]<- -30
# 
# points(meshpts.z,p_bz(meshpts.z,params),type="l",col="green")
# 
# params["flow.int"]<- -20
# 
# points(meshpts.z,p_bz(meshpts.z,params),type="l",col="blue")

###################################################################################################
#Approximate dynamics using Iwasa et al. 1991 Evolution
####################################################################################################

R0_calc <- function (params) {
	
	#params["flow.int"] <- fl.int
	#params["flow.z"] <- fl.int
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

	params["p.r"] <- 1
	equ.p.r <- 1/R0_calc(params)
	params["p.r"] <- equ.p.r
	
	params["flow.int"] <- params["flow.int"] + delta
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	lambda.up <- Re(eigen(IPM.kernel$K,only.values = TRUE)$values[1])
	params["flow.int"] <- params["flow.int"] - delta
	
	params["flow.int"] <- params["flow.int"] - delta
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	lambda.down <- Re(eigen(IPM.kernel$K,only.values = TRUE)$values[1])
	params["flow.int"] <- params["flow.int"] + delta

	d.lambda.d.beta <- (lambda.up - lambda.down)/(2*delta)
	
	beta.t.mean[gen] <- beta.t.mean[gen-1] + d.lambda.d.beta * var.beta.t[gen-1]
		
	params["flow.int"] <- beta.t.mean[gen]
	
	cat(gen,"  ",beta.t.mean[gen],"  ",d.lambda.d.beta,"\n")
		
}

return(beta.t.mean)
}

nBigMatrix <- nBigMatrix.z; 

beta.minus20 <- Approx_dynamics(m.par.true,-20,n.yrs,0.1,sol.2D.minus20$var.beta)
beta.minus30 <- Approx_dynamics(m.par.true,-30,n.yrs,0.1,sol.2D.minus20$var.beta)


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







nBigMatrix <- 100; 

U.z <- 1.1 * max.z.minus20
L.z <- 1.1 * min.z.minus20

# IPM.true <- mk_K(nBigMatrix, m.par.true, -3.5, 5.0)
# 
# meshpts <- IPM.true$meshpts
# h <- diff(meshpts)[1]
# 
# lambda0 <- Re(eigen(IPM.true$K,only.values = TRUE)$values[1])

params <- m.par.true
n.ints <- 200

R0s <- rep(NA,n.ints)
fl.ints <- seq(-24,-26,length=n.ints)
#fl.ints <- seq(2,10,length=100)

for(i in 1:n.ints) {
	R0s[i] <- R0_calc(fl.ints[i])
	}
	
ESS.fl.int <- fl.ints[which(R0s==max(R0s))]

ESS.fl.int

plot(fl.ints, R0s,type="l")
abline(v=m.par.true["flow.int"])
abline(v=ESS.fl.int,col="red")


par(mfrow=c(2,2),pty="s",bty="l")

t.plot <- 1:100

plot(t.plot,mean.flow.int.minus20[t.plot],type="l",ylim=c(-30,-20))
points(t.plot,mean.beta[t.plot],type="l",col="green")
points(t.plot,beta.t.mean[t.plot],type="l",col="blue")
abline(h=ESS.fl.int,col="red") #adding ESS from calculation below


plot(t.plot,var.flow.int.minus20[t.plot],type="l",ylim=c(0,6))
points(t.plot,var.beta[t.plot],type="l",col="green")
points(t.plot,beta.t.var[t.plot],type="l",col="blue")

plot(meshpts.beta,apply(nt1,2,sum)/sum(nt1),type="l")

params <- m.par.true

params["flow.int"]<- ESS.fl.int

plot(meshpts.z,p_bz(meshpts.z,params),type="l",col="red")

params["flow.int"]<- mean.beta[n.iter]

points(meshpts.z,p_bz(meshpts.z,params),type="l",col="green")

params["flow.int"]<- beta.t[n.iter]

points(meshpts.z,p_bz(meshpts.z,params),type="l",col="blue")

##########################################################################################
#Why does it work?
##########################################################################################

nt1.beta <- apply(nt1,2,sum)

nt1.z <- apply(nt1,1,sum)/sum(nt1)

lambdas <- rep(NA,nBigMatrix.beta)

for(i in 1:nBigMatrix.beta){
 params["flow.int"] <- meshpts.beta[i]
 IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
 lambdas[i] <- Re(eigen(IPM.kernel$K,only.values = TRUE)$values[1])
}

plot(meshpts.beta,f.recruits.with.betai /nt1.beta)
plot(meshpts.beta,lambdas )




























































