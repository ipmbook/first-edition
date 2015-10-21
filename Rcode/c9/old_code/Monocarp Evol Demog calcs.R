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
init.pop.size <- 250
n.yrs <-5000
init.mean.flow.int <- -20


setwd("~/Repos/ipm_book/Rcode/c10")
source("Monocarp Simulate Evol IBM.R") 

mean.flow.int.minus20 <- mean.flow.int
mean.fl.z.t.minus20 <- mean.fl.z.t 
var.flow.int.minus20 <- var.flow.int
min.flow.int.minus20 <- min.flow.int
max.flow.int.minus20 <- max.flow.int

init.mean.flow.int <- -30

source("Monocarp Simulate Evol IBM.R") 

mean.flow.int.minus30 <- mean.flow.int
mean.fl.z.t.minus30 <- mean.fl.z.t 
var.flow.int.minus30 <- var.flow.int
min.flow.int.minus30 <- min.flow.int
max.flow.int.minus30 <- max.flow.int

par(mfrow=c(1,2),pty="s",bty="l")

plot(1:n.yrs,mean.flow.int.minus20,type="l",ylim=c(-30,-20))
points(1:n.yrs,mean.flow.int.minus30,type="l",col="grey")
abline(h=-24.54545,col="red") #adding ESS from calculation below

plot(1:n.yrs,mean.fl.z.t.minus20,type="l",ylim=c(20,100))
points(1:n.yrs,mean.fl.z.t.minus30,type="l",col="grey")
abline(h=-24.54545,col="red") #adding ESS from calculation below

# 
# plot(1:n.yrs,var.flow.int.minus20,type="l")
# points(1:n.rs,var.flow.int.minus30,type="l",col="grey")
# #abline(h=-24.54545,col="red") #adding ESS from calculation below

#Do some ESS calculations

nBigMatrix <- 250; 

# IPM.true <- mk_K(nBigMatrix, m.par.true, -3.5, 5.0)
# 
# meshpts <- IPM.true$meshpts
# h <- diff(meshpts)[1]
# 
# lambda0 <- Re(eigen(IPM.true$K,only.values = TRUE)$values[1])

params <- m.par.true

R0_calc <- function (fl.int) {
	
	params["flow.int"] <- fl.int
	#params["flow.z"] <- fl.int
	IPM.true <- mk_K(nBigMatrix, params, -3.5, 5.0)
	# to keep close to the formulae in the text next we define the F and P iteration matricies
    P <- IPM.true$P;  F <- IPM.true$F;

    # Fundamental operator 
     N <- solve(diag(nBigMatrix)-P); 

    # Compute R0 as dominant eigenvalue of FN
    R <- F %*% N
    R0 <- abs(eigen(R)$values[1])
    return(R0)
}

R0s <- rep(NA,100)
fl.ints <- seq(-5,-50,length=100)
#fl.ints <- seq(2,10,length=100)

for(i in 1:100) {
	R0s[i] <- R0_calc(fl.ints[i])
	}
	
ESS.fl.int <- fl.ints[which(R0s==max(R0s))]

ESS.fl.int

plot(fl.ints, R0s,type="l")
abline(v=m.par.true["flow.int"])
abline(v=ESS.fl.int,col="red")

#Implement a 2-dimensional IPM 

init.mean.flow.int <- -20
init.pop.size <- 250

nBigMatrix <- 100; 

U.beta <- max.flow.int.minus20 * 0.9
L.beta <- min.flow.int.minus20 * 1.1

h.beta <- (U.beta - L.beta)/nBigMatrix 
meshpts.beta <- L.beta + ((1:nBigMatrix) - 1/2) * h.beta

U.z <- 5.0
L.z <- -3.5

h.z <- (U.z - L.z)/nBigMatrix
meshpts.z <- L.z + ((1:nBigMatrix) - 1/2) * h.z

#nt each column is the density of z for a given flowering intercept

z.s <- dnorm(meshpts.z, mean = 2.9, sd = m.par.true["rcsz.sd"])
flow.int.s <- dnorm(meshpts.beta, mean = init.mean.flow.int, sd = 0.1)

sum(z.s * h.z)

sum(flow.int.s * h.beta)

#nt each column is the density of z for a given flowering intercept

nt <- matrix(NA,nrow=nBigMatrix,ncol=nBigMatrix)
# 
# for(i in 1:nBigMatrix){
#  for(j in 1:nBigMatrix){
#  nt[i,j] <- flow.int.s[j] * z.s[i]
#  }
#  }
#  

nt <- init.pop.size*matrix(outer(z.s,flow.int.s),nrow=nBigMatrix)

sum(apply(nt,2,sum)*meshpts.beta/sum(nt))

nt1 <- nt

#array of P matricies for each beta
P.beta <- array(NA,c(nBigMatrix,nBigMatrix,nBigMatrix))

Seeds.beta <- array(NA,c(nBigMatrix,nBigMatrix))

Seeds_z <- function ( z, m.par) {

    return( p_bz(z, m.par) * b_z(z, m.par) )

}

params <- m.par.true

for(i in 1:nBigMatrix){
	params["flow.int"] <- meshpts.beta[i]
	P.beta[i,1:nBigMatrix,1:nBigMatrix] <- h.z * (outer(meshpts.z, meshpts.z, P_z1z, m.par = params))
	Seeds.beta[i,1:nBigMatrix] <- Seeds_z(meshpts.z, params)
}

R.beta <- nt


for(i in 1:nBigMatrix){
 R.beta[,i] <- dnorm(meshpts.beta,mean=meshpts.beta[i],sd=0.1)
 }

seeds.from.betai <- rep(NA,nBigMatrix)

off.pdf <- c_0z1(meshpts.z,m.par.true)

n.iter <- 5000

mean.beta <- rep(NA,n.iter)
var.beta <- rep(NA,n.iter)

index <- 1:nBigMatrix

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
var.beta[gen] <- sum(apply(nt1,2,sum)*meshpts.beta*meshpts.beta/sum(nt1)) - mean.beta[gen]*mean.beta[gen]


cat(gen,"  ",mean.beta[gen],"\n")


}

par(mfrow=c(1,2),pty="s",bty="l")

plot(1:n.yrs,mean.flow.int.minus20,type="l",ylim=c(-30,-20))
points(1:n.iter,mean.beta,type="l",col="green")
abline(h=-24.54545,col="red") #adding ESS from calculation below


plot(1:n.yrs,var.flow.int.minus20,type="l",)
points(1:n.iter,var.beta,type="l",col="green")























