## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Carlina IBM with individual variation in flowering intercept to illustrate evolutionary dynamics in a 
## stochastic environment
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls(all=TRUE))

library(doBy)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c9",sep="")); 

source("../utilities/Standard Graphical Pars.R");
source("Carlina Demog Funs.R");

#######################################################################################################
#  Simulate IBM with flowering intercept variation 
#######################################################################################################
# Set simulation parameters
init.pop.size <- 10000
init.mean.z <- 3
init.sd.z <- 2
n.yrs <-500
beta.off.sd <- 0.5
init.beta.sd <- 0.1
#Rec.mult multiples the actual number of recruits observed so we have a large population
Rec.mult <- 500

# this generates the simulation data, it takes a while so we use 
# a saved version - also note that the combination of demographic and 
# environmental stochasticity leads to err.... extinction at times.

mean.flow.ints <- matrix(NA,nrow=n.yrs,ncol=10)
var.flow.ints <- matrix(NA,nrow=n.yrs,ncol=10)
min.z <- rep(NA,10)
max.z <- rep(NA,10)
min.beta <- rep(NA,10)
max.beta <- rep(NA,10)

# for(reps in 1:10){
	# init.mean.flow.int <- ifelse(reps<6, -22, -5)
	# source("Carlina Simulate Evol IBM.R") 
    # mean.flow.ints[,reps] <- mean.flow.int
    # var.flow.ints[,reps] <- var.flow.int
    # min.z[reps] <- min.z.t
    # max.z[reps] <- max.z.t
    # min.beta[reps] <- min.flow.int
    # max.beta[reps] <- max.flow.int
    # cat(reps,"   ",init.mean.flow.int,"\n")
# }
# save(mean.flow.ints,var.flow.ints,min.z,max.z,min.beta,max.beta,file="SimDataStoc.Rdata")

load("SimDataStoc.Rdata")

set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Mean flowering intercept")
add_panel_label("a")
#Add ESS
abline(h= -14.34,col="red", lwd=2)
#Add simulation mean after the 1st 100 years
abline(h= mean(mean.flow.ints[100:500,],na.rm=TRUE),col="turquoise", lwd=2)

matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Variance flowering intercept")
add_panel_label("b")

#dev.copy2eps(file="~/Repos/ipm_book/c9/figures/CarlinaStocIBM.eps")

########################################################################################### 
# Simulate the 2-dimensional IPM with (size x intercept) structure 
###########################################################################################
m.par.sim <- m.par.true
#set mean flow.int <- 0 as each individual has it's own beta_0 and so this is the yearly deviation
m.par.sim["flow.int"] <- 0

nBigMatrix.z <- 60; 
nBigMatrix.beta <- 50; 

#######################  solve first for beta(0) = -5
U.beta <- max(max.beta[1:5]) +  1 
L.beta <-  min(min.beta[1:5]) - 1
h.beta <- (U.beta - L.beta)/nBigMatrix.beta 
meshpts.beta <- L.beta + ((1:nBigMatrix.beta) - 1/2) * h.beta

U.z <- max(max.z[1:5]) + 0.5 
L.z <- min(min.z[1:5]) - 0.5 
h.z <- (U.z - L.z)/nBigMatrix.z
meshpts.z <- L.z + ((1:nBigMatrix.z) - 1/2) * h.z

source("Carlina Simulate Evol IPM.R") # creates function to simulate the IPM 
sol.2D.minus20 <- iterate_model(m.par.true,n.yrs,-5,init.beta.sd,meshpts.beta,meshpts.z)

########################   solve beta(0) = -22 next
U.beta <- max(max.beta[6:10]) +  1 
L.beta <-  min(min.beta[6:10]) - 1
h.beta <- (U.beta - L.beta)/nBigMatrix.beta 
meshpts.beta <- L.beta + ((1:nBigMatrix.beta) - 1/2) * h.beta

U.z <- max(max.z[6:10]) + 0.5 
L.z <- min(min.z[6:10]) - 0.5 
h.z <- (U.z - L.z)/nBigMatrix.z
meshpts.z <- L.z + ((1:nBigMatrix.z) - 1/2) * h.z

sol.2D.minus30 <- iterate_model(m.par.true,n.yrs,-22,init.beta.sd,meshpts.beta,meshpts.z)

#### Plot the results 
set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Mean flowering intercept")
points(1:n.yrs,sol.2D.minus20$mean.beta,col="green",type="l")
points(1:n.yrs,sol.2D.minus30$mean.beta,col="blue",type="l")
abline(h=-14.34,col="red")
abline(h= mean(mean.flow.ints[100:500,],na.rm=TRUE),col="turquoise")

add_panel_label("a")

matplot(,var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab="Variance flowering intercept")
points(1:n.yrs,sol.2D.minus20$var.beta,col="green",type="l")
points(1:n.yrs,sol.2D.minus30$var.beta,col="blue",type="l")
add_panel_label("b")

#dev.copy2eps(file="~/Repos/ipm_book/c9/figures/CarlinaStocIBM.eps")

# ## plot beta distributions over time 
# beta20 <- sol.2D.minus20$betaDist;
# matplot(meshpts.beta,beta20[,c(1,100,300,450:500)], xlab="Beta", ylab="Relative Frequency",type="l"); 

# beta30 <- sol.2D.minus30$betaDist;
# matplot(meshpts.beta,beta30[,c(1,100,300,450:500)], xlab="Beta", ylab="Relative Frequency",type="l"); 

###################################################################################################
# Approximate dynamics, using Iwasa et al. 1991 Evolution
####################################################################################################
## simulate some trait yearly changes and store yearly parameters
init.pop.size <- 10000
Rec.mult <- 800
init.mean.z <- 3
init.sd.z <- 2
n.yrs <-500
beta.off.sd <- 0.5
init.beta.sd <- 3
init.mean.flow.int <- -14

nBigMatrix <- 100;

store.flow.int.sd <- m.par.sd.true["flow.int.sd"]
m.par.sd.true["flow.int.sd"] <- 0

source("Carlina Simulate Evol IBM.R") 

#m.par.sd.true["flow.int.sd"] <- store.flow.int.sd

# store.params.yr has the yearly parameters
# mean.flow.int the mean flowering intercepts
# var.flow.int the variance in flowering intercepts

U.z <- max.z.t + 0.5 
L.z <- min.z.t - 0.5 
delta <- 0.001

sel.t <- rep(NA,n.yrs-1)

#Calculate the actual changes 

IBM.sel <- mean.flow.int[2:n.yrs] - mean.flow.int[1:(n.yrs-1)]

for(gen in 1:(n.yrs-1)){
    params <- store.params.yr[gen,]
    IPM.kernel.parts <- mk_K(nBigMatrix, params, L.z, U.z)
    IPM.kernel <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
#  lambda <- Re(eigen(IPM.kernel,only.values = TRUE)$values[1])
#  wt <- Re(abs(eigen(IPM.kernel)$vectors[,1]))
#  vt <- Re(abs(eigen(t(IPM.kernel))$vectors[,1]))
#  sum(vt * IPM.kernel %*% wt)/sum(vt*wt)
#  lambda

#Calculate the derivative of log lambda with respect to beta
    params["flow.int"] <- params["flow.int"] + delta
    IPM.kernel.parts <- mk_K(nBigMatrix, params, L.z, U.z)
    IPM.kernel <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
    lambda.up <- log(Re(eigen(IPM.kernel,only.values = TRUE)$values[1]))
    params["flow.int"] <- params["flow.int"] - delta
    
	params["flow.int"] <- params["flow.int"] - delta
    IPM.kernel.parts <- mk_K(nBigMatrix, params, L.z, U.z)
	IPM.kernel <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
	lambda.down <- log(Re(eigen(IPM.kernel,only.values = TRUE)$values[1]))
	params["flow.int"] <- params["flow.int"] + delta

	d.log.lambda.d.beta <- (lambda.up - lambda.down)/(2*delta)

#Calculate the predicted change 
	sel.t[gen] <-  d.log.lambda.d.beta * var.flow.int[gen] 
	if(gen%%10==0) cat(gen,"  ",sel.t[gen],"  ",IBM.sel[gen],"\n")
}

e <- c(152:n.yrs)-1; 
par(mfrow=c(1,2),bty="l",pty="s")
plot(IBM.sel)
plot(sel.t[e],IBM.sel[e],xlab="Predicted change in mean(beta)",ylab="Observed change in mean(beta)")
IPMselfit <- lm(IBM.sel[e]~sel.t[e])
abline(IPMselfit); 
summary(IPMselfit)
cor.test(IBM.sel[e],sel.t[e]); # highly significant but rho is only 0.49.  
mean(IBM.sel[e]); mean(sel.t[e]);

#this doesn't work too well presumably as the population is a long way from the stable size distribution
#Let's compute the time dependent wt and vt and see if that works better...

### Get wt and Rt time series ###
	wt<-matrix(1/nBigMatrix, nrow=n.yrs+1, ncol=nBigMatrix);
	Rt<-rep(NA, n.yrs);
	for (i in 1:n.yrs) {
		params <- store.params.yr[i,]
		K<-mk_K(nBigMatrix, params, L.z, U.z);
		wt[i+1,]<-(K$P+params["p.r"]*K$F) %*% wt[i,]
		Rt[i]<-sum(wt[i+1,]);
		wt[i+1,]<-wt[i+1,]/Rt[i];
		if(i%%250==0) cat("wt and Rt ",i,"\n")

	}

	
### Get vt time series ###
	vt<-matrix(1/nBigMatrix, nrow=n.yrs+1, ncol=nBigMatrix);
	for (i in (n.yrs+1):2) {
		params <- store.params.yr[i-1,]
		K<-mk_K(nBigMatrix, params, L.z, U.z);
		vt[i-1,]<-vt[i,] %*% (K$P+params["p.r"]*K$F);
		vt[i-1,]<-vt[i-1,]/sum(vt[i-1,]);
		if(i%%250==0) cat("vt  ",i,"\n")

	}
	
for(gen in 1:(n.yrs-1)){

params <- store.params.yr[gen,]

IPM.kernel.parts <- mk_K(nBigMatrix, params, L.z, U.z)
IPM.kernel <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F

#Calculate the derivative of Kt with respect to beta

	params["flow.int"] <- params["flow.int"] + delta
	IPM.kernel.parts <- mk_K(nBigMatrix, params, L.z, U.z)
	IPM.kernel.up <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
	params["flow.int"] <- params["flow.int"] - delta
	
	params["flow.int"] <- params["flow.int"] - delta
	IPM.kernel.parts <- mk_K(nBigMatrix, params, L.z, U.z)
	IPM.kernel.down <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
	params["flow.int"] <- params["flow.int"] + delta

	d.Kt.d.beta <- (IPM.kernel.up - IPM.kernel.down)/(2*delta)

#Selection in year t
	sel.t[gen] <-  sum(vt[gen+1,] * (d.Kt.d.beta %*% wt[gen,]))/sum(vt[gen+1,] * (IPM.kernel %*% wt[gen,]))
	
	sel.t[gen] <- sel.t[gen] * var.flow.int[gen]
	
	
	if(gen%%10==0) cat(gen,"  ",sel.t[gen],"  ",IBM.sel[gen],"\n")
		
}

e <- c(152:n.yrs)-1; 
par(mfrow=c(1,2),bty="l",pty="s")
plot(IBM.sel)
plot(sel.t[e],IBM.sel[e],xlab="Predicted change in mean(beta)",ylab="Observed change in mean(beta)")
IPMselfit <- lm(IBM.sel[e]~sel.t[e])
abline(IPMselfit); 
summary(IPMselfit)
cor.test(IBM.sel[e],sel.t[e]); # highly significant but rho is only 0.49.  
mean(IBM.sel[e]); mean(sel.t[e]);

m.par.sd.true["flow.int.sd"] <- store.flow.int.sd


######################################################################
# End section on approximating stochastic Carlina evol dynamics(?) Yes this is old code, sorry Steve
#####################################################################

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

#dev.copy2eps(file="~/Repos/ipm_book/c9/figures/OenotheraIBMevolPreds.eps")





























































