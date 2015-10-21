## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Carlina IBM with individual variation in flowering intercept to test the "Stochastic
## Iwasa" approximation.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls(all=TRUE))

# set.seed(53241986)

library(smatr)

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c9",sep="")); 

source("../utilities/Standard Graphical Pars.R");
source("Carlina Demog Funs.R");
source("Carlina Simulate Evol IPM.R") # creates function to simulate the IPM 

# Set simulation parameters
init.pop.size <- 10000
init.mean.z <- 3
init.sd.z <- 2
n.yrs <-500
beta.off.sd <- 0.25
init.beta.sd <- 1
#Rec.mult multiples the actual number of recruits observed so we have a large population
Rec.mult <- 500

########################################################################################### 
# Simulate the 2-dimensional IPM with (size x intercept) structure 
###########################################################################################
m.par.sim <- m.par.true
#set mean flow.int <- 0 as each individual has it's own beta_0 and so this is the yearly deviation
m.par.sim["flow.int"] <- 0

nBigMatrix.z <- 80; 
nBigMatrix.beta <- 60; 

#######################  solve for beta(0) = -10
U.beta <- 5 
L.beta <-  -25 
h.beta <- (U.beta - L.beta)/nBigMatrix.beta 
meshpts.beta <- L.beta + ((1:nBigMatrix.beta) - 1/2) * h.beta

U.z <- 9 
L.z <- -3 
h.z <- (U.z - L.z)/nBigMatrix.z
meshpts.z <- L.z + ((1:nBigMatrix.z) - 1/2) * h.z

IPM.sol <- iterate_model(m.par.true,n.yrs,-10,init.beta.sd,meshpts.beta,meshpts.z)

#### Plot the results 
dev.new(width=8,height=4); 
set_graph_pars("panel2"); 
plot(1:n.yrs,IPM.sol$mean.beta,col="green",type="l",xlab="Time",ylab="Mean flowering intercept")
abline(h=-14.34,col="red")
add_panel_label("a")

plot(1:n.yrs,sqrt(IPM.sol$var.beta),col="green",type="l",xlab="Time",ylab="Std Dev flowering intercept")
add_panel_label("b")

###################################################################################################
# Approximate dynamics, using Iwasa et al. 1991 Evolution
####################################################################################################
IBM.sel <- IPM.sol$mean.beta[2:n.yrs] - IPM.sol$mean.beta[1:(n.yrs-1)]
params.yr <- IPM.sol$params.yr; store.params.yr <- params.yr; 
delta <- .Machine$double.eps^(1/3); sel.t <- numeric(n.yrs); 

nBigMatrix<- nBigMatrix.z; 

for(gen in 1:(n.yrs-1)){
    params <- params.yr[gen,]

#Calculate the derivative of log lambda with respect to beta
	#int.gen <- m.par.true["flow.int"] + params["flow.int"]
	int.gen <- IPM.sol$mean.beta[gen] + params["flow.int"]
    params["flow.int"] <- int.gen + delta
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
    IPM.kernel <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
    lambda.up <- log(Re(eigen(IPM.kernel,only.values = TRUE)$values[1]))
    
	params["flow.int"] <- int.gen - delta
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
	IPM.kernel <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
	lambda.down <- log(Re(eigen(IPM.kernel,only.values = TRUE)$values[1]))

	d.log.lambda.d.beta <- (lambda.up - lambda.down)/(2*delta)

#Calculate the predicted change 
	sel.t[gen] <-  d.log.lambda.d.beta * IPM.sol$var.beta[gen] 
	if(gen%%10==0) cat(gen,"  ",sel.t[gen],"  ",IBM.sel[gen],"\n")
}

#e <- c(102:n.yrs)-1; 
e <- c(102:(n.yrs-100))-1;

#e <- c(10:200); 

dev.new(width=8,height=8); 
par(mfrow=c(2,2),bty="l",pty="s"); 
#plot(IBM.sel)
plot(sel.t[e],IBM.sel[e],xlab="Predicted change in mean(beta)",ylab="Observed change in mean(beta)")
summary(lm(IBM.sel[e]~sel.t[e]))
IPMselfit <- sma(IBM.sel[e]~sel.t[e])
out <- cor.test(IBM.sel[e],sel.t[e]);   
out; title(main=out$estimate); 
abline(IPMselfit); abline(0,1,lty=2,col="red"); 
summary(IPMselfit)
mean(IBM.sel[e]); mean(sel.t[e]);

#############################################################################################################
#this doesn't work too well presumably as the population is a long way from the stable size distribution
#Let's compute the time dependent wt and vt and see if that works better...
#############################################################################################################

### Get wt and Rt time series ###
	wt<-matrix(1/nBigMatrix.z, nrow=n.yrs+1, ncol=nBigMatrix.z);
	Rt<-rep(NA, n.yrs);
	for (i in 1:n.yrs) {
        params <- params.yr[i,]
        int.gen <- IPM.sol$mean.beta[i] + params["flow.int"]
        params["flow.int"] <- int.gen 
    	K<-mk_K(nBigMatrix.z, params, L.z, U.z);
		wt[i+1,]<-(K$P+params["p.r"]*K$F) %*% wt[i,]
		Rt[i]<-sum(wt[i+1,]);
		wt[i+1,]<-wt[i+1,]/Rt[i];
		if(i%%250==0) cat("wt and Rt ",i,"\n")

	}

### Get vt time series ###
	vt<-matrix(1/nBigMatrix.z, nrow=n.yrs+1, ncol=nBigMatrix.z);
	for (i in (n.yrs):2) {
        params <- params.yr[i-1,]
        int.gen <- IPM.sol$mean.beta[i-1]+ params["flow.int"]
        params["flow.int"] <- int.gen 
    
		K<-mk_K(nBigMatrix.z, params, L.z, U.z);
		vt[i-1,]<-vt[i,] %*% (K$P+params["p.r"]*K$F);
		vt[i-1,]<-vt[i-1,]/sum(vt[i-1,]);
		if(i%%250==0) cat("vt  ",i,"\n")

	}
	
for(gen in 1:(n.yrs-1)){
    params <- params.yr[gen,]

#Calculate the derivative of log lambda with respect to beta
	int.gen <- IPM.sol$mean.beta[gen] + params["flow.int"]
    
    params["flow.int"] <- int.gen
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
    IPM.kernel <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
    
    params["flow.int"] <- int.gen + delta
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
    IPM.kernel.up <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
 
    params["flow.int"] <- int.gen - delta
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
	IPM.kernel.down <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F

	d.Kt.d.beta <- (IPM.kernel.up - IPM.kernel.down)/(2*delta)

#Selection in year t
	sel.t[gen] <-  sum(vt[gen+1,] * (d.Kt.d.beta %*% wt[gen,]))/sum(vt[gen+1,] * (IPM.kernel %*% wt[gen,]))
	sel.t[gen] <- sel.t[gen] * IPM.sol$var.beta[gen] 
	if(gen%%10==0) cat(gen,"  ",sel.t[gen],"  ",IBM.sel[gen],"\n")
		
}

#e <- c(102:(n.yrs-100))-1; 
#plot(IBM.sel); 
plot(sel.t[e],IBM.sel[e],xlab="Predicted change in mean(beta)",ylab="Observed change in mean(beta)")
summary(lm(IBM.sel[e]~sel.t[e]))
IPMselfit <- sma(IBM.sel[e]~sel.t[e])
abline(IPMselfit);abline(0,1,lty=2,col="red"); 
out <- cor.test(IBM.sel[e],sel.t[e]); # highly significant but rho is only 0.49.  
out; title(main=out$estimate); 

summary(IPMselfit)
mean(IBM.sel[e]); mean(sel.t[e]);


#############################################################################################################
# Third version: just time-dependent w, constant v  
#############################################################################################################

vt[,]=1; Kbar <- 0*IPM.kernel; 

for(gen in 1:(n.yrs-1)){
    params <- params.yr[gen,]

#Calculate the derivative of log lambda with respect to beta
	int.gen <- IPM.sol$mean.beta[gen] + params["flow.int"]
    
    params["flow.int"] <- int.gen
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
    IPM.kernel <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
	Kbar=Kbar+IPM.kernel; 
    
    params["flow.int"] <- int.gen + delta
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
    IPM.kernel.up <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
 
    params["flow.int"] <- int.gen - delta
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
	IPM.kernel.down <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F

	d.Kt.d.beta <- (IPM.kernel.up - IPM.kernel.down)/(2*delta)

#Selection in year t
	sel.t[gen] <-  sum(vt[gen+1,] * (d.Kt.d.beta %*% wt[gen,]))/sum(vt[gen+1,] * (IPM.kernel %*% wt[gen,]))
	sel.t[gen] <- sel.t[gen] * IPM.sol$var.beta[gen] 
	if(gen%%10==0) cat(gen,"  ",sel.t[gen],"  ",IBM.sel[gen],"\n")
		
}
Kbar = Kbar/(n.yrs-1); 

#dev.new(); par(mfrow=c(2,2)); 
#e <- c(102:(n.yrs-100))-1; 
#plot(IBM.sel); 
plot(sel.t[e],IBM.sel[e],xlab="Predicted change in mean(beta)",ylab="Observed change in mean(beta)")
summary(lm(IBM.sel[e]~sel.t[e]))
IPMselfit <- sma(IBM.sel[e]~sel.t[e])
abline(IPMselfit); abline(0,1,lty=2,col="red"); 
out <- cor.test(IBM.sel[e],sel.t[e]); # highly significant but rho is only 0.49.  
out; title(main=out$estimate); 

summary(IPMselfit)
mean(IBM.sel[e]); mean(sel.t[e]);

#############################################################################################################
# Fourth version: constant w (from mean kernel), v==1  
#############################################################################################################

vt[,]=1; w=Re(eigen(Kbar)$vectors[,1]); w=w/sum(w); 

for(gen in 1:(n.yrs-1)){
    params <- params.yr[gen,]

#Calculate the derivative of log lambda with respect to beta
	int.gen <- IPM.sol$mean.beta[gen] + params["flow.int"]
    
    params["flow.int"] <- int.gen
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
    IPM.kernel <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
	Kbar=Kbar+IPM.kernel; 
    
    params["flow.int"] <- int.gen + delta
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
    IPM.kernel.up <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F
 
    params["flow.int"] <- int.gen - delta
    IPM.kernel.parts <- mk_K(nBigMatrix.z, params, L.z, U.z)
	IPM.kernel.down <- IPM.kernel.parts$P + params["p.r"]*IPM.kernel.parts$F

	d.Kt.d.beta <- (IPM.kernel.up - IPM.kernel.down)/(2*delta)

#Selection in year t
	sel.t[gen] <-  sum(vt[gen+1,] * (d.Kt.d.beta %*% w))/sum(vt[gen+1,] * (IPM.kernel %*% w))
	sel.t[gen] <- sel.t[gen] * IPM.sol$var.beta[gen] 
	if(gen%%10==0) cat(gen,"  ",sel.t[gen],"  ",IBM.sel[gen],"\n")
		
}
Kbar = Kbar/(n.yrs-1); 

#plot(IBM.sel); 
plot(sel.t[e],IBM.sel[e],xlab="Predicted change in mean(beta)",ylab="Observed change in mean(beta)")
summary(lm(IBM.sel[e]~sel.t[e]))
IPMselfit <- sma(IBM.sel[e]~sel.t[e])
abline(IPMselfit); abline(0,1,lty=2,col="red"); 
out <- cor.test(IBM.sel[e],sel.t[e]); # highly significant but rho is only 0.49.  
out; title(main=out$estimate); 

summary(IPMselfit)
mean(IBM.sel[e]); mean(sel.t[e]);




