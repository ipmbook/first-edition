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
init.beta.sd <- 0.7
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

for(reps in 1:10){
	init.mean.flow.int <- ifelse(reps<6, -22, -5)
	source("Carlina Simulate Evol IBM.R") 
    mean.flow.ints[,reps] <- mean.flow.int
    var.flow.ints[,reps] <- var.flow.int
    min.z[reps] <- min.z.t
    max.z[reps] <- max.z.t
    min.beta[reps] <- min.flow.int
    max.beta[reps] <- max.flow.int
    cat(reps,"   ",init.mean.flow.int,"\n")
}
# save(mean.flow.ints,var.flow.ints,min.z,max.z,min.beta,max.beta,file="SimDataStoc.Rdata")

load("SimDataStoc.Rdata")
dev.new(height=4,width=8);
set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab=expression(bar(beta)[0]))
add_panel_label("a")
#Add ESS
abline(h= -14.34,col="red", lwd=2)
#Add simulation mean after the 1st 100 years
abline(h= mean(mean.flow.ints[100:500,],na.rm=TRUE),col="turquoise", lwd=2)

matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab=expression("Var("* beta[0] *")"))
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

dev.new(height=4,width=8);
set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab=expression(bar(beta)[0]))
points(1:n.yrs,sol.2D.minus20$mean.beta,col="green",type="l")
points(1:n.yrs,sol.2D.minus30$mean.beta,col="blue",type="l")
abline(h=-14.34,col="red")
abline(h= mean(mean.flow.ints[100:500,],na.rm=TRUE),col="turquoise")

add_panel_label("a")

matplot(,var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",ylab=expression("Var("* beta[0] *")"))
points(1:n.yrs,sol.2D.minus20$var.beta,col="green",type="l")
points(1:n.yrs,sol.2D.minus30$var.beta,col="blue",type="l")
add_panel_label("b")

#dev.copy2eps(file="~/Repos/ipm_book/c9/figures/CarlinaStocIBM.eps")





























