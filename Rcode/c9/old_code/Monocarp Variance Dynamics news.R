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

########### Run the IBM simulation 
if(FALSE) {
  mean.flow.ints <- matrix(NA,nrow=n.yrs,ncol=10)
  var.flow.ints <- matrix(NA,nrow=n.yrs,ncol=10)
  min.z <- rep(NA,10)
  max.z <- rep(NA,10)
  min.beta <- rep(NA,10)
  max.beta <- rep(NA,10)

  for(reps in 1:10){
     init.mean.flow.int <- ifelse(reps<6, -20, -30)
	 source("../c9/Monocarp Simulate Evol IBM.R") 
     mean.flow.ints[,reps] <- mean.flow.int
     var.flow.ints[,reps] <- var.flow.int
     min.z[reps] <- min.z.t
     max.z[reps] <- max.z.t
     min.beta[reps] <- min.flow.int
     max.beta[reps] <- max.flow.int
     cat(reps,"   ",init.mean.flow.int,"\n")
  }

save(mean.flow.ints,var.flow.ints,min.z,max.z,min.beta,max.beta,
n.yrs,beta.off.sd,init.beta.sd,file="../c9/SimDataSmallVar.Rdata")
}


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

#nt each column is the density of z for a given flowering intercept
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

#nt each column is the density of z for a given flowering intercept
sol.2D.minus30 <- iterate_model(m.par.true,n.yrs,-30,init.beta.sd,meshpts.beta,meshpts.z)

## plot beta distributions over time 
beta20 <- sol.2D.minus20$betaDist;
matplot(meshpts.beta^2,beta20[,c(1,100,300,500)], xlab="Beta", ylab="Relative Frequency",type="l"); 

beta30 <- sol.2D.minus30$betaDist;
matplot(meshpts.beta^2,beta30[,c(1,100,300,500)], xlab="Beta", ylab="Relative Frequency",type="l"); 
######################### END FIRST PLOT 

#################################################################
#  SECOND PLOT 
#  Compare IBM to IPM; compare IPM to Gaussian moment-closure 
##################################################################
dev.new(); 
set_graph_pars("panel4"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",
ylab="Mean flowering intercept",xlim=c(0,n.yrs))
points(1:n.yrs,sol.2D.minus20$mean.beta,col="blue",type="l",lwd=2)
points(1:n.yrs,sol.2D.minus30$mean.beta,col="turquoise",type="l",lwd=2)
abline(h=-24.92462,col="red")
legend("topright",legend=c("IBM -20","IBM -30", "IPM -20", "IPM -30"),
	col=c("black","grey","blue","turquoise"),lty=1,lwd=2,bty="n");  

add_panel_label("a")

matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",
ylab="Variance flowering intercept",,xlim=c(0,n.yrs))
points(1:n.yrs,sol.2D.minus20$var.beta,col="blue",type="l",lwd=2)
points(1:n.yrs,sol.2D.minus30$var.beta,col="turquoise",type="l",lwd=2)
add_panel_label("b")

plot(1:n.yrs,sol.2D.minus20$mean.beta2,col="blue",type="l",ylim=c(0,1000),lwd=2,ylab="Mean of beta^2")
points(1:n.yrs,sol.2D.minus30$mean.beta2,col="turquoise",type="l",lwd=2)
add_panel_label("c")

#calculate variance of S assuming Gaussian distribution
minus20.var <- sol.2D.minus20$mean.beta^4 + 6*(sol.2D.minus20$mean.beta^2)*sol.2D.minus20$var.beta +
                 3*sol.2D.minus20$var.beta^2 - sol.2D.minus20$mean.beta2^2

minus30.var <- sol.2D.minus30$mean.beta^4 + 6*(sol.2D.minus30$mean.beta^2)*sol.2D.minus30$var.beta +
                 3*sol.2D.minus30$var.beta^2 - sol.2D.minus30$mean.beta2^2

plot(1:n.yrs,sol.2D.minus20$var.beta2,col="blue",type="l",ylim=c(0,2000),lwd=4,ylab="Variance of beta^2")
points(1:n.yrs,minus20.var,col="red",type="l",lty=2)  # Gaussian approx 

points(1:n.yrs,sol.2D.minus30$var.beta2,col="turquoise",type="l",lwd=4)
points(1:n.yrs,minus30.var,col="red",type="l",lty=2) # Gaussian approx 
legend("right", c("Gaussian Approx","delta method"), lty=2,col=c("red","purple"),bty="n"); 
add_panel_label("d")

minus20.dvar <- 4*(sol.2D.minus20$mean.beta^2)*sol.2D.minus20$var.beta; 
minus30.dvar <- 4*(sol.2D.minus30$mean.beta^2)*sol.2D.minus30$var.beta; 
matpoints(1:n.yrs,cbind(minus20.dvar,minus30.dvar),type="l",lty=3,col="purple") 

#the Gaussian assumption looks pretty good 

#####################################################################################
# FOURTH PLOT: NEW VARIANCE APPROX
#####################################################################################
nBigMatrix <- nBigMatrix.z; 

### params,init.mean.flow.int,n.iter,delta,var.beta.t,prop.Recr
beta.minus20 <- Approx_dynamics_Var(m.par.true,-20,n.yrs,0.00001,sol.2D.minus20$var.beta,sol.2D.minus20$prop.Recr)
beta.minus30 <- Approx_dynamics_Var(m.par.true,-30,n.yrs,0.00001,sol.2D.minus30$var.beta,sol.2D.minus30$prop.Recr)

#quick check it all adds up
set_graph_pars("panel2");

sum.t <- beta.minus20$change.P.w+beta.minus20$change.F.w +beta.minus20$P.change.w+beta.minus20$F.change.w
plot(beta.minus20$dlam,sum.t)

cor.test(beta.minus20$dlam,sum.t)
max(abs(beta.minus20$dlam-sum.t),na.rm=TRUE)

sum.t <- beta.minus30$change.P.w+beta.minus30$change.F.w +beta.minus30$P.change.w+beta.minus30$F.change.w
plot(beta.minus30$dlam,sum.t)

cor.test(beta.minus30$dlam,sum.t)
max(abs(beta.minus30$dlam-sum.t),na.rm=TRUE)

#looks OK

dev.new(); ### start plot 
set_graph_pars("panel2"); 
matplot(mean.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",
	ylab="Mean beta",xlim=c(1,n.yrs)) 
points(1:n.yrs,sol.2D.minus20$mean.beta,col="blue",type="l",lty=1,lwd=4)
points(1:n.yrs,sol.2D.minus30$mean.beta,col="turquoise",type="l",lty=1,lwd=4)
points(1:n.yrs,beta.minus20$beta.t.mean,col="red",type="l",lty=2,lwd=2)
points(1:n.yrs,beta.minus30$beta.t.mean,col="red",type="l",lty=2,lwd=2)
legend("topright",legend=c("IPM -20","IPM -30","Approx"),col=c("blue","turquoise","red"),
    lty=1,bty="n",lwd=c(3,3,2));  

add_panel_label("a")


matplot(var.flow.ints,type="l",col=rep(c("black","grey"),rep(5,2)),xlab="Time",
	ylab="Variance beta",xlim=c(1,n.yrs),ylim=c(0,1))
points(1:n.yrs,sol.2D.minus20$var.beta,col="blue",type="l",lwd=3)
points(1:n.yrs,sol.2D.minus30$var.beta,col="turquoise",type="l",lwd=3)

points(1:n.yrs,beta.minus20$beta.t.var,col="purple",type="l",lwd=2,lty=2)
points(1:n.yrs,beta.minus30$beta.t.var,col="red",type="l",lwd=2,lty=2)
add_panel_label("b")

set_graph_pars("panel4"); 

plot(1:n.yrs,beta.minus20$dlam,type="l",ylim=c(-0.02,0.05),lwd=2,xlab="Time",
	ylab="Selection",main="-20 start")
points(1:n.yrs,beta.minus20$dsurv,col="red",type="l")
points(1:n.yrs,beta.minus20$dfec,col="turquoise",type="l")
abline(h=0)

legend("topright",legend=c("Total","Surv","Fec"),col=c("black","red","turquoise"),
    lty=1,bty="n",lwd=c(2,1,1))

plot(1:n.yrs,beta.minus30$dlam,type="l",ylim=c(-0.02,0.05),lwd=2,xlab="Time",
	ylab="Selection",main="-30 start")
points(1:n.yrs,beta.minus30$dsurv,col="red",type="l")
points(1:n.yrs,beta.minus30$dfec,col="turquoise",type="l")
abline(h=0)

legend("topright",legend=c("Total","Surv","Fec"),col=c("black","red","turquoise"),
    lty=1,bty="n",lwd=c(2,1,1))

plot(1:n.yrs,beta.minus20$dlam,type="l",ylim=c(-0.15,0.15),lwd=2,xlab="Time",
	ylab="Selection",main="-20 start")
points(1:n.yrs,beta.minus20$change.P.w+beta.minus20$change.F.w,col="red",type="l")
points(1:n.yrs,beta.minus20$P.change.w+beta.minus20$F.change.w,col="turquoise",type="l")
abline(h=0)

legend("topright",legend=c("Total","Kernel","w"),col=c("black","red","turquoise"),
    lty=1,bty="n",lwd=c(2,1,1))
    
 plot(1:n.yrs,beta.minus30$dlam,type="l",ylim=c(-0.15,0.15),lwd=2,xlab="Time",
	ylab="Selection",main="-30 start")
points(1:n.yrs,beta.minus30$change.P.w+beta.minus30$change.F.w,col="red",type="l")
points(1:n.yrs,beta.minus30$P.change.w+beta.minus30$F.change.w,col="turquoise",type="l")
abline(h=0)

legend("topright",legend=c("Total","Kernel","w"),col=c("black","red","turquoise"),
    lty=1,bty="n",lwd=c(2,1,1))

set_graph_pars("panel2"); 

plot(1:n.yrs,beta.minus20$dlam,type="l",ylim=c(-0.15,0.15),lwd=2,xlab="Time",
	ylab="Selection",main="-20 start")
points(1:n.yrs,beta.minus20$change.P.w,col="red",type="l",lwd=2)
points(1:n.yrs,beta.minus20$change.F.w,col="turquoise",type="l",lwd=2)
points(1:n.yrs,beta.minus20$P.change.w,col="red",type="l")
points(1:n.yrs,beta.minus20$F.change.w,col="turquoise",type="l")
abline(h=0)

legend("topright",legend=c("Total","partial-P","partial-F","P-partial.w","F.partial.w"),col=c("black","red","turquoise","red","turquoise"),
    lty=1,bty="n",lwd=c(2,2,2,1,1))

    
 plot(1:n.yrs,beta.minus30$dlam,type="l",ylim=c(-0.15,0.15),lwd=2,xlab="Time",
	ylab="Selection",main="-30 start")
points(1:n.yrs,beta.minus30$change.P.w,col="red",type="l",lwd=2)
points(1:n.yrs,beta.minus30$change.F.w,col="turquoise",type="l",lwd=2)
points(1:n.yrs,beta.minus30$P.change.w,col="red",type="l")
points(1:n.yrs,beta.minus30$F.change.w,col="turquoise",type="l")
abline(h=0)

legend("topright",legend=c("Total","partial-P","partial-F","P-partial.w","F.partial.w"),col=c("black","red","turquoise","red","turquoise"),
    lty=1,bty="n",lwd=c(2,2,2,1,1))





# plot(1:n.yrs,sol.2D.minus20$mean.beta2,col="blue",type="l",ylim=c(-1000,1000),
	# ylab="Mean of beta^2",lwd=4,lty=1)
# points(1:n.yrs,beta.minus20$S.t.mean,col="red",type="l",lwd=2,lty=2)

# points(1:n.yrs,sol.2D.minus30$mean.beta2,col="turquoise",type="l",lwd=4,lty=1)
# points(1:n.yrs,beta.minus30$S.t.mean,col="red",type="l",lwd=2,lty=2)
# add_panel_label("c")

#Again starting at -30 and it's fine at -20 it's not!

# #calculate variance of S assuming Gaussian distribution
                 
# plot(1:n.yrs,sol.2D.minus20$var.beta2,col="blue",type="l",ylim=c(0,1000),ylab="Variance of beta^2",
	# lty=1,lwd=4)
# points(1:n.yrs,beta.minus20$beta.t.var,col="red",type="l",lty=2)
# points(1:n.yrs,sol.2D.minus30$var.beta2,col="turquoise",type="l",lty=1,lwd=4)
# points(1:n.yrs,beta.minus30$beta.t.var,col="red",type="l",lty=2)
# add_panel_label("d")

#dev.copy2eps(file="~/Repos/ipm_book/c9/figures/OenotheraVarDyn.eps")

#




