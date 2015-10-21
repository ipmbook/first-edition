##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IPM to illustrate and test "Events in the Life cycle" analyses
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls(all=TRUE))

library(doBy)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
setwd("~/Repos/ipm_book/Rcode/c2")


source("Monocarp Demog Funs.R");
source("Monocarp Simulate IBM.R") 
cat(pop.size.t,"\n"); 

## trim an initial transient off the simulation
sim.data <- sim.data[sim.data$yr > 10,];
sim.data$yr <- sim.data$yr-10;

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Construct IPM kernels; assumes the IBM has been run and results stored in sim.data 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nBigMatrix <- 250; 
m.par.true <- m.par; 

min.size <- with(sim.data, min(z))
max.size <- with(sim.data, max(z))
IPM.true <- mk_K(nBigMatrix, m.par.true)

meshpts <- IPM.true$meshpts
h <- diff(meshpts)[1]

lambda0 <- Re(eigen(IPM.true$K)$values[1])
fit.pop.growth <- lm(log(pop.size.t[-(1:10)])~c(11:yr))
lambda.hat <- exp(coef(fit.pop.growth)[2])
cat("Growth rate: IPM ",lambda0,"  sim ",lambda.hat, "\n")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## IPM kernels for the calculations 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# to keep close to the formulae in the text next we define the F and P iteration matricies
P <- IPM.true$P;  F <- IPM.true$F;

# Fundamental operator 
N <- solve(diag(nBigMatrix)-P); 

# Compute R0 as dominant eigenvalue of FN
R0 <- abs(eigen(F %*% N)$values[1])
cat("R0=",R0,"\n");

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Some useful vectors 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define a probability distribution for offspring sizes
n_0 <- h*c_0z1(meshpts,m.par); 
sum(n_0); # should be 1, and it is. 

#Define the e vector which we will use for summing down columns
e=matrix(1,nrow=1,ncol=nBigMatrix)

# With mixing at birth, we can calculate R0 directly as the 
# mean per-capita lifetime fecundity across a birth cohort
R0 <- sum((e %*% F %*% N) * n_0)
cat("R0=",R0,"\n");


#Lots of the calculations involve multiplying e by P let's see what this does
cat(round((e %*% P),4)[1:10],"\n");
# So we're summing down the columns of P, which is equivalent of summing 
# over all the states an individual can move to. 
# As everyone that survives goes somewhere 
# this must equal the probability of survival (remembering that 
# flowering results in death), so
cat(round(s_z(meshpts,m.par.true)*(1-p_bz(meshpts,m.par.true)),4)[1:10],"\n");
#does the same.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Age-specific vital rates: 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#alive means you didn't flower or die (Repr==0 & Surv==1), so this gives the probability of
#survival at each age

# summaryBy(alive~age,data=sim.data)
# 
# #How many individuals in each age class?
# summaryBy(alive~age,FUN=length,data=sim.data)
# #So we have plenty up to about age 5

# vector to hold survivorship curve 
# la is the probability of surviving to age a
la <- rep(NA,20) 

# We can calculate la[1], survival to age 1, as
la[1] <- sum((e %*% P)*n_0);

#Later survivorships require P^a so let's do the calculation recursively
Pa <- P
for(a in 2:20){
	Pa=Pa %*% P
	la[a]= sum((e %*% Pa)*n_0)
}
la

# to calculate the probability of survival from age a to age a+1 
# we just calculate la[a+1]/la[a]
pa <- la; pa[2:20] <- la[2:20]/la[1:19]
pa
# which compares well with the simulation data (plotted below) 
sim.pa <- summaryBy(alive~age,data=sim.data,keep.names=TRUE)

# Let's calculate the age-specific expected fecundity, fa
Pa <- P
fa=rep(NA,20)
fa.0=sum((e %*% F)*n_0)
fa[1]=sum((e %*% F %*% P)*n_0)
for(a in 2:20){
	Pa=Pa %*% P
	fa[a]= sum((e %*% F %*% Pa)*n_0)
}
fa <- fa/la
round(fa,4)

# Compute age-specific fecundity in simulation 
sim.data$Recruits = m.par.true["p.r"]*sim.data$Seeds
sim.data$Recruits[is.na(sim.data$Seeds)] <- 0
sim.fa <- summaryBy(Recruits~age,data=sim.data,keep.names=TRUE)
sim.fa
summaryBy(Repr~age,FUN=sum,data=sim.data)

### Compare theoretical pa and fa with simulation results 
dev.new(width=6,height=4)
par(mfrow=c(1,2), bty="l", pty="s", pch=19,mgp=c(2.5,1,0),mar=c(4,4,1,1))
plot(alive ~ age, data = sim.pa[1:6,],ylim=c(0,1),
	xlab=expression("Age"),
	ylab=expression("p"[a]))
lines(0:5,pa[1:6])
mtext(side = 3, line=0.5, adj = 0, text = "A)")

plot(Recruits ~ age, data = sim.fa[1:9,],
	xlab=expression("Age"),
	ylab=expression("f"[a]))

lines(0:8,c(fa.0,fa[1:8]))
mtext(side = 3, line=0.5, adj = 0, text = "B)")

#dev.copy2eps(file="~/Repos/ipm_book/c3/figures/Oenotherapafa.eps");

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mean and variance in lifespan
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
#First calculate the fundamental matrix N
N <- solve(diag(nBigMatrix)-P)

# Then for the offspring distribution calculate the expected lifespan
# To get the mean age at death we subtract 1 from mean lifespan
# because by our convention lifespan is the number of censuses at which an
# individual is alive, so if you die with lifespan 1 your age at death is 0.    
mean.age.death <- round(sum((e %*% N)*n_0) -1,4)

#check with the simulation, by calculating the mean age at death
mean.age.death.sim <- round(with(sim.data,mean(age[alive==0])),4)

# Now the variance
Var.nu <- round(sum((e %*% (2 * N %*% N - N))* n_0)  - sum((e %*% N * n_0 ))^2,4)


#check with the simulation - adding a constant doesn't change the variance so we don't need to
#subtract 1
Var.nu.sim <- round(with(sim.data,var(age[alive==0])),4)

cat("Age at death","\n")
cat("Theory: mean=",mean.age.death, " variance=",Var.nu,"\n"); 
cat("Simulation: mean=",mean.age.death.sim, " variance=",Var.nu.sim,"\n");

# ymax <- max(sim.data$yr); 
# Vnu.sim <- numeric(ymax)
# for(j in 1:ymax) {
#   Vnu.sim[j] <- with(sim.data,var(age[(alive == 0)& (yr==j)]))
# }  
# 
# par(mfrow=c(1,1),mar=c(5,4,3,1),mgp=c(2.5,1,0)); 
# plot(2:ymax, Vnu.sim[-1],xlab="Length of simulation",ylab="Var(age at death)",
#      ylim=c(0,1.1*Var.nu),bty="l",pch=1);
# abline(h=Var.nu,lty=2); 
# title(main="Monocarp Var(lifespan): sim vs. theory"); 
# 
# py <- Vnu.sim[10:ymax]; 
# px <- 1/(10:ymax); 
# fitVnu <- lm(py~px); 
# points(10:ymax,fitVnu$fitted,type="l",col="red"); 
# cat(Var.nu,fitVnu$coef[1],"\n"); 
# 
# legend("bottomright",legend=c("Simulation","Theory"),
#        lty=c(1,2),pch=c(1,NA),bty="n"); 
# 


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
#  On to the size at death calculations - remember there are 2 ways of dying! 
#  3hrs later I (MR) remember in monocarps flowering is fatal!
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
omega <-  (meshpts*(p_bz(meshpts,m.par.true) 
          + (1-p_bz(meshpts,m.par.true))*(1-s_z(meshpts,m.par.true)))) %*% N

#  how can you die? You flower with probability p_bz, and if you don't 
#  flower (1-p_bz) you die with probability (1-s_z)
mean.size.death <- sum(omega * n_0)
round(mean.size.death,4)
with(sim.data,mean(z[alive==0])) 

# Let's calculate the size at death kernel
Omega <- (p_bz(meshpts,m.par.true) + (1-p_bz(meshpts,m.par.true))*(1-s_z(meshpts,m.par.true))) * N

# each of the columns defines a probability distribution and so should sum to 1, let's check
round(e %*% Omega,4)
#all 1's as it should be.

# then the distribution of sizes at death for the offspring distribution is
dist.size.death <-(Omega %*% n_0)

# again this should be a probability distribution and so sum to 1, let's check
round(e %*% dist.size.death,4)
# so it's all good, let's calculate some moments

# 1st the mean is
mean.size.death     <-round(sum(dist.size.death*meshpts),4)
mean.size.death.sim <- round(with(sim.data,mean(z[alive==0])),4)    

# 2nd the variance is
Var.size.death <- round(sum(dist.size.death*meshpts*meshpts) - 
          sum(dist.size.death*meshpts)*sum(dist.size.death*meshpts),4)
          
Var.size.death.sim <- round(with(sim.data,var(z[alive==0])),4)

cat("Size at death","\n")
cat("Theory: mean=",mean.size.death, " variance=",Var.size.death,"\n"); 
cat("Simulation: mean=",mean.size.death.sim, " variance=",Var.size.death.sim,"\n");


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    Reproduction: who, when and how much? 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#As reproduction is fatal the P kernel is the required P0 kernel
P0 <- P
N0 <- solve(diag(nBigMatrix)-P0)

#So B is given by
B <- p_bz(meshpts,m.par.true) %*% N0

# plot(meshpts,B,type="l")
# points(meshpts,p_bz(meshpts,m.par.true),type="l",col="red")


# Breeding-weighted offspring size distribution
breed.n_0 = B*n_0/sum(B*n_0)

B.m.c <- matrix(B,nrow=nBigMatrix,ncol=nBigMatrix)
B.m.r <- matrix(B,nrow=nBigMatrix,ncol=nBigMatrix,byrow=TRUE)
P.b <- (P0 * B.m.c ) / B.m.r
N.b <- solve(diag(nBigMatrix)-P.b)
mean.Repr.age <- sum((e %*% N.b) * breed.n_0)-1

# Compare with simulations
cat(mean.Repr.age,"\n"); 
cat(with(sim.data,mean(age[Repr == 1])),"\n"); 

ymax <- max(sim.data$yr); 
abarR <- numeric(ymax)
for(j in 1:ymax) {
  abarR[j] <- with(sim.data,mean(age[(Repr == 1) & (yr<=j)]))
}  

plot(2:ymax,abarR[-1],xlab="Length of simulation",ylab="Mean age at reproduction",
     ylim=c(0,1.1*mean.Repr.age),bty="l");
abline(h=mean.Repr.age,lty=2); 
title(main="Monocarp Eg1: sim vs. theory"); 

py <- abarR[6:ymax]; 
px <- 1/(6:ymax); 
fitAbar <- lm(py~px); 
points(6:ymax,fitAbar$fitted,type="l"); 
cat(mean.Repr.age,fitAbar$coef[1],"\n"); 

legend("bottomright",legend=c("Simulation","Theory"),
       lty=c(1,2),pch=c(1,NA),bty="n"); 


#Distribution of sizes at reproduction    
Omega.b <- as.vector(1-(e %*% P.b)) * N.b
dist.size.repr <- (Omega.b %*% n_0)

#mean size at reproduction
mean.size.flowering <- sum(h*dist.size.repr*meshpts)
mean.size.flowering.sim <- with(sim.data,mean(z[Repr == 1]))

#variance in size at reproduction
var.size.flowering <- sum(h*dist.size.repr*meshpts*meshpts)-sum(h*dist.size.repr*meshpts)^2
var.size.flowering.sim <- with(sim.data,var(z[Repr == 1]))

cat("Size at reproduction","\n")
cat("Theory: mean=",mean.size.flowering, " variance=",var.size.flowering,"\n"); 
cat("Simulation: mean=",mean.size.flowering.sim, " variance=",var.size.flowering.sim,"\n");

#How often do they reproduce?
breedingFreq <- p_bz(meshpts,m.par.true) %*% N0/B
cat("Mean breeding frequency (theory, should =1):", range(breedingFreq),"\n"); 
# well they're monocarpic and reproduction is fatal...


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## Plots 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
dev.new(6,6)
par(mfrow=c(2,2), bty="l", pty="s", pch=19)

## 1 - plot population density versus time...
plot(1:yr, mean.z.death.t [1:yr], type="l",xlab="Time",ylab="Mean size at death")
abline(h=mean.size.death, col="red")
mtext(side = 3, line=0.5, adj = 0, text = "a)")
## ...roughly linear for log(Nt) vs time so exponential growth


## 2 - plot mean size versus time...
plot(1:yr, mean.age.death.t[1:yr], type="l",xlab="Time",ylab="Mean age at death")
abline(h=mean.age.death, col="red")
mtext(side = 3, line=0.5, adj = 0, text = "b)")
## ...mean size seems to settle down after a while

## 3 - plot mean flowering size versus time...
plot(1:yr, mean.fl.z.t[1:yr], type="l",xlab="Time",ylab="Mean size at flowering",ylim=c(0,4))
abline(h=mean.size.flowering, col="red")
mtext(side = 3, line=0.5, adj = 0, text = "c)")
## ...mean flowering size seems to settle down after a while

## 4 - plot of density estimates at time 50 and the end 
plot(1:yr, mean.fl.age.t[1:yr], type="l",xlab="Time",ylab="Mean age at flowering",ylim=c(0,4))
abline(h=mean.Repr.age, col="red")
mtext(side = 3, line=0.5, adj = 0, text = "d)")


} #END IF(FALSE) 












