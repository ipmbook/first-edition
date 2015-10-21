##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IPM to illustrate and test "Events in the Life cycle" analyses
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls(all=TRUE))

library(doBy)
library(parallel)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
mainDir<-paste0(root,"/ipm_book/Rcode/c3")
setwd(mainDir); 

source("../c2/Monocarp Demog Funs.R");
source("../utilities/Standard Graphical Pars.R");


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## Generate artificial data by simulating IBM from Chapter 2 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m.par.true["p.r"] <- 0.005325442  #changed p.r so R0=1
init.pop.size <- 100000
n.yrs <-50
source("../c2/Monocarp Simulate IBM.R") 
cat(pop.size.t,"\n"); 

## trim an initial transient off the simulation
sim.data <- sim.data[sim.data$yr > 10,];
sim.data$yr <- sim.data$yr-10;

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Construct IPM kernels; assumes the IBM has been run and results stored in sim.data 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nBigMatrix <- 250; 

min.size <- with(sim.data, min(z))
max.size <- with(sim.data, max(z))
IPM.true <- mk_K(nBigMatrix, m.par.true, -3.5, 5.5)

meshpts <- IPM.true$meshpts
h <- diff(meshpts)[1]

lambda0 <- Re(eigen(IPM.true$K,only.values = TRUE)$values[1])
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
R <- F %*% N
R0 <- abs(eigen(R)$values[1])
cat("R0=",R0,"\n");

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Some useful vectors 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define a probability distribution for offspring sizes
offspring.prob <- h*c_0z1(meshpts,m.par.true); 
sum(offspring.prob); # should be 1, and it is. 

#Define the e vector which we will use for summing down columns
e=matrix(1,nrow=1,ncol=dim(P)[1])

# With mixing at birth, we can calculate R0 directly as the 
# mean per-capita lifetime fecundity across a birth cohort
R0 <- sum((e %*% R) * offspring.prob)
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

 summaryBy(alive~age,data=sim.data)
# 
# #How many individuals in each age class?
 summaryBy(alive~age,FUN=length,data=sim.data)
# #So we have plenty up to about age 5

# vector to hold survivorship curve 
# la is the probability of surviving to age a
la <- rep(NA,20)

# We can calculate la[1], survival to age 1, as
la[1] <- sum((e %*% P)*offspring.prob);

#Later survivorships require P^a so let's do the calculation recursively
Pa <- P
for(a in 2:20){
	Pa=Pa %*% P
	la[a]= sum((e %*% Pa)*offspring.prob)
}

la <- c(1,la)

# to calculate the probability of survival from age a to age a+1 
# we calculate la[a+1]/la[a]
pa<- la[2:21]/la[1:20]

# which compares well with the simulation data (plotted below) 
sim.pa <- summaryBy(alive~age,data=sim.data,keep.names=TRUE)

# Let's calculate the age-specific expected fecundity, fa
Pa <- P
fa=rep(NA,20)
fa.0=sum((e %*% F)*offspring.prob)
fa[1]=sum((e %*% F %*% P)*offspring.prob)
for(a in 2:20){
	Pa=Pa %*% P
	fa[a]= sum((e %*% F %*% Pa)*offspring.prob)
}

fa <- c(fa.0,fa)

fa <- fa/la
round(fa,4)

# Compute age-specific fecundity in simulation 
sim.data$Recruits = m.par.true["p.r"]*sim.data$Seeds
sim.data$Recruits[is.na(sim.data$Seeds)] <- 0
sim.fa <- summaryBy(Recruits~age,data=sim.data,keep.names=TRUE)
sim.fa
#summaryBy(Repr~age,FUN=sum,data=sim.data)

### Compare theoretical pa and fa with simulation results 
dev.new(height=4,width=8); set_graph_pars("panel2"); 
plot(alive ~ age, data = sim.pa[1:10,],ylim=c(0,1),pch=19,
	xlab=expression("Age"),
	ylab=expression("p"[a]))
lines(0:9,pa[1:10])
add_panel_label("a")

plot(Recruits ~ age, data = sim.fa[1:10,],pch=19,
	xlab=expression("Age"), 
	ylab=expression("f"[a]))

lines(0:9,fa[1:10])
add_panel_label("b")

#dev.copy2eps(file="../../c3/figures/Oenotherapafa.eps",colormodel="cmyk");

#Generation time calculations

T1 <- log(R0)/log(lambda0)
T2 <- sum((1:21)*fa*la)/sum(fa*la)
cat("Generation time log(R0)/log(lambda) = ",T1,"\n")
cat("Generation time from age-specific rates = ",T2,"\n")

w <- Re(abs(eigen(IPM.true$K)$vectors[,1]))

off.dist <- F %*% w / (sum(F %*% w))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mean and variance in lifespan
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# For the offspring distribution calculate the expected lifespan
# To get the mean age at death we subtract 1 from mean lifespan
# because by our convention lifespan is the number of censuses at which an
# individual is alive, so if you die with lifespan 1 your age at death is 0.    
mean.age.death <- round(sum((e %*% N)*offspring.prob) -1,4)

#check with the simulation, by calculating the mean age at death
mean.age.death.sim <- round(with(sim.data,mean(age[alive==0])),4)

# Now the variance
Var.nu <- round(sum((e %*% (2 * N %*% N - N))* offspring.prob)  - sum((e %*% N * offspring.prob ))^2,4)

# check with the simulation - adding a constant doesn't change the variance so we don't need to
# subtract 1
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
#  3hrs later I (MR) remembers in monocarps flowering is fatal!
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
omega <-  (meshpts*(p_bz(meshpts,m.par.true) 
          + (1-p_bz(meshpts,m.par.true))*(1-s_z(meshpts,m.par.true)))) %*% N

#  how can you die? You flower with probability p_bz, and if you don't 
#  flower (1-p_bz) you die with probability (1-s_z)
mean.size.death <- sum(omega * offspring.prob)
round(mean.size.death,4)
with(sim.data,mean(z[alive==0])) 

# Let's calculate the size at death kernel
Omega <- (p_bz(meshpts,m.par.true) + (1-p_bz(meshpts,m.par.true))*(1-s_z(meshpts,m.par.true))) * N

# each of the columns defines a probability distribution and so should sum to 1, let's check
round(e %*% Omega,4)
#all 1's as it should be.

# then the distribution of sizes at death for the offspring distribution is
dist.size.death <- (Omega %*% offspring.prob)

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


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    Reproduction: who, when and how much? 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#As reproduction is fatal the P kernel is the required P0 kernel
P0 <- P
N0 <- solve(diag(nBigMatrix)-P0)

#So B is given by
B <- p_bz(meshpts,m.par.true) %*% N0

# plot(meshpts,B,type="l")
# points(meshpts,p_bz(meshpts,m.par.true),type="l",col="red")


# Breeding-weighted offspring size distribution
breed.offspring.prob <- B*offspring.prob/sum(B*offspring.prob)

B.m.c <- matrix(B,nrow=nBigMatrix,ncol=nBigMatrix)
B.m.r <- matrix(B,nrow=nBigMatrix,ncol=nBigMatrix,byrow=TRUE)
P.b <- (P0 * B.m.c ) / B.m.r
N.b <- solve(diag(nBigMatrix)-P.b)
mean.Repr.age <- sum((e %*% N.b) * breed.offspring.prob)-1

# Compare with simulations
cat("Age at reproduction","\n")
cat("Theory: mean=",mean.Repr.age,"\n"); 
cat("Simulation: mean=",with(sim.data,mean(age[Repr == 1])),"\n");

# 
# ymax <- max(sim.data$yr); 
# abarR <- numeric(ymax)
# for(j in 1:ymax) {
#   abarR[j] <- with(sim.data,mean(age[(Repr == 1) & (yr<=j)]))
# }  
# 
# plot(2:ymax,abarR[-1],xlab="Length of simulation",ylab="Mean age at reproduction",
#      ylim=c(0,1.1*mean.Repr.age),bty="l");
# abline(h=mean.Repr.age,lty=2); 
# title(main="Monocarp Eg1: sim vs. theory"); 
# 
# py <- abarR[6:ymax]; 
# px <- 1/(6:ymax); 
# fitAbar <- lm(py~px); 
# points(6:ymax,fitAbar$fitted,type="l"); 
# cat(mean.Repr.age,fitAbar$coef[1],"\n"); 
# 
# legend("bottomright",legend=c("Simulation","Theory"),
#        lty=c(1,2),pch=c(1,NA),bty="n"); 


# Distribution of sizes at reproduction    
Omega.b <- as.vector(1-(e %*% P.b)) * N.b
dist.size.repr <- (Omega.b %*% t(breed.offspring.prob))

#mean size at reproduction
mean.size.flowering <- sum(dist.size.repr*meshpts)
mean.size.flowering.sim <- with(sim.data,mean(z[Repr == 1]))

#variance in size at reproduction
var.size.flowering <- sum(dist.size.repr*meshpts*meshpts)-sum(dist.size.repr*meshpts)^2
var.size.flowering.sim <- with(sim.data,var(z[Repr == 1]))

cat("Size at reproduction","\n")
cat("Theory: mean=",mean.size.flowering, " variance=",var.size.flowering,"\n"); 
cat("Simulation: mean=",mean.size.flowering.sim, " variance=",var.size.flowering.sim,"\n");

#How often do they reproduce?
breedingFreq <- p_bz(meshpts,m.par.true) %*% N0/B
cat("Mean breeding frequency (theory, should =1):", range(breedingFreq),"\n"); 
# well they're monocarpic and reproduction is fatal...

#Variance in life-time reproductive output
#For a Poisson distribution the mean equals the variance, so for the number of recruits (not seeds)

mean.recs <-  p_bz(meshpts,m.par.true) * m.par.true["p.r"] *  b_z(meshpts,m.par.true)
#sigma.b.2 <-  m.par.true["p.r"] * p_bz(meshpts,m.par.true) * b_z(meshpts,m.par.true)
sigma.b.2 <-  mean.recs + (1 - p_bz(meshpts,m.par.true)) * mean.recs * mean.recs / p_bz(meshpts,m.par.true)

rbar <- e %*% R
r2 <- (sigma.b.2 + (m.par.true["p.r"] * p_bz(meshpts,m.par.true) * b_z(meshpts,m.par.true))^2) %*% N
var.b <- r2 - rbar*rbar

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# now calculate mean and variance in reproductive output using an IBM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
#  Version 1: use a for loop 
 mean.var.Repr <- matrix(NA,nrow=nBigMatrix,ncol=4)
 for(i in 1:nBigMatrix){
 	mean.var.Repr[i,] <- get_mean_var_Repr(meshpts[i],100000)
    cat(i,"\n"); 
}

# Version 2: use sapply
# mean.var.Repr <- t(sapply(1:nBigMatrix, function(i) get.mean.var.Repr(meshpts[i],100000)))

# And to avoid all the waiting - multicore! 
# OS X, mclapply makes it easy  
 mean.var.Repr <- t(simplify2array(mclapply(seq.int(1,nBigMatrix), function(i) get_mean_var_Repr(meshpts[i],100000), 
         mc.cores = detectCores() ) ) )

##### Any OS, foreach works but needs some setup 
# load parallel, foreach, and iterators 
#require(doParallel); require(parallel); 
#
# setup the cluster of cores 
#c1<- makeCluster(detectCores()-2); registerDoParallel(c1);  
#
## export the data to the clusters; sending it all is quite lazy but works 
#clusterExport(cl=c1, varlist=ls(all=TRUE), envir = .GlobalEnv)
#
## do the work, finally. Each i computes a row of the matrix.  
#mean.var.Repr <- foreach(i=1:nBigMatrix,.combine=rbind) %dopar% {
#       out <- get_mean_var_Repr(meshpts[i],100000)
#}
#stopCluster(c1); 
###### close down the cluster, done with foreach  
#~~~~~~~ Done with IBM calculation of mean and variance ~~~~~~~~~~~~~~~~

dev.new(height=8,width=8) #; par(mar=c(4,4,1,1),mgp=c(2.5,1,0)); 
set_graph_pars("panel4"); 

plot(meshpts,rbar,type="l",xlab=expression("Size t, "*italic(z)), #log="y",
ylab="Mean lifetime reproduction")
points(meshpts,mean.var.Repr[,1],pch=16,cex=0.45)

add_panel_label("a")

plot(meshpts,var.b,type="l",xlab=expression("Size t, "*italic(z)), #log="y",
ylab="Variance lifetime reproduction")
points(meshpts,sigma.b.2,type="l",col="red")
points(meshpts,mean.var.Repr[,2],pch=16,cex=0.45)

add_panel_label("b")

# Didn't expect that, but I think it makes sense, when very large p_b = 1 so you get the Poisson fecundity variance. 
# When very small everyone dies so small variance, in the middle variability builds up as some individuals make it to reproduce
# then drop as everyone reproduces immediately before building up again. 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Next, look at mean and variance in number of reproductive events IBM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


############### IPM calculations 
# Redefine b_z so it equals 1 - remember that in F b_z is multiplied by m.par.true["p.r"]
# This means that when breeding occurs, the number of offspring is 1
# so we can compute mean, var in #offspring for this modified model 
# to get the mean and variance in the number of reproductive events 
b_z <- function(z, m.par){
    N <- 1/m.par.true["p.r"]     # seed production of a size z plant
    return(N)
}

sigma.b.2 <-   p_bz(meshpts,m.par.true) * (1-p_bz(meshpts,m.par.true))
IPM.true <- mk_K(nBigMatrix, m.par.true, -3.5, 5.5)
P <- IPM.true$P;  F <- IPM.true$F;

# Fundamental operator 
N <- solve(diag(nBigMatrix)-P); 
rbar <- e %*% (F %*% N)
r2 <- (sigma.b.2 + (m.par.true["p.r"] * p_bz(meshpts,m.par.true) * b_z(meshpts,m.par.true))^2) %*% N
var.b <- r2 - rbar*rbar

############### IBM calculations  

# Another cheat: modify b_z so that breeding results in at least one seed
# with probability nearly 1. We can then use our same get_mean_var_Repr
# function, with (Seeds>0) now being equivalent to breeding. 
b_z <- function(z,m.par) 5000 

mean.var.Repr2 <- matrix(NA,nrow=nBigMatrix,ncol=4)
 for(i in 1:nBigMatrix){
 	mean.var.Repr2[i,] <- get_mean_var_Repr(meshpts[i],100000)
    cat(i,"\n"); 
}

# as above, this could be done on multiple cores
# mean.var.Repr2 <- t(simplify2array(mclapply(seq.int(1,nBigMatrix), function(i) get_mean_var_Repr(meshpts[i],100000), mc.cores =detectCores())))

plot(meshpts,rbar,type="l",xlab=expression("Size t, "*italic(z)), #log="y",
ylab="Mean reproductive events")
points(meshpts,mean.var.Repr2[,3],pch=16,cex=0.45)

add_panel_label("c")

plot(meshpts,var.b,type="l",xlab=expression("Size t, "*italic(z)), #log="y",
ylab="Variance reproductive events")
points(meshpts,mean.var.Repr2[,4],pch=16,cex=0.45)

add_panel_label("d")

#dev.copy2eps(file="../../c3/figures/OenotheraMeanVarRepr.eps",colormodel="cmyk");

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Next, look at mean and variance of age at last reproduction. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Right, set the b_z function equal to it's original form - Poisson seed production
b_z <- function(z, m.par) {
    N <- exp(m.par["seed.int"] + m.par["seed.z"] * z)    # seed production of a size z plant
    return(N)
}

IPM.true <- mk_K(nBigMatrix, m.par.true, -3.5, 5.5)
P <- IPM.true$P;  F <- IPM.true$F;

# Fundamental operator 
N <- solve(diag(nBigMatrix)-P); 
# 
# matrix.p.bz <- matrix(1-p_bz(meshpts,m.par.true),nrow=nBigMatrix,ncol=nBigMatrix,byrow=TRUE)
# 
# pie.0 <- P/matrix.p.bz
# 
# d0 <- 1-e %*% pie.0
# 
# A0 <- pie.0 * t(matrix.p.bz)
# 
# f0 <- d0 %*% (solve(diag(nBigMatrix)-A0))
# 
# fb <- rep(1,nBigMatrix)

B <- p_bz(meshpts,m.par.true) %*% N

E.bS <- (p_bz(meshpts,m.par.true) %*% (N %*% N - N))/B

mean.age.last.breed <- sum(breed.offspring.prob*E.bS)

mean.age.last.breed.sim <- with(sim.data,mean(age[Repr == 1]))

E.bS2 <- (p_bz(meshpts,m.par.true) %*% (P + P %*% P) %*% (N %*% N %*% N)) / B

var.age.last.breed <- sum(breed.offspring.prob*E.bS2) - (sum(breed.offspring.prob*E.bS)^2)

var.age.last.breed.sim <- with(sim.data,var(age[Repr == 1]))

cat("Age last reproduction","\n")
cat("Theory: mean=",mean.age.last.breed, " variance=",var.age.last.breed,"\n"); 
cat("Simulation: mean=",mean.age.last.breed.sim, " variance=",var.age.last.breed.sim,"\n");
 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## Plots - just checking 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
dev.new(height=8,width=8) 
set_graph_pars("panel4")

## 1 - plot population density versus time...
plot(1:yr, mean.z.death.t [1:yr], type="l",xlab="Time",ylab="Mean size at death")
abline(h=mean.size.death, col="red")
mtext(side = 3, line=0.5, adj = 0, text = "A)")
## ...roughly linear for log(Nt) vs time so exponential growth

## 2 - plot mean size versus time...
plot(1:yr, mean.age.death.t[1:yr], type="l",xlab="Time",ylab="Mean age at death")
abline(h=mean.age.death, col="red")
mtext(side = 3, line=0.5, adj = 0, text = "B)")
## ...mean size seems to settle down after a while

## 3 - plot mean flowering size versus time...
plot(1:yr, mean.fl.z.t[1:yr], type="l",xlab="Time",ylab="Mean size at flowering",ylim=c(2,3.5))
abline(h=mean.size.flowering, col="red")
mtext(side = 3, line=0.5, adj = 0, text = "C)")
## ...mean flowering size seems to settle down after a while

## 4 - plot of density estimates at time 50 and the end 
plot(1:yr, mean.fl.age.t[1:yr], type="l",xlab="Time",ylab="Mean age at flowering",ylim=c(0,4))
abline(h=mean.Repr.age, col="red")
mtext(side = 3, line=0.5, adj = 0, text = "D)")












