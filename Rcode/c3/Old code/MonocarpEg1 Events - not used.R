## IBM to generate the data for the simple recruitment limited IPM - there's no recruitment limitation 1st
## Size is on a log scale so don't worry about the -ive ones!

## The parameters and life cycle follow Kachi and Hirose J.Ecol. 1985, see Rees and Rose PRSB 2002, but for
## simplicity we have used a logistic function for survival rather than the linear function with an upper bound
## as in the previous papers

## MORE COMMENTS STILL TO BE ADDED

rm(list=ls(all=TRUE))

library(doBy)
set.seed(53241986)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##
## Define the true parameter vector, each parameter is given a name so the formulae, below, are easier to
## read. We'll use 'model.term' scheme to name elements of the vector
##

m.par.true <- c(surv.int  =  -0.65,
                surv.z    =   0.75,
                flow.int  = -18.00,
                flow.z    =   6.9,
                grow.int  =   0.96,
                grow.z    =   0.59,
                grow.sd   =   0.67,
                rcsz.int  =  -0.08,
                rcsz.sd   =   0.76,
                seed.int  =   1.00,
                seed.z    =   2.20,
                p.r       =   0.005325442)  #changed p.r so R0=1

##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets, say the true values and
## estimated ones
##

## Growth function, given you are size z now returns the pdf of size z1 next time

g_z1z <- function(z1, z, m.par)
{
    mean <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year
    sd <- m.par["grow.sd"]                                    # sd about mean
    p.den.grow <- dnorm(z1, mean = mean, sd = sd)             # pdf that you are size z1 given you were size z
    return(p.den.grow)
}

## Survival function, logistic regression

s_z <- function(z, m.par)
{
    linear.p <- m.par["surv.int"] + m.par["surv.z"] * z  # linear predictor
    p <- 1/(1+exp(-linear.p))                            # logistic transformation to probability
    return(p)
}

## Probability of flowering function, logistic regression

p_bz <- function(z, m.par)
{
    linear.p <- m.par["flow.int"] + m.par["flow.z"] * z      # linear predictor
    p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
    return(p)
}

## Seed production function

b_z <- function(z, m.par)
{
    N <- exp(m.par["seed.int"] + m.par["seed.z"] * z)    # seed production of a size z plant
    return(N)
}

## Recruit size pdf

c_z1 <- function(z1, m.par)
{
    mean <- m.par["rcsz.int"]
    sd <- m.par["rcsz.sd"]
    p.deRecr <- dnorm(z1, mean = mean, sd = sd)              # pdf of a size z1 recruit
    return(p.deRecr)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##
## Parameters controlling the IBM simulation
##

m.par <- m.par.true

init.pop.size <- 1000000
n.yrs <-100
z <- rnorm(init.pop.size, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])
age <- rep(0,init.pop.size)

## Calculate initial pop size and mean size

pop.size.t <- init.pop.size
mean.z.t <- mean(z)
mean.age.t <- mean(age)


## iterate the model using the "true" parameters and store data in a data.frame
yr <- 1
while(yr != n.yrs & length(z) < 1500000) {

    ## calculate population size
    pop.size <- length(z)

    ## generate binomial random number for the probability of flowering, where the probability of flowering
    ## depends on your size z, this is a vector of 0's and 1's, you get a 1 if you flower
    Repr <- rbinom(n=pop.size, prob=p_bz(z, m.par), size=1)

    ## number of plants that flowered
    num.Repr <- sum(Repr)

    ## calculate seed production
    Seeds <- rep(NA, pop.size)

    ## we'll assume plant make a Poisson distributed number of seeds with a mean given by
    ## exp(params["seed.int"]+params["seed.size"] * z)
    ## rpois generated Poisson distributed random numbers
    Seeds[Repr==1] <- rpois(num.Repr, b_z(z[Repr==1], m.par))

    ## generate the number of recruits
    Recr <- if (num.Repr==0) 0 else rbinom(1,sum(Seeds, na.rm=TRUE),m.par["p.r"])

    ## generate new recruit sizes
    ## rnorm generated normally distributed random numbers
    Rcsz <- rnorm(Recr, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])

    ## for the non-reproductive plants generate random number for survival
    Surv <- rep(NA, pop.size)
    Surv[Repr==0] <- rbinom(n = pop.size - num.Repr, prob = s_z(z[Repr==0], m.par), size = 1)
    num.die <- sum(Surv==0, na.rm=TRUE)

    ## index for individuals that did not flower and survived
    i.subset <- which(Repr==0 & Surv==1)

    ## let them grow
    z1 <- rep(NA, pop.size)
    E.z1 <- m.par["grow.int"]+m.par["grow.z"]*z[i.subset]
    z1[i.subset] <- rnorm(n = pop.size - num.Repr - num.die, mean = E.z1, sd = m.par["grow.sd"])

    ## store the simulation data, we'll use this later
    sim.data <- data.frame(z=z, Repr, Seeds, Surv, z1=z1, age=age, alive=Repr==0 & Surv==1)

    ## create new population size vector
    z1 <- c(Rcsz, z1[i.subset])
    age1 <- c(rep(0,Recr),age[i.subset]+1)

    pop.size.t          <- c(pop.size.t, length(z1))
    mean.z.death.t      <- if (yr==1) mean(z[!(Repr==0 & Surv==1)]) else 
                              c(mean.z.death.t,mean(z[!(Repr==0 & Surv==1)]))
    mean.age.death.t    <- if (yr==1) mean(age[!(Repr==0 & Surv==1)]) else 
                              c(mean.age.death.t,mean(age[!(Repr==0 & Surv==1)]))
    mean.fl.z.t         <- if (yr==1) mean(z[Repr==1]) else c(mean.fl.z.t,mean(z[Repr==1])) 
	mean.fl.age.t       <- if (yr==1) mean(age[Repr==1]) else c(mean.fl.age.t,mean(age[Repr==1]))
	
	
    z <- z1
    age <- age1
    
    cat(paste(yr, mean.age.death.t [yr], mean.fl.age.t[yr], "\n", sep=" "))

    yr <- yr+1

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 4 - Build an IPM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
p_z1z <- function (z1, z, m.par) {

    return((1 - p_bz(z, m.par)) * s_z(z, m.par) * g_z1z(z1, z, m.par))

}

## Define the fecundity kernel
f_z1z <- function (z1, z, m.par) {

    return( p_bz(z, m.par) * b_z(z, m.par) * m.par["p.r"] * c_z1(z1, m.par))

}

mk_K <- function(m, m.par) {

	# upper and lower integration limits
	L <- min.size - 0.2
	U <- max.size + 0.2

	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	P <- h * (outer(meshpts, meshpts, p_z1z, m.par = m.par))
	F <- h * (outer(meshpts, meshpts, f_z1z, m.par = m.par))
	K <- P + F
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 4 - Projection
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

min.size <- with(sim.data, min(z))
max.size <- with(sim.data, max(z))
IPM.true <- mk_K(100, m.par.true)

Re(eigen(IPM.true$K)$values[1])

fit.pop.growth <- lm(log(pop.size.t)~c(1:yr))

exp(coef(fit.pop.growth)[2])



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 5 - Event calculations
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

R0=abs(Re(eigen(IPM.true$F %*% solve(diag(100)-IPM.true$P))$values[1]))

meshpts <- IPM.true$meshpts
h <- diff(meshpts)[1]

#to keep close to the formulae in the text next we define the F and P iteration matricies

P <- IPM.true$P
F <- IPM.true$F

offspring.dist=c_z1(meshpts,m.par.true)

#Define the e vector which we will use for summing down columns
e=matrix(1,nrow=1,ncol=100)

#As there are 2 ways of dying create a vector describing those that are alive - with 1 for those that didn't die or flower and 0 for those that died

sim.data$alive <- sim.data$Repr==0 & sim.data$Surv==1

#Let's calculate the probability of survival for each age class

summaryBy(alive~age,data=sim.data)

#How many individuals in each age class?

summaryBy(alive~age,FUN=length,data=sim.data)

#So we have plenty up to about age 5

#Let's calculate la for the a cohort starting with the offspring distribution c_z1
#so let's check the offspring distribution sums to 1

h*round((e %*% offspring.dist),4)

#looks good
#Lots of the calculations involve multiplying e by P let's see what this does

round((e %*% P),4)[1:10]

# So we're summing down the columns of P, which is equivalent of summing 
# over all the states an individual can 
# move to, and as everyone that survives goes somewhere 
# this must equal the probability of survival (remembering that 
# flowering results in death), so

round(s_z(meshpts,m.par.true)*(1-p_bz(meshpts,m.par.true)),4)[1:10]

#does the same.


#We can calculate the first la as

sum((e %*% P)*offspring.dist)*h

sum(s_z(IPM.true$meshpts,m.par.true)*(1-p_bz(IPM.true$meshpts,m.par.true))*offspring.dist)*h

#The next terms require P^a so let's do the calculation recursively

Pa <- P

la=rep(NA,12)

la[1]=sum((e %*% P)*offspring.dist)*h

for(a in 2:12){
	Pa=Pa %*% P
	la[a]= sum((e %*% Pa)*offspring.dist)*h
}

la

# la is the probability of surviving to age a, 
# to calculate the probability of survival for a particular age class we
# just calculate la[a+1]/la[a]

qa <- la
qa[2:12] <- la[2:12]/la[1:11]

qa

#which compares well with the simulation data

sim.qa <- summaryBy(alive~age,data=sim.data,keep.names=TRUE)

sim.qa

#Let's calculate the age-specific expected fecundity, fa

Pa <- P

fa=rep(NA,12)

fa.0=sum((e %*% F)*offspring.dist)*h

fa[1]=sum((e %*% F %*% P)*offspring.dist)*h

for(a in 2:12){
	Pa=Pa %*% P
	fa[a]= sum((e %*% F %*% Pa)*offspring.dist)*h
}

fa <- fa/la

round(fa,4)

sim.data$Recruits = m.par.true["p.r"]*sim.data$Seeds

sim.data$Recruits[is.na(sim.data$Seeds)] <-0

sim.fa <- summaryBy(Recruits~age,data=sim.data,keep.names=TRUE)

sim.fa

summaryBy(Repr~age,FUN=sum,data=sim.data)

dev.new(6,6)

par(mfrow=c(1,2), bty="l", pty="s", pch=19)

plot(alive ~ age, data = sim.qa[1:6,],ylim=c(0,1),
	xlab=expression("Age"),
	ylab=expression("q"[a]))

lines(0:5,qa[1:6])

mtext(side = 3, line=0.5, adj = 0, text = "a)")

plot(Recruits ~ age, data = sim.fa[1:9,],
	xlab=expression("Age"),
	ylab=expression("f"[a]))

lines(0:8,c(fa.0,fa[1:8]))

mtext(side = 3, line=0.5, adj = 0, text = "b)")

#Mean and variance in lifespan
#First calculate the fundamental matrix

N <- solve(diag(100)-P)

#then for the offspring distribution calculate the expected lifespan

mean.age.death <- sum((e %*% N)*offspring.dist)*h -1

round(mean.age.death,4)

#check with the simulation, by calculating the mean age at death

round(with(sim.data,mean(age[alive==0])),4)

#we subtract 1 on because by our convention lifespan is the number of censuses at which an
#individual is alive, so if you die with lifespan 1 your age is 0.

#Now the variance

#First the Var(nu) equation from the text
Var.nu <- e %*% (2 * N %*% N - N) - (e %*% N) ^ 2

#then average with respect to the offspring distribution

sum(Var.nu * offspring.dist)*h

#check with the simulation - adding a constant doesn't change the variance so we don't need to
#subtract 1

round(with(sim.data,var(age[alive==0])),4)

#Onto the size at death calculations - remember there are 2 ways of dying! 3hrs later I (MR) remember in monocarps flowering is fatal!

omega <-  (meshpts * (p_bz(meshpts,m.par.true) + (1-p_bz(meshpts,m.par.true))*(1-s_z(meshpts,m.par.true)))) %*% N

#how can you die? You flower with probability p_bz, and if you don't flower (1-p_bz) you die with probability (1-s_z)

mean.size.death <- sum(omega * offspring.dist)*h

round(mean.size.death,4)

with(sim.data,mean(z[alive==0])) 

#Let's calculate the size at death kernel

Omega <- (p_bz(meshpts,m.par.true) + (1-p_bz(meshpts,m.par.true))*(1-s_z(meshpts,m.par.true))) * N

#each of the columns defines a probability distribution and so should sum to 1, let's check

round(e %*% Omega,4)

#all 1's as it should be.

#then the distribution of sizes at death for the offspring distribution is

dist.size.death <-(Omega %*% offspring.dist)

#again this should be a probability distribution and so sum to 1, let's
#check

round(e %*% dist.size.death,4)*h

#so it's all good, let's calculate some moments
#1st the mean is

round(sum(dist.size.death*meshpts)*h,4)

#2nd the variance is

round(sum(dist.size.death*meshpts*meshpts*h)-sum(dist.size.death*meshpts*h)*sum(dist.size.death*meshpts*h),4)

#check with the simulation

round(with(sim.data,var(z[alive==0])),4)

#Reproduction: who, when and how much? 

#As reproduction is fatal the P kernel is the required P0 kernel

P0 <- P

N0 <- solve(diag(100)-P0)

#So B is given by

B <- p_bz(meshpts,m.par.true) %*% N0

# plot(meshpts,B,type="l")
# points(meshpts,p_bz(meshpts,m.par.true),type="l",col="red")

B.m.c <- matrix(B,nrow=100,ncol=100)

B.m.r <- matrix(B,nrow=100,ncol=100,byrow=TRUE)

P.b <- (P0 * B.m.c ) / B.m.r

N.b <- solve(diag(100)-P.b)

mean.Repr.age <- sum((e %*% N.b) * offspring.dist * h)-1

mean.Repr.age
with(sim.data,mean(age[Repr == 1]))

#these don't match and I don't know why it must be a fault with the simulation, as the next calculation seems fine

Omega.b <- as.vector(1-(e %*% P.b)) * N.b

#Distribution of sizes at reproduction

dist.size.repr <- (Omega.b %*% offspring.dist)

#mean size at reproduction

mean.size.flowering <- sum(h*dist.size.repr*meshpts)

with(sim.data,mean(z[Repr == 1]))

#variance in size at reproduction

sum(dist.size.repr*meshpts*meshpts)-sum(dist.size.repr*meshpts)*sum(dist.size.repr*meshpts)

with(sim.data,var(z[Repr == 1]))

#How often do they reproduce?

p_bz(meshpts,m.par.true) %*% N0 / B

#well they're monocarpic and reproduction is fatal

sum((e %*% F %*% N) * offspring.dist)

R0=abs(Re(eigen(F %*% solve(diag(100)-P))$values[1]))

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















