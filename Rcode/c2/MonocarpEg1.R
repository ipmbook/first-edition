## IBM to generate the data for the simple recruitment limited IPM -
## there's no recruitment limitation 1st Size is on a log scale so don't
## worry about the -ive ones!

## The parameters and life cycle follow Kachi and Hirose J.Ecol. 1985,
## see Rees and Rose PRSB 2002, but for simplicity we have used a
## logistic function for survival rather than the linear function with
## an upper bound as in the previous papers

rm(list = ls(all = TRUE))

library(doBy)
set.seed(53241986)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so
## the formulae, below, are easier to read. We'll use 'model.term'
## scheme to name elements of the vector

m.par.true <- c(surv.int = -0.65, surv.z = 0.75, flow.int = -18, flow.z = 6.9, 
                grow.int = 0.96, grow.z = 0.59, grow.sd = 0.67, rcsz.int = -0.08, rcsz.sd = 0.76, 
                seed.int = 1, seed.z = 2.2, p.r = 0.007)

## Define the various demographic functions, we pass a vector of
## parameters 'm.par' to each function this makes it easier to compare
## the results of different paramter sets, say the true values and
## estimated ones

## Growth function, given you are size z now returns the pdf of size z1
## next time

g_z1z <- function(z1, z, m.par) {
  mean <- m.par["grow.int"] + m.par["grow.z"] * z  # mean size next year
  sd <- m.par["grow.sd"]  # sd about mean
  p.den.grow <- dnorm(z1, mean = mean, sd = sd)  # pdf that you are size z1 given you were size z
  return(p.den.grow)
}

## Survival function, logistic regression

s_z <- function(z, m.par) {
  linear.p <- m.par["surv.int"] + m.par["surv.z"] * z  # linear predictor
  p <- 1/(1 + exp(-linear.p))  # logistic transformation to probability
  return(p)
}

## Probability of flowering function, logistic regression

p_bz <- function(z, m.par) {
  linear.p <- m.par["flow.int"] + m.par["flow.z"] * z  # linear predictor
  p <- 1/(1 + exp(-linear.p))  # logistic transformation to probability
  return(p)
}

## Seed production function

b_z <- function(z, m.par) {
  N <- exp(m.par["seed.int"] + m.par["seed.z"] * z)  # seed production of a size z plant
  return(N)
}

## Recruit size pdf

c_z1 <- function(z1, m.par) {
  mean <- m.par["rcsz.int"]
  sd <- m.par["rcsz.sd"]
  p.deRecr <- dnorm(z1, mean = mean, sd = sd)  # pdf of a size z1 recruit
  return(p.deRecr)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Parameters controlling the IBM simulation

m.par <- m.par.true

init.pop.size <- 250
n.yrs <- 50
z <- rnorm(init.pop.size, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])

## Calculate initial pop size and mean size

pop.size.t <- init.pop.size
mean.z.t <- mean(z)
mean.fl.z.t <- NA

## iterate the model using the 'true' parameters and store data in a
## data.frame
yr <- 1
while (yr != n.yrs & length(z) < 5000) {
  
  ## calculate population size
  pop.size <- length(z)
  
  ## generate binomial random number for the probability of flowering,
  ## where the probability of flowering depends on your size z, this is a
  ## vector of 0's and 1's, you get a 1 if you flower
  Repr <- rbinom(n = pop.size, prob = p_bz(z, m.par), size = 1)
  
  ## number of plants that flowered
  num.Repr <- sum(Repr)
  
  ## calculate seed production
  Seeds <- rep(NA, pop.size)
  
  ## we'll assume plant make a Poisson distributed number of seeds with a
  ## mean given by exp(params['seed.int']+params['seed.size'] * z) rpois
  ## generated Poisson distributed random numbers
  Seeds[Repr == 1] <- rpois(num.Repr, b_z(z[Repr == 1], m.par))
  
  ## generate the number of recruits
  Recr <- if (num.Repr == 0) 
    0 else rbinom(1, sum(Seeds, na.rm = TRUE), m.par["p.r"])
  
  ## generate new recruit sizes rnorm generated normally distributed
  ## random numbers
  Rcsz <- rnorm(Recr, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])
  
  ## for the non-reproductive plants generate random number for survival
  Surv <- rep(NA, pop.size)
  Surv[Repr == 0] <- rbinom(n = pop.size - num.Repr, prob = s_z(z[Repr == 
                                                                    0], m.par), size = 1)
  num.die <- sum(Surv == 0, na.rm = TRUE)
  
  ## index for individuals that did not flower and survived
  i.subset <- which(Repr == 0 & Surv == 1)
  
  ## let them grow
  z1 <- rep(NA, pop.size)
  E.z1 <- m.par["grow.int"] + m.par["grow.z"] * z[i.subset]
  z1[i.subset] <- rnorm(n = pop.size - num.Repr - num.die, mean = E.z1, 
                        sd = m.par["grow.sd"])
  
  ## store the simulation data, we'll use this later
  sim.data <- data.frame(z = z, Repr, Seeds, Surv, z1 = z1)
  
  ## create new population size vector
  z1 <- c(Rcsz, z1[i.subset])
  
  pop.size.t <- c(pop.size.t, length(z1))
  mean.z.t <- c(mean.z.t, mean(z1))
  mean.fl.z.t <- c(mean.fl.z.t, mean(z[Repr == 1]))
  
  z <- z1
  
  if (yr == n.yrs) 
    z.time50 <- z.sim
  
  cat(paste(yr, pop.size.t[yr], mean.z.t[yr], Recr/sum(Seeds, na.rm = TRUE), 
            "\n", sep = " "))
  
  yr <- yr + 1
  
}

# Extract 500 observations to use as our data set

sample.index <- sample(1:nrow(sim.data), size = 1000, replace = FALSE)
sim.data <- sim.data[sample.index, ]

### CHANGED DOWN TO HERE, beyond this won't work.

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 3 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Let's use sim data to construct an IPM, assuming we know the order of
## events first a bit of house keeping for plotting - chop up the sizes
## into classes and sort the data by size

sim.data <- transform(sim.data, z.classes = cut(z, 20))
sim.data <- sim.data[order(sim.data$z), ]

## fit the functions and plot some graphs

dev.new(6, 6)

par(mfrow = c(2, 2), bty = "l", pty = "s", pch = 19)

## 1 - growth

plot(z1 ~ z, data = sim.data, pch = 16, cex = 0.25, 
     xlab = expression("Size t, " * italic(z)), ylab = expression("Size t+1, " * italic(z) * "'"))

mod.Grow <- lm(z1 ~ z, data = sim.data)

abline(mod.Grow, col = "red")

mtext(side = 3, line = 0.5, adj = 0, text = "a)")


## 2 - flowering

Repr.ps <- summaryBy(z + Repr ~ z.classes, data = sim.data)

plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, 
     xlab = expression("Size t, " * italic(z)), ylab = "Probability of flowering")

mod.Repr <- glm(Repr ~ z, family = binomial, data = sim.data)

lines(fitted(mod.Repr) ~ z, data = sim.data, col = "red")

mtext(side = 3, line = 0.5, adj = 0, text = "b)")

## 3 - survival

sim.data.noRepr <- subset(sim.data, Repr == 0)

surv.ps <- summaryBy(z + Surv ~ z.classes, data = sim.data.noRepr, na.rm = TRUE)

plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19, xlab = expression("Size t, " * 
                                                                       italic(z)), ylab = "Probability of survival")

mod.Surv <- glm(Surv ~ z, family = binomial, data = sim.data.noRepr)

lines(fitted(mod.Surv) ~ z, data = sim.data.noRepr, col = "red")

mtext(side = 3, line = 0.5, adj = 0, text = "c)")

## 4 - seed production

sim.data.Repr <- subset(sim.data, Repr == 1)

plot(log(Seeds) ~ z, data = sim.data.Repr, xlab = expression("Size t, " * 
                                                               italic(z)), ylab = "Seeds production")

mod.Seeds <- glm(Seeds ~ z, family = poisson, data = sim.data.Repr)

abline(mod.Seeds, col = "red")

mtext(side = 3, line = 0.5, adj = 0, text = "d)")

dev.copy2eps(file='../../figures/c2/OenotheraDemog.eps');

## 5 - recruit size

mod.Rcsz <- lm(Rcsz ~ 1)

## 6 - establishment probability WE SHOULD USE sim.data TO DO THIS CALC

p.r.est <- as.numeric(Recr)/sum(Seeds, na.rm = TRUE)

## Finally, store the estimated parameters

m.par.est <- c(surv = coef(mod.Surv), flow = coef(mod.Repr), grow = coef(mod.Grow), 
               grow.sd = summary(mod.Grow)$sigma, rcsz = coef(mod.Rcsz), rcsz.sd = summary(mod.Rcsz)$sigma, 
               seed = coef(mod.Seeds), p.r = p.r.est)

# par.names <- names(m.par.est) names(m.par.est) <- sub('(Intercept)',
# 'int', par.names, fixed=TRUE)

names(m.par.est) <- names(m.par.true)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 4 - Build and IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
p_z1z <- function(z1, z, m.par) {
  
  return((1 - p_bz(z, m.par)) * s_z(z, m.par) * g_z1z(z1, z, m.par))
  
}

## Define the fecundity kernel
f_z1z <- function(z1, z, m.par) {
  
  return(p_bz(z, m.par) * b_z(z, m.par) * m.par["p.r"] * c_z1(z1, m.par))
  
}

mk_K <- function(m, m.par) {
  
  # upper and lower integration limits
  L <- 0.9 * min.size
  U <- 1.1 * max.size
  
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

IPM.est <- mk_K(100, m.par.est)

Re(eigen(IPM.true$K)$values[1])

Re(eigen(IPM.est$K)$values[1])

fit.pop.growth <- lm(log(pop.size.t) ~ c(1:yr))

exp(coef(fit.pop.growth)[2])

meshpts <- IPM.true$meshpts

w.est <- Re(eigen(IPM.est$K)$vectors[, 1])
stable.z.dist.est <- w.est/sum(w.est)
mean.z.est <- sum(stable.z.dist.est * meshpts)
mean.z.est

w.true <- Re(eigen(IPM.true$K)$vectors[, 1])
stable.z.dist.true <- w.true/sum(w.true)
mean.z.true <- sum(stable.z.dist.true * meshpts)
mean.z.true

wb.est <- p_bz(meshpts, m.par.est) * w.est
stable.flowering.dist.est <- wb.est/sum(wb.est)
mean.flowering.z.est <- sum(stable.flowering.dist.est * meshpts)
mean.flowering.z.est


wb.true <- p_bz(meshpts, m.par.true) * w.true
stable.flowering.dist.true <- wb.true/sum(wb.true)
mean.flowering.z.true <- sum(stable.flowering.dist.true * meshpts)
mean.flowering.z.true

par(mfrow = c(2, 2), bty = "l", pty = "s")

## 1 - plot population density versus time...

plot(1:yr, log(pop.size.t[1:yr]), type = "l", xlab = "Time", ylab = "Population size")
mod.pop.growth <- lm(log(pop.size.t) ~ seq.int(1, yr))
abline(mod.pop.growth, col = "blue")

mtext(side = 3, line = 0.5, adj = 0, text = "a)")
## ...roughly linear for log(Nt) vs time so exponential growth


## 2 - plot mean size versus time...
plot(1:yr, mean.z.t[1:yr], type = "l", xlab = "Time", ylab = "Mean plant size")
abline(h = mean.z.true, col = "red")
mtext(side = 3, line = 0.5, adj = 0, text = "b)")
## ...mean size seems to settle down after a while

## 3 - plot mean flowering size versus time...
plot(1:yr, mean.fl.z.t[1:yr], type = "l", xlab = "Time", ylab = "Mean flowering plant size")
abline(h = mean.flowering.z.true, col = "red")
mtext(side = 3, line = 0.5, adj = 0, text = "c)")
## ...mean flowering size seems to settle down after a while

## 4 - plot of density estimates at time 50 and the end
plot(density(sim.data$z), main = "", ylim = c(0, 0.4))
# lines(density(sim.data$z))
lines(IPM.est$meshpts, stable.z.dist.est/diff(IPM.est$meshpts)[1], col = "red")
mtext(side = 3, line = 0.5, adj = 0, text = "d)")

# dev.copy2eps(file='../../figures/c2/OenotheraSim.eps');

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 5 - Example calculations
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

