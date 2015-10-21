## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IBM to illustrate bootstrap uncertainty analysis for
## an IPM We assume the simplest situation, an unstructured bootstrap
## corresponding to data where new individuals are marked at random each
## year, and their fate the following year is observed, but individuals
## are not followed for life. We therefore construct bootstrap samples
## by resampling totally at random (with replacement) from the data.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
library(doBy)
library(boot)

set.seed(53241986)

source("Monocarp Demog Funs.R")
source("../utilities/Standard Graphical Pars.R")

## General IPM parameters
nBigMatrix <- 250
p.r.est <- p.r <- 0.007

IPM.est <- mk_K(nBigMatrix, m.par.true, -2.65, 4.5)
true.lambda <- Re(eigen(IPM.est$K)$values[1])

## Simulate IBM to get a data set.
init.pop.size <- 1000
n.yrs <- 50
source("Monocarp Simulate IBM.R")

## Use 49th year and select 1000 observations to be our dataset
sim.data <- sim.data[sim.data$yr == 49, ]
indices <- sample(1:nrow(sim.data), 1000, replace = FALSE)
sim.data <- sim.data[indices, ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to compute lambda from a bootstrapped data set in the format
# that the boot library requires.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boot.lam <- function(dataset, sample.index) {
  
  ## extract the data used to make this fit
  boot.data <- dataset[sample.index, ]
  
  ## fit the functions
  mod.Grow <- lm(z1 ~ z, data = boot.data)
  mod.Repr <- glm(Repr ~ z, family = binomial, data = boot.data)
  sim.data.noRepr <- subset(boot.data, Repr == 0)
  mod.Surv <- glm(Surv ~ z, family = binomial, data = sim.data.noRepr)
  no.surv <- sum(sim.data.noRepr$Surv)
  sim.data.Repr <- subset(boot.data, Repr == 1)
  mod.Seeds <- glm(Seeds ~ z, family = poisson, data = sim.data.Repr)
  no.repr <- sum(sim.data.Repr$Repr)
  sim.data.Rec <- subset(boot.data, age == 0)
  mod.Rcsz <- lm(z ~ 1, sim.data.Rec)
  no.recs <- dim(sim.data.Rec)[1]
  
  ## Store the estimated parameters
  m.par.est <- c(surv = coef(mod.Surv), flow = coef(mod.Repr), grow = coef(mod.Grow), 
                 grow.sd = summary(mod.Grow)$sigma, rcsz = coef(mod.Rcsz), rcsz.sd = summary(mod.Rcsz)$sigma, 
                 seed = coef(mod.Seeds), p.r = p.r)
  names(m.par.est) <- names(m.par.true)
  
  IPM.est <- mk_K(nBigMatrix, m.par.est, -2.65, 4.5)
  lam.boot <- Re(eigen(IPM.est$K, only.values = TRUE)$values[1])
  cat(lam.boot, "\n")
  return(lam.boot)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Do the bootstrap
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

boot.out <- boot(data = sim.data, statistic = boot.lam, simple = TRUE, 
                 R = 999)

boot.ci(boot.out, type = c("norm", "basic", "perc"))

#### Bias-corrected percentile intervals;
bcpi <- function(t0, t, alpha) {
  B <- length(t)
  z0 <- qnorm(mean(t < t0))
  a1 <- pnorm(2 * z0 + qnorm(alpha/2))
  a2 <- pnorm(2 * z0 + qnorm(1 - alpha/2))
  c1 <- quantile(t, a1)
  c2 <- quantile(t, a2)
  return(as.numeric(c(c1, c2)))
}
bcpi(boot.out$t0, boot.out$t, 0.05)
