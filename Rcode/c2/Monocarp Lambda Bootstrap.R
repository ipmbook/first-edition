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
set.seed(53241986)

source("Monocarp Demog Funs.R")
source("../utilities/Standard Graphical Pars.R")

# General IPM parameters
nBigMatrix <- 250
p.r.est <- p.r <- 0.007

IPM.est <- mk_K(nBigMatrix, m.par.true, -2.65, 4.5)
true.lambda <- Re(eigen(IPM.est$K)$values[1])

# Simulate IBM to get very large data set.
init.pop.size <- 1000
n.yrs <- 75
source("Monocarp Simulate IBM.R")

# Drop the first 10 years of transients
sim.data.big <- sim.data[sim.data$yr > 10, ]
sim.data.big$yr <- sim.data.big$yr - 10

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 1: the unknowable 'TRUTH': How much does estimated lamba vary
# across many replicate data sets?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No.sims <- 250
Lambdas <- rep(NA, No.sims)

indices <- sample(1:(No.sims * 1000), replace = FALSE)
indices <- matrix(indices, ncol = No.sims)

for (i in 1:No.sims) {
  
  # Extract 1000 observations to use as our data set
  
  sample.index <- indices[, i]
  sim.data <- sim.data.big[sample.index, ]
  
  ## Fit statistical models to simulated data
  
  mod.Grow <- lm(z1 ~ z, data = sim.data)
  mod.Repr <- glm(Repr ~ z, family = binomial, data = sim.data)
  sim.data.noRepr <- subset(sim.data, Repr == 0)
  mod.Surv <- glm(Surv ~ z, family = binomial, data = sim.data.noRepr)
  no.surv <- sum(sim.data.noRepr$Surv)
  sim.data.Repr <- subset(sim.data, Repr == 1)
  mod.Seeds <- glm(Seeds ~ z, family = poisson, data = sim.data.Repr)
  no.repr <- sum(sim.data.Repr$Repr)
  sim.data.Rec <- subset(sim.data, age == 0)
  mod.Rcsz <- lm(z ~ 1, sim.data.Rec)
  no.recs <- dim(sim.data.Rec)[1]
  
  ## Finally, store the estimated parameters
  m.par.est <- c(surv = coef(mod.Surv), flow = coef(mod.Repr), grow = coef(mod.Grow), 
                 grow.sd = summary(mod.Grow)$sigma, rcsz = coef(mod.Rcsz), rcsz.sd = summary(mod.Rcsz)$sigma, 
                 seed = coef(mod.Seeds), p.r = p.r.est)
  names(m.par.est) <- names(m.par.true)
  
  ## Construct Kernels and projection population size
  
  IPM.est <- mk_K(nBigMatrix, m.par.est, -2.65, 4.5)
  Lambdas[i] <- Re(eigen(IPM.est$K, only.values = TRUE)$values[1])
  cat(i, "   ", Lambdas[i], "  ", no.surv, "  ", no.repr, "  ", no.recs, "\n")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 2: What you can do with one real data set, is bootstrap by
# resampling from the data in hand.  We do this 5 times to see that the
# results are consistent
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~Function to do the bootstrap Inputs: a dataset, number
# of bootstrap samples, and value of parameter p.r
boot_lambda <- function(dataset, n.samp = 250, p.r) {
  
  lam.boot <- rep(NA, n.samp)
  for (i in 1:n.samp) {
    sample.index <- sample(1:nrow(dataset), size = dim(dataset)[1], 
                           replace = TRUE)
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
    lam.boot[i] <- Re(eigen(IPM.est$K, only.values = TRUE)$values[1])
    cat(i, "   ", lam.boot[i], "  ", no.surv, "  ", no.repr, "  ", no.recs, "\n")
  }
  return(lam.boot)
}

# Now we do 5 replicates of bootstrapping from one data set
n.pop.samples <- 20
for (pop in 1:n.pop.samples) {
  
  # Generate a data set
  init.pop.size <- 500
  n.yrs <- 75
  source("Monocarp Simulate IBM.R")
  sim.data <- sim.data[sim.data$yr > 10, ]
  sim.data$yr <- sim.data$yr - 10
  
  # Extract 1000 observations to use as our data set
  sample.index <- sample(1:nrow(sim.data), size = 1000, replace = FALSE)
  sim.data <- sim.data[sample.index, ]
  
  ## Bootstrap from the data set
  boot.strap <- boot_lambda(sim.data, n.samp = 250, p.r = p.r.est)
  
  if (pop == 1) {
    boot.data <- data.frame(boot.lam = boot.strap, pop = 1)
  } else {
    boot.data <- rbind(boot.data, data.frame(boot.lam = boot.strap, 
                                             pop = pop))
  }
}

# Plot the replicate data set results: distribution of lambda values
# across the replicate data sets
graphics.off()
dev.new(width = 8, height = 4)
set_graph_pars("panel2")
par(xpd = TRUE, yaxs = "i")
hist(Lambdas, xlim = c(0.8, 1.3), xlab = "Estimated lambda", ylab = "Frequency", 
     main = "", freq = FALSE, ylim = c(0, 7.75), col = "grey60", border = "grey35")
out <- density(Lambdas, bw = "SJ")
points(out$x, out$y, type = "l", lwd = 2)
points(x = rep(true.lambda, 2), y = c(0, 0), type = "p", pch = 16, lwd = 2, 
       cex = 1.5)
add_panel_label("a")


# ## Plot the bootstrap results. Histogram for the first out of the 5,
# ## and just the density for the other 4
# Lam=boot.data$boot.lam[boot.data$pop==1];
# hist(Lam,xlim=c(0.8,1.3),xlab='Estimated
# lambda',ylab='Frequency',main='',
# freq=FALSE,ylim=c(0,7.75),col='grey60',border='grey35');
# out=density(Lam,bw='SJ');points(out$x,out$y,type='l',lwd=2);
# add_panel_label('b')

# for(j in 2:5){ Lam=boot.data$boot.lam[boot.data$pop==j];;
# out=density(Lam,bw='SJ'); points(out$x,out$y,type='l',lty=2); }
# points(x=rep(true.lambda,2),y=c(0,0),type='p',pch=16,lwd=2,cex=1.5)
# dev.copy2eps(file='../../c2/figures/MonocarpLambdaBootstrap.eps')

## Plot the bootstrap results. Histogram for the first out of the 5, and
## just the density for the other 4

plot(out$x, out$y, type = "l", lwd = 2, xlab = "Estimated lambda", ylab = "Frequency", 
     ylim = c(0, 9.5))

add_panel_label("b")

for (j in 1:n.pop.samples) {
  Lam <- boot.data$boot.lam[boot.data$pop == j]
  out <- density(Lam, bw = "SJ")
  points(out$x, out$y, type = "l", lty = 2)
}
points(x = rep(true.lambda, 2), y = c(0, 0), type = "p", pch = 16, lwd = 2, 
       cex = 1.5)

# dev.copy2eps(file = "../../c2/figures/MonocarpLambdaBootstrap.eps")



