## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Monocarp IBM to illustrate the construction of an IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))

library(doBy)
set.seed(53241986)

# Set simulation parameters
init.pop.size <- 250
n.yrs <- 100

source("../utilities/Standard Graphical Pars.R")

source("Monocarp Demog Funs.R")
source("Monocarp Simulate IBM.R")
cat(pop.size.t, "\n")

## trim an initial transient off the simulation
sim.data <- sim.data[sim.data$yr > 10, ]
sim.data$yr <- sim.data$yr - 10

# Extract 1000 observations to use as our data set

sample.index <- sample(1:nrow(sim.data), size = 1000, replace = FALSE)
sim.data <- sim.data[sample.index, ]

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Fit statistical models to simulated data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Let's use sim data to construct an IPM, assuming we know the order of
## events first a bit of house keeping for plotting - chop up the sizes
## into classes and sort the data by size

sim.data <- transform(sim.data, z.classes = cut(z, 20))
sim.data <- sim.data[order(sim.data$z), ]

## fit the functions and plot some graphs

dev.new(); 
set_graph_pars("panel4")

## 1 - growth

plot(z1 ~ z, data = sim.data, pch = 16, cex = 0.25, 
     xlab = expression("Size t, " * italic(z)), ylab = expression("Size t+1, " * italic(z) * "'"))

mod.Grow <- lm(z1 ~ z, data = sim.data)

abline(mod.Grow, col = "red")

add_panel_label("a")

## 2 - flowering

Repr.ps <- summaryBy(z + Repr ~ z.classes, data = sim.data)

plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlab = expression("Size t, " * 
                                                                                  italic(z)), ylab = "Probability of flowering")

mod.Repr <- glm(Repr ~ z, family = binomial, data = sim.data)

lines(fitted(mod.Repr) ~ z, data = sim.data, col = "red")

add_panel_label("b")

## 3 - survival

sim.data.noRepr <- subset(sim.data, Repr == 0)

surv.ps <- summaryBy(z + Surv ~ z.classes, data = sim.data.noRepr, na.rm = TRUE)

plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19, xlab = expression("Size t, " * 
                                                                       italic(z)), ylab = "Probability of survival")

mod.Surv <- glm(Surv ~ z, family = binomial, data = sim.data.noRepr)

lines(fitted(mod.Surv) ~ z, data = sim.data.noRepr, col = "red")

add_panel_label("c")

## 4 - seed production

sim.data.Repr <- subset(sim.data, Repr == 1)

plot(log(Seeds) ~ z, data = sim.data.Repr, pch = 19, xlab = expression("Size t, " * 
                                                                         italic(z)), ylab = "Seeds production")

mod.Seeds <- glm(Seeds ~ z, family = poisson, data = sim.data.Repr)

abline(mod.Seeds, col = "red")

add_panel_label("d")

# dev.copy2eps(file='~/Repos/ipm_book/c2/figures/OenotheraDemog.eps');

## 5 - recruit size

sim.data.Rec <- subset(sim.data, age == 0)

mod.Rcsz <- lm(z ~ 1, sim.data.Rec)


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
## Section 2 - Construct Kernels and projection population size
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nBigMatrix <- 250
min.size <- with(sim.data, min(z))
min.size
max.size <- with(sim.data, max(z))
max.size

# so let's set the lower and upper limits for the size range at -2.65
# and 4.5 so slightly smaller/bigger than observed

IPM.true <- mk_K(nBigMatrix, m.par.true, -2.65, 4.5)

IPM.est <- mk_K(nBigMatrix, m.par.est, -2.65, 4.5)

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

dev.new(); 
set_graph_pars("panel4")

## 1 - plot population density versus time...

plot(1:yr, log(pop.size.t[1:yr]), type = "l", xlab = "Time", ylab = "Population size")
mod.pop.growth <- lm(log(pop.size.t) ~ seq.int(1, yr))
abline(mod.pop.growth, col = "blue")

add_panel_label("a")

## ...roughly linear for log(Nt) vs time so exponential growth

## 2 - plot mean size versus time...
plot(1:yr, mean.z.t[1:yr], type = "l", xlab = "Time", ylab = "Mean plant size")
abline(h = mean.z.true, col = "red")
add_panel_label("b")

## ...mean size seems to settle down after a while

## 3 - plot mean flowering size versus time...
plot(1:yr, mean.fl.z.t[1:yr], type = "l", xlab = "Time", ylab = "Mean flowering plant size")
abline(h = mean.flowering.z.true, col = "red")
add_panel_label("c")

## ...mean flowering size seems to settle down after a while

## 4 - plot of density estimates at time 50 and the end
plot(density(sim.data$z), main = "", ylim = c(0, 0.4), xlab = "Plant size")
# lines(density(sim.data$z))
lines(IPM.est$meshpts, stable.z.dist.est/diff(IPM.est$meshpts)[1], col = "red")
add_panel_label("d")

# dev.copy2eps(file = "../../figures/c2/OenotheraSim.eps")
