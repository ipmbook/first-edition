### This script assumes that you have just source'd the Monocarp model
### using MonocarpSimulateIBM.R 
  
# or, load an .Rdata file with saved simulation results 
load("MonocarpSimData.Rdata")

require(car)
require(mgcv)

source("../utilities/Standard Graphical Pars.R")

## Construct a data set of plausible size
pick.data <- seq(1, nrow(sim.data), length = 300)
test.data <- sim.data[round(pick.data), ]
test.data <- na.omit(subset(test.data, select = c(size, size1)))
e <- order(test.data$size)
test.data <- test.data[e, ]

## refit models to the reduced data set
mod.grow <- lm(size1 ~ size, data = test.data)
cat(length(mod.grow$fitted))

set_graph_pars("panel4")

# Plot residuals versus fitted for growth model
zhat <- fitted(mod.grow)
resid <- residuals(mod.grow)
plot(zhat, resid, xlab = "Fitted values", ylab = "Residuals")
gam.resid <- gam(resid ~ s(zhat), method = "REML")
rhat <- predict(gam.resid, type = "response")
points(zhat, rhat, type = "l")
add_panel_label("a")

# Normal qq-plot for growth model
sresid <- rstandard(mod.grow)
qqPlot(sresid, main = "", xlab = "Normal quantiles", ylab = "Standardized residual quantiles", 
       col.lines = "black", lwd = 1)
add_panel_label("b")

# Absolute residuals versus fitted
plot(zhat, sqrt(abs(sresid)), xlab = "Fitted values", ylab = "sqrt(|Std Residuals|)")
gam.sresid <- gam(sqrt(abs(sresid)) ~ s(zhat), method = "REML")
rhat <- predict(gam.sresid, type = "response")
points(zhat, rhat, type = "l")
add_panel_label("c")

# compare to a gam fit
gam.grow <- gam(size1 ~ s(size), data = test.data, method = "REML")
AIC(gam.grow, mod.grow)
gam.grow.fitted <- predict(gam.grow, type = "response")
matplot(test.data$size, cbind(fitted(mod.grow), gam.grow.fitted), type = "l", 
        lty = c(1, 2), lwd = 2, xlab = "Size t", ylab = "Fitted size t+1")
add_panel_label("d")

# dev.copy2eps(file = "../../figures/c2/DiagnoseMonocarp1.eps")
