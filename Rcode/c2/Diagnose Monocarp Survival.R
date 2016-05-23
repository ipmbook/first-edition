### Assumes that you have just source'd the Monocarp model or loaded the
### results saved in MonocarpSimData.Rdata
require(car)
require(mgcv)
require(doBy)

load("MonocarpSimData.Rdata")
source("../utilities/Standard Graphical Pars.R"); 

## Construct a data set of plausible size
pick.data <- seq(1, nrow(sim.data), length = 200)
test.data <- sim.data[round(pick.data), ]
test.data <- subset(test.data, flow == 0)
e <- order(test.data$size)
test.data <- test.data[e, ]
attach(test.data)

# fit models to the reduced data set
mod.surv <- glm(surv ~ size, family = binomial, data = test.data)
gam.surv <- gam(surv ~ s(size), family = binomial, data = test.data, method = "REML")

dev.new(height = 4, width = 8)
set_graph_pars("panel2")

# compare fitted glm, grouped data, and fitted gam
sizerank <- seq_along(size)/nrow(test.data)
# make size classes based on quantiles
test.data$sizeclass <- round(9 * sizerank) + 1
surv.ps <- summaryBy(size + surv ~ sizeclass, data = test.data, na.rm = TRUE)
attach(surv.ps)
plot(size.mean, surv.mean, xlab = "Size z", ylab = "Survival probability", 
     pch = 1, xlim = range(size), ylim = c(0, 1))
points(size, fitted(mod.surv), type = "l", lty = 1, lwd = 2)
svals <- seq(min(size), max(size), length = 12)
ghat <- predict(gam.surv, newdata = data.frame(size = svals), type = "response")
points(svals, ghat, type = "p", col = "red", pch = 16)
add_panel_label("a")

# Show the fitted GAM's linear predictor
plot(gam.surv, seWithMean = TRUE, xlab = "Size z", ylab = "Spline(z,edf=1.09)")
add_panel_label("b")

# dev.copy2eps(file = "../../figures/c2/DiagnoseMonocarp2.eps")

