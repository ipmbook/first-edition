### Assumes that you have just source'd the Monocarp model
### or loaded the results saved in MonocarpSimData.Rdata
require(car); require(mgcv); require(doBy); 
load("~/Repos/ipm_book/Rcode/MonocarpSimData.Rdata")

## Construct a data set of plausible size
pick.data <- seq(1,nrow(sim.data),length=200); 
test.data <- sim.data[round(pick.data),];
test.data<- subset(test.data, flow==0)
e <- order(test.data$size); test.data <- test.data[e,]; 
attach(test.data); 

# fit models to the reduced data set 
mod.surv <- glm(surv ~ size, family = binomial, data = test.data)
gam.surv <- gam(surv ~ s(size), family = binomial, data = test.data,method="REML")

dev.new(height=4,width=8);
set_graph_pars("panel2"); 
#graphics.off(); quartz(w=6.5,h=3.5); 
#par(mfrow=c(1,2),mgp=c(2.5,1,0),bty="l",yaxs="i",mar=c(4,4,3,1),cex.main=1.1); 

# compare fitted glm, grouped data, and fitted gam 
sizerank <- seq_along(size)/nrow(test.data); 
test.data$sizeclass <- round(9*sizerank)+1; # make size classes based on quantiles
surv.ps <- summaryBy(size + surv ~ sizeclass, data = test.data, na.rm = TRUE)
attach(surv.ps); 
plot(size.mean,surv.mean,xlab="Size z", ylab="Survival probability",pch=1,
     xlim=range(size),ylim=c(0,1));
#,main="GLM and GAM models"
points(size,fitted(mod.surv),type="l",lty=1,lwd=2); 
svals <- seq(min(size),max(size),length=12); 
ghat <- predict(gam.surv,newdata=data.frame(size=svals),type="response"); 
points(svals,ghat,type="p",col="red",pch=16) 

add_panel_label("a")
#mtext("A)",side=3,adj=0,cex=1.3);


# Show the fitted GAM's linear predictor
plot(gam.surv,seWithMean=TRUE,xlab="Size z",ylab="Spline(z,edf=1.09)"); 
#title(main="GAM linear predictor"); 

add_panel_label("b")
#mtext("B)",side=3,adj=0,cex=1.3);

dev.copy2eps(file="~/Repos/ipm_book/c2/figures/DiagnoseMonocarp2.eps");

