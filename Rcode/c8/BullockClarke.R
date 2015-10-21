# Analysis of Calluna north-east transect, using data tabulated
# in Bullock and Clarke 2000. 

source("../utilities/Standard Graphical Pars.R")

# The data 
distance <- c(0.6,0.8,1,1.5,2,3,4,6,8,10,15,20,30,40,60,80);
seedNumber <- c(2004,1323,369,149,86,69,27,17,6,11,3,2,0,0,0,0);
trapNumber <- c(2,2,2,3,4,6,8,12,16,20,30,40,60,80,120,160);

# reality check: direct calculation here vs. numbers in their Table
# for the estimated number of seeds per m^2 as a function of distance
trapSize <- pi*(0.09/2)^2
trapArea <- trapNumber*trapSize;
seedDensity <- seedNumber/trapArea;

# check that our calculation matches what's in the paper 
seedDensity1 <- c(157504,103981,29002,7807,3380,1808,531,223,59,86,16,8,0,0,0,0); 
plot(seedDensity,seedDensity1); abline(0,1)  #Good!

dev.new(); 
set_graph_pars("panel1"); par(cex.lab=1.4,cex.axis=1.3,pch=19)
plot(log(distance[1:12]),log(seedDensity[1:12]),xlab="Log distance from shrub center",
	 ylab="Log seeds/m2");
fit0 <- lm(log(seedDensity[1:12])~log(distance[1:12])); abline(fit0);


nllPower <- function(a,b) {
	mean=a*(distance^(-b))*trapArea
	d = dpois(seedNumber,mean,log=TRUE);
	return(sum(-d))
}
require(stats4);
a0 <- exp(10.4); b0 <-2.8; 
fit <- mle(nllPower,start=list(a=a0,b=b0),method="Nelder-Mead",control=list(maxit=5000,trace=4))

dev.copy2eps(file="../../c8/figures/BullockClarke.eps")