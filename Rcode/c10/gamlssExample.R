require(gamlss);
require(gamlss.add);

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 

graphics.off();

## run the utility functions
source("~/repos/ipm_book/Rcode/utilities/Standard Graphical Pars.R")

set_graph_pars(ptype = "panel4")
# par(cex.axis=1.25,cex.lab=1.4); 

# Generate data with t distribution 
z <- seq(1,50,length=150); z <- sort(z); 
z1bar <- 40*z/(15+z);  
z1=z1bar+rt(150,df=5)*4*exp(-0.5+z/50)

plot(z,z1,xlab="Size at time t",
ylab="Size at time t+1"); 
points(z,z1bar,col="black",type="l",lwd=2); 

add_panel_label("a")

# Fit mean as penalized regression spline with 3rd derivative penalty
# scale parameter as penalized regression spline with 2nd deriv penalty (default)
# and df constant
fit1=gamlss(z1~pb(z,control=pb.control(degree=4,order=3)),sigma.fo=~pb(z),nu.fo=~1,family=TF);

# Same except scale parameter is given a parametric model 
# Default link function for scale parameter is the log, so the 
# fit below specifies that log(scale parameter)=linear function of z,
# which is the correct form for the data 
fit2=gamlss(z1~pb(z,control=pb.control(degree=4,order=3)),sigma.fo=~z,nu.fo=~1,family=TF);

# Same as the first fit, but the penalized regression splines are replaced
# with thin plate regression splines using ga() which provides an interface
# from gamlss to the smoothers in mgcv. 
# This is not plotted, but you can verify that predict(fit3) is almost identical
# to predict(fit1) and predict(fit2)
X=data.frame(z=z,z1=z1)
fit3=gamlss(z1~ga(~s(z,m=3)),sigma.fo=~ga(~s(z)),nu.fo=~1,family=TF,data=X);

# matpoints(z,cbind(predict(fit1),predict(fit2)),type="l",
# col=c("red","blue"), lty=c(1,2,3),lwd=c(2,3,4)); 

Runit=FALSE
if(Runit) {
nReps=250; 
yhat1 <- yhat2 <- sighat1 <- sighat2<- matrix(0,150,nReps) 
dfhat1 <- dfhat2 <- numeric(nReps);
zj <- seq(1,50,length=150);
for(j in 1:nReps) {
	z1bar <- 40*zj/(15+zj);  
	z1=z1bar+rt(150,df=5)*4*exp(-0.5+zj/50)
	fit1j=gamlss(z1~pb(zj,control=pb.control(degree=4,order=3)),sigma.fo=~pb(zj),nu.fo=~1,family=TF);
	fit2j=gamlss(z1~pb(zj,control=pb.control(degree=4,order=3)),sigma.fo=~zj,nu.fo=~1,family=TF);
	sighat1[,j]=fit1j$sigma.fv; sighat2[,j]=fit2j$sigma.fv; 
	yhat1[,j]=predict(fit1j); yhat2[,j]=predict(fit2j); 
	dfhat1[j]=fit1j$nu.fv[1]; 
	dfhat2[j]=fit2j$nu.fv[1]; 
	cat(j,dfhat1[j],dfhat2[j],"\n")
	
}

}# end if Runit

# plot pointwise 5th and 95th percentiles of the plot
r90 <- function(x) quantile(x,probs=c(0.05,0.95))
yr1 <- apply(yhat1,1,r90)
yr2 <- apply(yhat2,1,r90)

matpoints(z,cbind(t(yr1),t(yr2)),type="l",lwd=c(2,3),
col=c("red","red","blue","blue"), lty=c(2,2,3,3)); 

matplot(z,4*exp(-0.5+z/50),
type="l",lwd=2,col=c("black","red","blue"),
lty=c(1,2,3),xlab="Individual size", ylab="Growth standard deviation");

legend("topleft",legend=c("pb(z) + pb(z)","pb(z) + z"),col=c("red","blue"),
bty="n",lty=c(2,3),lwd=c(2,3),cex=1.2)
add_panel_label("b")

ys1 <- apply(sighat1,1,r90)
ys2 <- apply(sighat2,1,r90)

matpoints(z,cbind(t(ys1),t(ys2)),type="l",lwd=2,
col=c("red","red","blue","blue"), lty=c(2,2,3,3)); 

e<-((dfhat2<20)&(dfhat1<20))
plot(dfhat2[e],dfhat1[e],xlab="Estimated df, true sd function",
ylab="Estimated df, spline sd function");
add_panel_label("c")

dev.copy2eps(file="~/Repos/ipm_book/HowTo/gamlssExample.eps")






