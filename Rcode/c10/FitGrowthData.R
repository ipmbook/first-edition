rm(list=ls(all=TRUE))
## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 

require(stats4); 
source("../Utilities/Standard Graphical Pars.R"); 

X=read.csv("CRLNAGRW.csv");
size1 <- X$lst; size2 <- X$lst1; #log-transformed sizes 

####################################################### 
## plot the data 
#######################################################
graphics.off; set_graph_pars("panel4"); 

plot(size1,size2,xlab="log Size t", ylab="log Size t+1")
add_panel_label("a"); 

####################################################### 
### Fit simple linear regression by maximum likelihood 
### size2 ~ Gaussian(mean=a+b*size1, sd=sigma)
####################################################### 

LinearLogLik1 <- function(a,b,sigma) {
	loglik <- sum(dnorm(size2,mean=a+b*size1,sd=sigma,log=TRUE))
	return(-loglik)
}
fit <- mle(LinearLogLik1,start=list(a=1,b=1,sigma=1), method="Nelder-Mead")	
summary(fit);
parms <- coef(fit); a<-parms[1]; b<-parms[2]; sigma<-parms[3]; 
abline(a,b,lwd=2); 

## show the estimated standard deviation on the plot 
s1<-seq(min(size1),1.1*max(size1),length=10); 
up<-a+b*s1+2*sigma; down<-a+b*s1-2*sigma; 
matpoints(s1,cbind(down,up),type="l",lty=2,lwd=1,col="black"); 

####################################################### 
## Scale-location plot
####################################################### 
# raw residuals 
err<-size2-(a+b*size1); 
# Standardize the residuals 
Xmat=cbind(rep(1,length(size1)),size1); # design matrix 
H=Xmat%*%solve(t(Xmat)%*%Xmat)%*%t(Xmat); 
sresid = err/sqrt(1-diag(H)); 

plot(size1,sqrt(abs(sresid)),xlab="log Size t", ylab="Sqrt |Std Residuals|",
type="p",cex=1.2);  
rfit <- lm(sqrt(abs(sresid))~size1); abline(rfit,lwd=2,lty=1)

add_panel_label("b"); 

####################################################### 
### Fit a linear model with nonconstant variance 
####################################################### 
# likelihood function 
LinearLogLik2=function(a,b,sigma0,d) {
	loglik=sum(dnorm(size2,mean=a+b*size1,sd=sigma0*exp(-d*size1),log=TRUE))
	return(-loglik)
}
# do the fit 	
fit2=mle(LinearLogLik2,start=list(a=1,b=1,sigma0=sqrt(mean(err^2)),d=0), method="Nelder-Mead")	
summary(fit2); confint(fit2); 

# Likelihood ratio test for constant vs. nonconstant variance 
2*(logLik(fit2)-logLik(fit)); 

## Plot the fit and the regression line 
plot(size1,size2,xlab="log Size t", ylab="log Size t+1")

## add the regression line to the plot
parms=coef(fit2); 
a=parms[1]; b=parms[2]; sigma0=parms[3]; d=parms[4];
abline(a,b,lty=1,lwd=2); 

## add the estimated standard deviation
s1=seq(min(size1),1.1*max(size1),length=50); 
up=a+b*s1+2*(sigma0*exp(-d*s1)); 
down=a+b*s1-2*(sigma0*exp(-d*s1)); 
matpoints(s1,cbind(down,up),type="l",lty=2,lwd=1,col="black"); 
add_panel_label("c"); 

####################################################### 
## Scale-location plot for nonconstant variance model
####################################################### 
err=size2-(a+b*size1); 
err=err/(sigma0*exp(-d*size1));  

# Standardize the residuals 
Winv <- diag((sigma0*exp(-d*size1))^(-2)); 
Xmat=cbind(rep(1,length(size1)),size1); 
H=Xmat%*%solve(t(Xmat)%*%Winv%*%Xmat)%*%t(Xmat)%*%Winv; 
stdresid = err/sqrt(1-diag(H)); 

plot(size1,sqrt(abs(stdresid)),xlab="log Size t",ylab="Sqrt |Std Residuals|"); 
rfit2 <- lm(sqrt(abs(stdresid))~size1); abline(rfit2,lwd=2,lty=1)

add_panel_label("d"); 
dev.copy2eps(file="../../c10/figures/FitGrowthData.eps"); 