# Negative log likelihood function for regression with mean and standard
# deviation both linear functions of the independent variable
require(stats4)
ncvNLL=function(a,b,c,d) {
	LL <- sum(dnorm(y,mean=a+b*x,sd=c+d*x,log=TRUE)); 
	return(-LL)
}

N=100;
x=runif(N,1,10); y=1 + 0.5*x + (0.1+0.02*x)*rnorm(N); 
fit=mle(ncvNLL,start=list(a=2,b=1,c=1,d=1)); 
summary(fit); confint(fit);


