# Profile negative log likelihood function for transformation parameter
bcNLL=function(lambda,x,y) {
	xl=x^lambda; yl=y^lambda;
	fit=lm(yl~xl); s2hat=mean(fit$residuals^2); 
	Lmax = -0.5*length(x)*log(s2hat/lambda^2) + (lambda-1)*sum(log(y));
	return(-Lmax)
}

# Do replicate trials to see how well lambda is estimated
# In this test case the true value is lambda=0.5
lvals=numeric(5000); N=100;
for(j in 1:5000) {
  # Data on the scale where a linear model is valid
  z=runif(N,1,10); z1=1 + 0.5*z + 0.2*rnorm(N); 

  # the 'data', created so square-root transformation is optimal 
  x=z^2; x1=z1^2; # plot(x,x1);

  lvals[j]=optimize(bcNLL,c(0.01,3),x=x,y=x1)$minimum 
  if(j%%250==0) cat(j,"\n"); 
}
 


