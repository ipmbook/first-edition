## Neubert-Caswell deterministic MPM in 1 dimension with Gaussian kernel
require(Matrix); 
setwd("~/Repos/ipm_book/Rcode/c8"); 
## Define parameters 
a11=0.5; a21=0.5; a12=2; a22=0.8; 
A=matrix(c(a11,a12,a21,a22),2,2,byrow=TRUE);
A;  
cat(eigen(A)$values[1],"\n"); 

nx=30000; x=seq(-1500,1500,length=nx); h=x[2]-x[1];
sigmaS=4; sigmaL=10;  

### Represent dispersal as sparse matrices acting on population vector
### using bandSparse in the Matrix package. Otherwise dispersal 
### requires a nested loop and is very slow.  
### Dispersal range is limited here to 6*sigma, so the dispersal
### matrix is zero except near the diagonal. 

# Make dispersal matrix for small individuals 
ncol=max(5,round(6*sigmaS/h)); #limit dispersal range
bMat=matrix(0,nx,ncol);
for(j in 1:ncol) {
	dij=h*j; 
	pij=h*dnorm(dij,mean=0,sd=sigmaS) 
	bMat[,j]=pij;
}	
bMat=as.data.frame(bMat); 
pijS <- bandSparse(nx, k = 1:ncol, diag = bMat,symmetric=TRUE)
diag(pijS) <- h*dnorm(0,mean=0,sd=sigmaS); 

# Make dispersal matrix for large individuals 
ncol=max(5,round(6*sigmaL/h)); #limit dispersal range
bMat=matrix(0,nx,ncol);
for(j in 1:ncol) {
	dij=h*j; 
	pij=h*dnorm(dij,mean=0,sd=sigmaL) 
	bMat[,j]=pij;
}	
bMat=as.data.frame(bMat); 
pijL <- bandSparse(nx, k = 1:ncol, diag = bMat,symmetric=TRUE)
diag(pijL) <- h*dnorm(0,mean=0,sd=sigmaL); 
sum(pijS[,100]); sum(pijL[,100]); 

n.yrs <-100
pop.size.t <- rcrit.t <- matrix(NA,n.yrs)

XS=XL=XSnew=XLnew=as.numeric(abs(x)<2); 
## iterate the model 
yr <- 1;
while(yr <= n.yrs) {

    ## Calculate local population density
    X = (XS+XL); 
	if(yr%%25==0) plot(x,X,yaxs="i",main=yr);
	
    xfar <- x[which(X > 0.01)];
    rcrit.t[yr]=diff(range(xfar))/2; 

    ## density-dependent mortality 
    survivalProb <- 1/(1+0.5*pmax(as.numeric(X)-3,0)); 
 
    ## local stage transitions 
	XSnew=survivalProb*(a11*XS + a12*XL);
	XLnew=survivalProb*(a21*XS + a22*XL); 
	
	## everybody move 
	XS=pijS%*%XSnew; XL=pijL%*%XLnew; 

    ## Store population size
    pop.size <- sum(X) 
    pop.size.t[yr] <- pop.size
    
    if(yr%%25==0) cat(paste(yr, pop.size.t[yr], "\n", sep=" "))
    yr <- yr+1
   
}

## Plot the growth of occupied area over time 
plot(1:n.yrs,rcrit.t,type="l",lty=c(1,2,3,4),col="black",xlab="Years",
ylab="Radius of area occupied"); 

# Estimate spread rate from the last 50 years of spread. 
t.end <- seq(n.yrs-50,n.yrs,by=1);
fit<- lm(rcrit.t[t.end]~t.end)
spread.rate <- fit$coef[2]

## moment generating function for a Gaussian distribution 
## with mean 0, variance sigma^2. 
Gaussmgf <- function(s,sigma) {  exp(sigma^2*s^2/2) }

## wave speed c(s)
cs <- function(s) { 
	MS=Gaussmgf(s,sigmaS);
	ML=Gaussmgf(s,sigmaL); 
	Ma11=MS*a11; Ma12=MS*a12;
	Ma21=ML*a21; Ma22=ML*a22;
	MA=matrix(c(Ma11,Ma12,Ma21,Ma22),2,2,byrow=TRUE); 
	L1 = max(abs(eigen(MA)$values)); 
    return((1/s)*log(L1)) 
   }
cs = Vectorize(cs,"s"); 
plot(cs,0.05,2);

## Find the asymptotic wave speed c*(s) 
out=optimize(cs,lower=0.01,upper=2); 
cat(out$objective,"\n"); 
cat(spread.rate,"\n"); 

dev.new(height=8,width=6); 
par(mfrow=c(3,1),xaxs="i",cex.axis=1.65,cex.lab=1.65,bty="l",mar=c(4,4,2,4),mgp=c(2.5,1,0),adj=0.5); 
e=(x>0); 
matplot(x[e],cbind(XS[e],XL[e]),type="l",col="black",lty=c(2,1),lwd=c(2,2),
xlab="Spatial location x",ylab="Population density")
legend("topright",legend=c("Small  ","Large  "),lty=c(2,1),lwd=c(2,2),bty="n",cex=1.6); 
mtext("A)",side=3,adj=0,cex=1.35); 


e=(x>750)&(x<830); 
matplot(x[e],log10(cbind(XS[e],XL[e])),type="l",col="black",lty=c(2,1),lwd=c(2,2),
xlab="Spatial location x",ylab="log10(density)")
mtext("B)",side=3,adj=0,cex=1.35)

px=x[-1]; ds=diff(log(XS))/h; dl=diff(log(XL))/h; 
py=cbind(as.numeric(ds),as.numeric(dl)); 
e=px>0; 
matplot(px[e],py[e,1:2],type="l",lty=c(2,1),lwd=c(2,2),col="black",
xlab="Spatial location x",ylab="Derivative of log(density)");
abline(h=-out$minimum,lty=3,lwd=3);
mtext("C)",side=3,adj=0,cex=1.35)

dev.copy2eps(file="~/Repos/ipm_book/c8/figures/NC1DGaussianIDE.eps"); 

