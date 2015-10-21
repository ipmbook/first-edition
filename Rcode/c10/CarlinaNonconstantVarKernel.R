rm(list=ls(all=TRUE))
## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 

require(stats4); 
source("../Utilities/Standard Graphical Pars.R"); 

X=read.csv("CRLNAGRW.csv");
size1 <- X$lst; size2 <- X$lst1; #log-transformed sizes 

graphics.off; set_graph_pars("panel2"); 

####################################################### 
### Pilot fit: a linear model with nonconstant variance 
####################################################### 
LinearLogLik=function(a,b,sigma0,d) {
	loglik=sum(dnorm(size2,mean=a+b*size1,sd=sigma0*exp(-d*size1),log=TRUE))
	return(-loglik)
}

fit1 <- mle(LinearLogLik,start=list(a=1,b=1,sigma0=1,d=0.2), method="Nelder-Mead",
    control=list(trace=4,maxit=2500));	
parms <- as.numeric(coef(fit1));     
a=parms[1]; b=parms[2]; sigma0=parms[3]; d=parms[4];  

sdG <- function(z) sigma0*exp(-d*z)
scaledErr <- (size2-(a+b*size1))/sdG(size1); 
h <- bw.SJ(scaledErr); sdhat=sd(scaledErr); 
alpha <- sdhat/sqrt(sdhat^2+h^2); 
hResids <- alpha*scaledErr;
         
dfun <- function(z1,z) (1/sdG(z))*mean(dnorm(z1/sdG(z),mean=hResids,sd=h)) 
dfun <- Vectorize(dfun,"z1");
G_z1z <- function(z1,z) dfun(z1-(a+b*z),z);
G_z1z_Gaussian <- function(z1,z) dnorm(z1,mean=a+b*z,sd=sigma0*exp(-d*z));

z1vals <- seq(1,6,length=200); 
Ky <- cbind(G_z1z(z1vals,2),G_z1z(z1vals,3),G_z1z(z1vals,4));
Gy <- G_z1z_Gaussian(z1vals,2);
graphics.off(); dev.new(height=5,width=7); 
set_graph_pars("panel1"); par(cex.lab=1.35,pty="m"); 
matplot(z1vals,cbind(Ky,Gy),type="l",lty=c(1,1,1,2),col="black",lwd=c(1,2,3,1),
    xlab="Final size z'", ylab="Probability density")
   
legend("topleft",c("G(z',2)","G(z',3)","G(z',4)","Gaussian"), lwd=c(1,2,3,1),
lty=c(1,1,1,2), bty="n",cex=1.3); 
 
dev.copy2eps(file="../../c10/figures/CarlinaNonParKernel.eps")




 