rm(list=ls(all=TRUE))

## Working directory must be set here
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 

require(subplex); require(maxLik); 

# "True" population parameters
a=0.9; mu = 2; sigma=0.8; q=0.85; Delta=0.6; nu=.5; tau=.5; Kf=8; 
parms=c(a,mu,sigma,q,Delta,nu,tau,Kf); 
names(parms)<-c("a","mu","sigma","q","Delta","nu","tau","Kf"); 

# Size range and bins for Poisson approximate likelihood  
L=0; U=30; 
breaks=seq(L,U,length=61); h=breaks[2]-breaks[1]; 
mids=breaks[-1]-0.5*h;
rights=breaks[-1]; lefts=breaks[-length(breaks)];

######## IPM functions 
# growth 
f=function(y,x,parms) dnorm(y, parms["a"]*x+parms["mu"], parms["sigma"]) 

# fecundity
b = function(x,parms) (x^4)/(parms["Kf"]^4+x^4)*parms["Delta"];

# offspring size distribution 
g=function(y,parms) dlnorm(y,meanlog=log(parms["nu"]),sdlog=parms["tau"]); 

# kernel: survival*(growth + fecundity*offspring.size)
K=function(y,x,parms) parms["q"]*(f(y,x,parms) + b(x,parms)*g(y,parms)) 

######## Function to simulate the population 
xnew=function(x,parms) {
        # random parameter variation, pt=parameters in year t (iid) 
		pt=exp(log(parms)+0.05*rnorm(length(parms))); 
        nx=length(x);

		live=runif(nx)<pt["q"]; # survival 
		x.survive=x[live];      # size of survivors 
		
        # offspring 
		nkids=rpois(1,sum(b(x.survive,parms)));
		kidsize=rlnorm(nkids,meanlog=log(pt["nu"]),sdlog=pt["tau"]);
        
        # growth of survivors; kill any with new size < 0
		newx.survive=rnorm(length(x.survive),mean=pt["mu"]+pt["a"]*x.survive,sd=pt["sigma"]);
		newx.survive[newx.survive>U] <- U; 
		newx.survive=newx.survive[newx.survive>0]; 
		
	 	# return everyone (don't know newborns from survivors) 
		return(sort(c(kidsize,newx.survive))); 
}		

######## Approximate likelihood function (independent Poisson bin counts)  

# Likelihood function for one time step  '
# x0 is population at time t, x1 is population at (t+1)
# Looks to global environment for bins and midpoints 
LogLik1=function(x0,x1,parms) {
		nbins=length(mids); means=rep(0,nbins);
		names(parms)<-c("a","mu","sigma","q","Delta","nu","tau","Kf"); 
		
		# expected number in each bin at (t+1), based on IPM and x0 
		# integral over bin is done with 3-point trapezoid rule 
		for(k in seq_along(x0)){
		    means=means+0.25*(K(lefts,x0[k],parms)+K(rights,x0[k],parms))+0.5*K(mids,x0[k],parms)
		}  
		# Poisson approximate likelihood 
		obs=hist(x1,breaks=breaks,plot=FALSE)$counts;
		sum(dpois(obs,lambda=h*means,log=TRUE))
}		

# Likelihood function for 7 time steps 
logLikfun=function(p) {
	parms=exp(p); LL=0;
	for(j in 1:7) {LL=LL+LogLik1(Pop[[j]],Pop[[j+1]],parms)}
	if(verbose.objfun) cat(rep,parms,LL,"\n"); 
	return(LL); 
}

########################################################## 
# Fit repeatedly and estimate parameters  
##########################################################
Pop=list(8); # hold population point-pattern at 8 times 

nrep=24; 
estimates=ses=matrix(0,nrep,length(parms)); converged=rep("0",nrep); 

# simulate population data from true parameters 
for(rep in 1:nrep) {
    Pop[[1]]=rnorm(1000,mean=3,sd=0.8); 
    if(rep<=12) {for(j in 1:25) Pop[[1]]=xnew(Pop[[1]],parms)}
    Pop[[1]]=sample(Pop[[1]],250,replace=FALSE); 

    plot(Pop[[1]],rep(0,length(Pop[[1]])),xlim=c(-1,30),ylim=c(0,4),pch="+");  
    for(j in 2:8) {
		Pop[[j]]=xnew(Pop[[j-1]],parms); P1=Pop[[j]]; 
		points(P1,rep(j-1,length(P1)),pch="+");  
		text(max(P1)+2,j-1,labels=as.character(length(P1)),col="blue",cex=1)
     }
	 
	cat("*** Version 2, Population range ", range(Pop), " $$rep ",rep,"\n"); 
     			

# maximize the likelihood, using repeated Nelder-Mead to hopefully find the global
# minimum, and then BFGS to hone in on the minimum and compute the Hessian 

verbose.objfun=FALSE; 
out=maxLik(logLikfun,start=log(parms),method="NM",print.level=1,iterlim=2500,finalHessian=FALSE);
#out=maxLik(logLikfun,start=out$estimate,method="NM",print.level=1,iterlim=2500,finalHessian=FALSE);
out=maxLik(logLikfun,start=out$estimate,method="NM",print.level=1,iterlim=2500,finalHessian=FALSE);
verbose.objfun=TRUE;
out=maxLik(logLikfun,start=out$estimate,method="BFGS",print.level=1,iterlim=2500,finalHessian=TRUE);

est.parms=exp(out$estimate); est.parms; parms;
se=sqrt(diag(vcov(out)));  	
estimates[rep,]=est.parms;
ses[rep,]=se; 
z=(log(est.parms)-log(parms))/se;  
converged[rep] = (summary(out)$returnMessage=="successful convergence ")
cat("Replicate ", rep,z,"\n"); 
}

save.image(file="CountIPMfit.Rdata")




