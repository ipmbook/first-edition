## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Evolutionary calculations for the Carlina stochastic IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
setwd("~/Repos/ipm_book/Rcode/c2")
source("Monocarp Demog Funs.R");

source("../utilities/Standard Graphical Pars.R");

# Set simulation parameters
setwd("~/Repos/ipm_book/Rcode/c10")
setwd("~/temp/karen"); 
rm(list=ls(all=TRUE))
library(MASS)
library(glmmML)
library(mgcv)

dataf<-read.csv("Carlina survive or flower.fre",header=T)

attach(dataf)

#fit model ignoring year effects, used in ESS analysis

fit.flow<-glm(flow~lst,family=binomial)

fit.flow.r<-glmmML(flow~lst,data=dataf,family=binomial,cluster=yeart)

s.dataf<-dataf[flow==1,]

size.fl<-lst[flow==1]

s.dataf$flow<-0
names(s.dataf)[5]<-"die"

detach(dataf)

dataf<-read.csv("Carlina survive or died.txt",header=T)

attach(dataf)

dataf<-rbind(dataf,s.dataf)
attach(dataf)

minsize<-1.4481; 
maxsize<-5.170528;

n.big.matrix<-50;

#actual numbers of recruits per year

nrec<-c(20,42,12,17,8,19,58,45,44,2,56,25,75,92,94,6,4,34,104)

store.nrec<-nrec;

#store parameters in p.vec

p.vec.names<-rep(NA,11)
p.vec<-rep(NA,11);
sd.vec<-rep(NA,11)

p.vec[1]<- -2.284390			; p.vec.names[1]<-"1st survival param";		sd.vec[1]<-1.161576
p.vec[2]<- 0.9040202 			; p.vec.names[2]<-"2nd survival param";		sd.vec[2]<-0.4136052
p.vec[3]<- -16.190				; p.vec.names[3]<-"1st flow param    ";		sd.vec[3]<-1.029
p.vec[4]<- 3.883 				; p.vec.names[4]<-"2nd flow param    ";		sd.vec[4]<-0.0
p.vec[5]<- 1.135735	 			; p.vec.names[5]<-"ag                ";		sd.vec[5]<-0.1928372
p.vec[6]<- 0.734684				; p.vec.names[6]<-"bg                ";		sd.vec[6]<-0.1329790
p.vec[7]<- 0.08103536			; p.vec.names[7]<-"sigma2 growth     ";		sd.vec[7]<-0.0
p.vec[8]<- 1					; p.vec.names[8]<-"intercept seeds   ";		sd.vec[8]<-0.0       
p.vec[9]<- 2					; p.vec.names[9]<-"slope seeds       ";		sd.vec[9]<-0.0

p.vec[10]<-3.162432 			; p.vec.names[10]<-"mean kids size   ";		sd.vec[10]<-0.2706014
p.vec[11]<-0.248814 			; p.vec.names[11]<-"sigma2 kids size ";		sd.vec[11]<-0.0

#parameters for the probability of flowering, from lmer, fitting both slope and intercept as random

#p.vec[3]<- -15.7621			; p.vec.names[3]<-"1st flow param    ";		sd.vec[3]<-2.2361e-05
#p.vec[4]<- 3.7931 			; p.vec.names[4]<-"2nd flow param    ";		sd.vec[4]<-2.7432e-01



store.p.vec<-p.vec       #stored p.vec don't change!!
store.sd.vec<-sd.vec     #stored sd.vec don't change!!


##########################################################################
#Variance-covariance matrix for random effects estimated simultaneously using WinBugs

VarCovar.grow.rec<-matrix(c(0.03718619,0.04022101,0.04022101,0.07322512),nrow=2)

chol.VarCovar.grow.rec=chol(VarCovar.grow.rec)

store.VarCovar.grow.rec=VarCovar.grow.rec

mean.grow.rec<-c(1.135735,3.162432)

##########################################################################
##apply bias correction to estimated parameter sds

#sd.vec<-sqrt(16/15)*sd.vec

#VarCovar.grow.rec<-(16/15)*VarCovar.grow.rec

p.est.DI<-0.00095;       #density independent estimate used to roughly match log(lamdda s)=0.13 from field population

#p.est.DI<-0.00065        #density independent estimate with log(lamdda s)~0 used for testing


# Part (II) ##############################################################
# Compute the kernel component functions from the fitted models

sx<-function(x,params) {
	u<-exp(params[1]+params[2]*x);
return(u/(1+u));
}

fx<-function(x,params) {
	u<-exp(params[3]+params[4]*x);
              
return(u/(1+u));
}

gyx<-function(y,x,params) {
	mux<-params[5]+params[6]*x;
	sigmax2<-params[7];
	sigmax<-sqrt(sigmax2);
	fac1<-sqrt(2*pi)*sigmax;
	fac2<-((y-mux)^2)/(2*sigmax2);
	return(exp(-fac2)/fac1);
}

#assumes death before flowering

pyx<-function(y,x,params) { return(sx(x,params)*(1-fx(x,params))*gyx(y,x,params)) }

fecyx<-function(y,x,params) {
	nkids<-exp(params[8]+params[9]*x);
	kidsize.mean<- params[10];
	kidsize.var<- params[11]; 
	fac1<-sqrt(2*pi)*sqrt(kidsize.var);
	fac2<-((y-kidsize.mean)^2)/(2*kidsize.var);
	f<-nkids*exp(-fac2)/fac1;
	return(f);
}

fdyx<-function(y,x,params) {
	kidsize.mean<- params[10];
	kidsize.var<- params[11]; 
	fac1<-sqrt(2*pi)*sqrt(kidsize.var);
	fac2<-((y-kidsize.mean)^2)/(2*kidsize.var);
	f<-exp(-fac2)/fac1;
	return(f);
}

pfyx<-function(y,x,params,p.est) {
	nkids<-exp(params[8]+params[9]*x);
	kidsize.mean<- params[10];
	kidsize.var<- params[11]; 
	fac1<-sqrt(2*pi)*sqrt(kidsize.var);
	fac2<-((y-kidsize.mean)^2)/(2*kidsize.var);
	f<-(p.est*nkids*exp(-fac2)/fac1)-gyx(y,x,params)
	return(f);
}


fyx<-function(y,x,params) {
	nkids<-exp(params[8]+params[9]*x);
	kidsize.mean<- params[10];
	kidsize.var<- params[11]; 
	fac1<-sqrt(2*pi)*sqrt(kidsize.var);
	fac2<-((y-kidsize.mean)^2)/(2*kidsize.var);
	f<-sx(x,params)*fx(x,params)*nkids*exp(-fac2)/fac1;
	return(f);
}

sens.pf=function(y,x,params) {
	nkids<-exp(params[8]+params[9]*x);
	kidsize.mean<- params[10];
	kidsize.var<- params[11]; 
	fac1<-sqrt(2*pi)*sqrt(kidsize.var);
	fac2<-((y-kidsize.mean)^2)/(2*kidsize.var);
	f<-sx(x,params)*fx(x,params)*((nkids*exp(-fac2)/fac1)-gyx(y,x,params));
	return(f);
}

sens.fn=function(y,x,params) {
	kidsize.mean<- params[10];
	kidsize.var<- params[11]; 
	fac1<-sqrt(2*pi)*sqrt(kidsize.var);
	fac2<-((y-kidsize.mean)^2)/(2*kidsize.var);
	f<-sx(x,params)*fx(x,params)*exp(-fac2)/fac1;
	return(f);
}
	

############## The 'big matrix' M of size n x n 
bigmatrix<-function(n,params) {
# upper and lower integration limits
	L<-minsize; U<-maxsize;h=(U-L)/n
	
# boundary points b and mesh points y
	b<-L+c(0:n)*h;
	y<-0.5*(b[1:n]+b[2:(n+1)]);

# construct the matrix
	P<-h*outer(y,y,pyx,params=params)
	B<-h*outer(y,y,fyx,params=params)
	spf<-h*outer(y,y,sens.pf,params=params)
	M<-P+B
	
	return(list(matrix=M,meshpts=y,Pmatrix=P,Bmatrix=B,spf=spf)); 
}


#iterate model

iterate.model<-function(params,n.est,store.env=F,DD=T) {

	if(store.env==T){
		params.yr<-matrix(NA,ncol=11,nrow=n.est)
		p.est.yr<-rep(NA,n.est)
	} else {
		params.yr<-NULL
		p.est.yr<-NULL
	}
	
	nt<-rep(1/n.big.matrix,n.big.matrix)
	Rt<-rep(NA,n.est)
	size.dist<-rep(0,n.big.matrix)
	size.dist.fl<-rep(0,n.big.matrix)
	N.year<-rep(NA,n.est)

	mean.sx=rep(0,n.big.matrix)
	mean.pfl=rep(0,n.big.matrix)
	mean.fn=rep(0,n.big.matrix)
	mean.gyx=matrix(0,n.big.matrix,n.big.matrix)

	var.sx=rep(0,n.big.matrix)
	var.pfl=rep(0,n.big.matrix)
	var.fn=rep(0,n.big.matrix)




	mean.K=matrix(0,nrow=n.big.matrix,ncol=n.big.matrix)
	mean.B=matrix(0,nrow=n.big.matrix,ncol=n.big.matrix)
	mean.P=matrix(0,nrow=n.big.matrix,ncol=n.big.matrix)

	for (year.t in 1:n.est){
		if(year.t%%250==0) cat("iterate: ", year.t,"\n");

		#params.t<-rnorm(11,params,sd.vec)

		params.t=qnorm(runif(11,0.001,0.999))*sd.vec+params

		#sample from multivariate normal distribution for growth intercepts and yearly mean recruit size
		#params.t[c(5,10)]<-mvrnorm(1,mean.grow.rec,VarCovar.grow.rec)

		params.t[c(5,10)]=params[c(5,10)]+matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 
		
		year.K<-bigmatrix(n.big.matrix,params.t)

		if(DD==T) {

		bt<-sum(year.K$Bmatrix %*% nt)

		p.est<-nrec[sample(1:16,1)]/bt

		} else {

		p.est=p.est.DI

		}
		
		if (store.env==T){
			params.yr[year.t,]<-params.t
			p.est.yr[year.t]<-p.est
		}

		nt1<-(year.K$Pmatrix+p.est*year.K$Bmatrix) %*% nt

		#calculates the mean kernel

		mean.K=mean.K+year.K$Pmatrix+p.est*year.K$Bmatrix
		mean.B=mean.K+year.K$Bmatrix
		mean.P=mean.K+year.K$Pmatrix
		
		sum.nt1<-sum(nt1)
		sum.nt<-sum(nt)

		Rt[year.t]<-log(sum.nt1/sum.nt)

		N.year[year.t]<-sum.nt

		dist.fl.year<-nt*sx(year.K$meshpts,params.t)*fx(year.K$meshpts,params.t)

		size.dist<-size.dist+nt

		size.dist.fl<-size.dist.fl+dist.fl.year

		#calculates means of the various functions used to construct the kernel

		mean.sx=mean.sx+sx(year.K$meshpts,params.t)

		var.sx=var.sx+(sx(year.K$meshpts,params.t)-sx(year.K$meshpts,params))^2

		mean.pfl=mean.pfl+fx(year.K$meshpts,params.t)

		var.pfl=var.pfl+(fx(year.K$meshpts,params.t)-fx(year.K$meshpts,params))^2

		mean.fn=mean.fn+exp(params.t[8]+params.t[9]*year.K$meshpts);	

		var.fn=var.fn+(exp(params.t[8]+params.t[9]*year.K$meshpts)-exp(params[8]+params[9]*year.K$meshpts))^2
	
		mean.gyx=mean.gyx+outer(year.K$meshpts,year.K$meshpts,gyx,params=params.t)

		nt<-nt1

		if(DD==F) nt=nt/sum.nt1

		#normalised here to prevent population explosion
	
		#cat("iteration ",year.t,"  ","\n")

	}

	mean.K=mean.K/n.est
	mean.B=mean.B/n.est
	mean.P=mean.P/n.est

	mean.sx=mean.sx/n.est
	var.sx=var.sx/n.est

	mean.pfl=mean.pfl/n.est
	var.pfl=var.pfl/n.est

	mean.fn=mean.fn/n.est
	var.fn=var.fn/n.est

	mean.gyx=mean.gyx/n.est

	return(list(size.dist=size.dist,size.dist.fl=size.dist.fl,params.yr=params.yr,p.est.yr=p.est.yr,
		    N.year=N.year,Rt=Rt,meshpts=year.K$meshpts,mean.K=mean.K,mean.B=mean.B,mean.P=mean.P,mean.sx=mean.sx,
		    mean.pfl=mean.pfl,mean.gyx=mean.gyx,mean.fn=mean.fn,
		    var.sx=var.sx,var.pfl=var.pfl,var.fn=var.fn))
}

#####################################################################################################
##Perturbation analysis

stoc.pert.analysis<-function(n.est){

### Get series of year types and p.est values that define the environment ###
	iter<-iterate.model(p.vec,n.est,store.env=T,DD=T)

### Get wt and Rt time series ###
	wt<-matrix(1/n.big.matrix, nrow=n.est+1, ncol=n.big.matrix);
	Rt<-rep(NA, n.est);
	for (i in 1:n.est) {
		
		K<-bigmatrix(n.big.matrix,iter$params.yr[i,]);
		wt[i+1,]<-(K$Pmatrix+iter$p.est.yr[i]*K$Bmatrix)%*%wt[i,]
		Rt[i]<-sum(wt[i+1,]);
		wt[i+1,]<-wt[i+1,]/Rt[i];
		if(i%%250==0) cat("wt and Rt ",i,"\n")

	}

	
### Get vt time series ###
	vt<-matrix(1/n.big.matrix, nrow=n.est+1, ncol=n.big.matrix);
	for (i in (n.est+1):2) {
		
		K<-bigmatrix(n.big.matrix,iter$params.yr[i-1,]);
		vt[i-1,]<-vt[i,]%*%(K$Pmatrix+iter$p.est.yr[i-1]*K$Bmatrix);
		vt[i-1,]<-vt[i-1,]/sum(vt[i-1,]);
		if(i%%250==0) cat("vt  ",i,"\n")

	}


### Generate sensitivity and elasticity matrices ###
	sens.s<-matrix(0,nrow=n.big.matrix,ncol=n.big.matrix);
	sens.s.v2.0=sens.s
	elas.s<-sens.s;
	elas.s.mean<-sens.s;
	elas.s.gyx<-sens.s;
	elas.s.gyx.mean<-sens.s
	elas.s.sx=rep(0,n.big.matrix)
	elas.s.sx.mean=rep(0,n.big.matrix)
	elas.s.pf=rep(0,n.big.matrix)
	elas.s.pf.mean=rep(0,n.big.matrix)
	elas.s.fn=rep(0,n.big.matrix)
	elas.s.fn.mean=rep(0,n.big.matrix)
	sens.s.fn.v2.0=rep(0,n.big.matrix)
        elas.s.sd=sens.s; 
	kernel.sd=sens.s; 


	y=iter$meshpts
	h=iter$meshpts[2]-iter$meshpts[1]

	for (i in 1:n.est) {

		#standard calculations needed for the various formulae
	
		K<-bigmatrix(n.big.matrix,iter$params.yr[i,]);
		vt1.wt=vt[i+1,]%*%t(wt[i,])
		Rt.vt1.wt1=(Rt[i]*t(vt[i+1,])%*%wt[i+1,])[1]

		Kt=K$Pmatrix+iter$p.est.yr[i]*K$Bmatrix

		#calculation of the standard sensitivities and elasticities
	
		sens.s<-sens.s+vt1.wt/Rt.vt1.wt1;

		elas.s<-elas.s+Kt*(vt1.wt/Rt.vt1.wt1);

		elas.s.mean<-elas.s.mean+iter$mean.K*(vt1.wt/Rt.vt1.wt1);

	        elas.s.sd = elas.s.sd+(Kt-iter$mean.K)*(vt1.wt/Rt.vt1.wt1); 

		kernel.sd=kernel.sd+(Kt-iter$mean.K)^2	
	
		#calculates equation 9 from "More on elasticities" actually equation 10 but that's not numbered.

		sens.s.v2.0<-sens.s.v2.0-(vt1.wt * vt1.wt)/(Rt.vt1.wt1*Rt.vt1.wt1)/2;

		#survival elasticities, this implements equation 5 in "Even more on Elasticities" 
		#the complicated transpose thing allows me to recreate various parts of the kernel

		p.surv=sx(y,iter$params.yr[i,])
		
		C.1=t(p.surv*(1-fx(y,iter$params.yr[i,]))*t(outer(y,y,function(y,x,params) {gyx(y,x,params)*h},params=iter$params.yr[i,])))
		
		C.2=t(p.surv*fx(y,iter$params.yr[i,])*t(outer(y,y,function(y,x,params) {fecyx(y,x,params)*h},params=iter$params.yr[i,])))

		Ct.sx=C.1+iter$p.est.yr[i]*C.2

		elas.s.sx=elas.s.sx+wt[i,]*(vt[i+1,] %*% Ct.sx)/Rt.vt1.wt1

		#this can be calculated more efficiently using
		#elas.sx=elas.sx+wt[i,]*(vt[i+1,] %*% Kt)/Rt.vt1.wt1
		#but I'll use the more explicit version

		#now set the probability of survival to the mean value, and repeat the calculation to get elasticity to mean

		p.surv=iter$mean.sx
				
		C.1=t(p.surv*(1-fx(y,iter$params.yr[i,]))*t(outer(y,y,function(y,x,params) {gyx(y,x,params)*h},params=iter$params.yr[i,])))
		
		C.2=t(p.surv*fx(y,iter$params.yr[i,])*t(outer(y,y,function(y,x,params) {fecyx(y,x,params)*h},params=iter$params.yr[i,])))

		Ct.mean.sx=C.1+iter$p.est.yr[i]*C.2

		elas.s.sx.mean=elas.s.sx.mean+wt[i,]*(vt[i+1,] %*% Ct.mean.sx)/Rt.vt1.wt1

		#probability of flowering elasticities, this implements equation 7 in "Even more on Elasticities" 

		p.flow=fx(y,iter$params.yr[i,])

		C.pfl=t(sx(y,iter$params.yr[i,])*p.flow*t(outer(y,y,pfyx,params=iter$params.yr[i,],p.est=iter$p.est.yr[i])))*h

		elas.s.pf=elas.s.pf+wt[i,]*(vt[i+1,] %*% C.pfl)/Rt.vt1.wt1

		#this can be done more efficiently with
		#elas.pf=elas.pf+wt[i,]*(vt[i+1,] %*% K$spf)/Rt.vt1.wt1
		#but I'll use the more explicit version

		#now set the probability of flowering to the mean value, and repeat the calculation to get elasticity to mean

		p.flow=iter$mean.pfl

		C.mean.pfl=t(sx(y,iter$params.yr[i,])*p.flow*t(outer(y,y,pfyx,params=iter$params.yr[i,],p.est=iter$p.est.yr[i])))*h

		elas.s.pf.mean=elas.s.pf.mean+wt[i,]*(vt[i+1,] %*% C.mean.pfl)/Rt.vt1.wt1

		#growth elasticities, implements equation 9. Taking g(y,x) back inside the expectation we calculate the 
		#expectation by element-wise multiplication of vt1.wt by the K$Pmatrix, an alternative is to use the approach used below
		#for the elasticity to the mean but with
		#p.grow=outer(y,y,function(y,x,params) {gyx(y,x,params)},params=iter$params.yr[i,])

		elas.s.gyx=elas.s.gyx+(vt1.wt*K$Pmatrix)/Rt.vt1.wt1

		#set p.grow equal to the mean and repeat calculation

		p.grow=iter$mean.gyx

		C.mean.gyx=t(sx(y,iter$params.yr[i,])*(1-fx(y,iter$params.yr[i,]))*t(p.grow))*h

		elas.s.gyx.mean=elas.s.gyx.mean+vt1.wt*C.mean.gyx/Rt.vt1.wt1

		#seed production elasticities, perturbs the seed production function, fn -> fn+e*fn
	
		elas.s.fn=elas.s.fn+wt[i,]*(vt[i+1,] %*% (iter$p.est.yr[i]*K$Bmatrix))/Rt.vt1.wt1

		#put seed production equal to mean
		
		seeds=iter$mean.fn

		#seeds=exp(iter$params.yr[i,8]+iter$params.yr[i,9]*y)

		C.mean.fn=t(sx(y,iter$params.yr[i,])*fx(y,iter$params.yr[i,])*seeds*t(outer(y,y,fdyx,params=iter$params.yr[i,])))*h

		elas.s.fn.mean=elas.s.fn.mean+wt[i,]*(vt[i+1,] %*% (iter$p.est.yr[i]*C.mean.fn))/Rt.vt1.wt1

		#fn is constant so we use equation 13

		Ht=t(outer(y,y,sens.fn,params=iter$params.yr[i,]))*h

		tmp=(wt[i,]*(vt[i+1,] %*% Ht))^2/(Rt.vt1.wt1*Rt.vt1.wt1)

		sens.s.fn.v2.0=sens.s.fn.v2.0+tmp
		
		
		if(i%%250==0) cat("Elasticities and sensitivities ",i,"\n")

	}

	lam.stoch=exp(mean(iter$Rt[n.runin:n.est],na.rm=T))

	sens.s<-lam.stoch*sens.s/n.est;
	sens.s.v2.0<-lam.stoch*sens.s.v2.0/n.est;
	elas.s.sx<-elas.s.sx/n.est;
	elas.s.sx.mean<-elas.s.sx.mean/n.est;
	elas.s.pf<-elas.s.pf/n.est;
	elas.s.pf.mean<-elas.s.pf.mean/n.est;
	elas.s.gyx<-elas.s.gyx/n.est;
	elas.s.gyx.mean<-elas.s.gyx.mean/n.est
	elas.s.fn<-elas.s.fn/n.est
	elas.s.fn.mean<-elas.s.fn.mean/n.est
	sens.s.fn.v2.0=-(lam.stoch/2)*sens.s.fn.v2.0/n.est
	elas.s=elas.s/n.est;
	elas.s.mean=elas.s.mean/n.est;
	elas.s.sd = elas.s.sd/n.est; 
	
	kernel.sd=sqrt(kernel.sd/(n.est-1)); 

#### SPE: added some new returns: elasticity to sd, kernel sd. 
return(list(sens.s=sens.s,elas.s=elas.s,elas.s.mean=elas.s.mean,meshpts=iter$meshpts,mean.K=iter$mean.K,
	    sens.s.v2.0=sens.s.v2.0,elas.s.sx=elas.s.sx,elas.s.pf=elas.s.pf,elas.s.fn=elas.s.fn,elas.s.gyx=elas.s.gyx,
	    elas.s.sx.mean=elas.s.sx.mean,elas.s.pf.mean=elas.s.pf.mean,elas.s.gyx.mean=elas.s.gyx.mean,	
	    elas.s.fn.mean=elas.s.fn.mean,elas.s.sd=elas.s.sd,kernel.sd=kernel.sd,sens.s.fn.v2.0=sens.s.fn.v2.0,
	    lam.s=lam.stoch,mean.gyx=iter$mean.gyx,mean.sx=iter$mean.sx,mean.pfl=iter$mean.pfl,mean.fn=iter$mean.fn,
	    var.sx=iter$var.sx,var.pfl=iter$var.pfl,var.fn=iter$var.fn))
}


#####################################################################################################
##Perturbation analysis

stoc.pert.analysis.lower.p<-function(n.est,no.p,delta){

### Get series of year types and p.est values that define the environment ###
	iter<-iterate.model(p.vec,n.est,store.env=T,DD=T)

	lam.stoch=exp(mean(iter$Rt[n.runin:n.est],na.rm=T))



### Get wt and Rt time series ###
	wt<-matrix(1/n.big.matrix, nrow=n.est+1, ncol=n.big.matrix);
	Rt<-rep(NA, n.est);
	for (i in 1:n.est) {

		K<-bigmatrix(n.big.matrix,iter$params.yr[i,]);
		wt[i+1,]<-(K$Pmatrix+iter$p.est.yr[i]*K$Bmatrix)%*%wt[i,]
		Rt[i]<-sum(wt[i+1,]);
		wt[i+1,]<-wt[i+1,]/Rt[i];
		#cat("wt and Rt ",i,"\n")

	}

	
### Get vt time series ###
	vt<-matrix(1/n.big.matrix, nrow=n.est+1, ncol=n.big.matrix);
	for (i in (n.est+1):2) {
		
		K<-bigmatrix(n.big.matrix,iter$params.yr[i-1,]);
		vt[i-1,]<-vt[i,]%*%(K$Pmatrix+iter$p.est.yr[i-1]*K$Bmatrix);
		vt[i-1,]<-vt[i-1,]/sum(vt[i-1,]);
		#cat("vt  ",i,"\n")

	}

elas.s.p=rep(NA,no.p)
elas.s.p.mean=rep(NA,no.p)
sens.s.p=rep(NA,no.p)
sens.s.p.v2.0=rep(NA,no.p)
sens.s.p.v2.0.old=rep(NA,no.p)
t1.s.p=rep(NA,no.p)
t2.s.p=rep(NA,no.p)
var.vt.partial.k.wt=rep(NA,no.p)


for(ps in 1:no.p){

	p.pert=rep(0,11)
	p.pert[ps]=delta



### Generate sensitivity and elasticity matrices ###
	
	elas.p=0;
	elas.p.mean=0;
	sens.p=0
	sens.p.v2.0=0
	sens.p.v2.0.non.lin=0
	
	store.partial.k=rep(0,n.est)
	

	for (i in 1:n.est) {
	
		K<-bigmatrix(n.big.matrix,iter$params.yr[i,]);
		
		K.ups=bigmatrix(n.big.matrix,(1+p.pert)*iter$params.yr[i,]);
		K.downs=bigmatrix(n.big.matrix,(1-p.pert)*iter$params.yr[i,]);

		K.0=K$Pmatrix+iter$p.est.yr[i]*K$Bmatrix
		K.up=K.ups$Pmatrix+iter$p.est.yr[i]*K.ups$Bmatrix
		K.down=K.downs$Pmatrix+iter$p.est.yr[i]*K.downs$Bmatrix

		partial.K=(K.up-K.down)/(2*delta*iter$params.yr[i,ps])
		
		partial.K2=(((K.up-K.0)/(delta*iter$params.yr[i,ps]))-
				   ((K.0-K.down)/(delta*iter$params.yr[i,ps])))/
				   (delta*iter$params.yr[i,ps])

		Rt.vt1.wt=(Rt[i]*t(vt[i+1,])%*%wt[i+1,])[1]

		vt.partial.k.wt=vt[i+1,] %*% (partial.K %*% wt[i,])
		
		store.partial.k[i]=vt.partial.k.wt/Rt.vt1.wt

		vt.partial.k2.wt=vt[i+1,] %*% (partial.K2 %*% wt[i,])
		
		sens.p=sens.p+(vt.partial.k.wt)/(Rt.vt1.wt)
		
		sens.p.v2.0=sens.p.v2.0+(vt.partial.k.wt)*(vt.partial.k.wt)/(Rt.vt1.wt*Rt.vt1.wt)
		
		sens.p.v2.0.non.lin=sens.p.v2.0.non.lin+vt.partial.k2.wt/Rt.vt1.wt

		elas.p<-elas.p+iter$params.yr[i,ps]*(vt.partial.k.wt)/(Rt.vt1.wt)

		elas.p.mean<-elas.p.mean+p.vec[ps]*(vt.partial.k.wt)/(Rt.vt1.wt)

		
		#cat("Elasticities  ",i,"\n")

	}

	t1.s.p[ps]=sens.p.v2.0/n.est
	t2.s.p[ps]=sens.p.v2.0.non.lin/n.est
	elas.s.p[ps]=elas.p/n.est;
	sens.s.p[ps]=lam.stoch*sens.p/n.est;
	elas.s.p.mean[ps]=elas.p.mean/n.est;
	sens.s.p.v2.0[ps]=(lam.stoch/2)*(sens.p.v2.0.non.lin-sens.p.v2.0)/n.est
	sens.s.p.v2.0.old[ps]=(lam.stoch/2)*(-sens.p.v2.0)/n.est
	var.vt.partial.k.wt[ps]=var(store.partial.k,na.rm=T)
	cat(sens.p.v2.0.non.lin/n.est)
	cat("para ",ps,"\n")
	}

return(list(elas.s.p=elas.s.p,elas.s.p.mean=elas.s.p.mean,sens.s.p=sens.s.p,sens.s.p.v2.0=sens.s.p.v2.0,
            sens.s.p.v2.0.old=sens.s.p.v2.0.old,lam.stoch=lam.stoch,t1.s.p=t1.s.p,t2.s.p=t2.s.p,mean.k=iter$mean.K,
            mean.B=iter$mean.B,mean.p=iter$mean.P,mean.p.est=mean(iter$p.est.yr),
            var.vt.partial.k.wt=var.vt.partial.k.wt))
}


det.pert.analysis<-function(A){

	w<-Re(eigen(A)$vectors[,1]); 
	v<-Re(eigen(t(A))$vectors[,1]);
	vw<-sum(v*w);
	s<-outer(v,w)
	sens<-s/vw
	lam<-Re(eigen(A)$values[1])
	elas<-(sens*A)/lam;
return(list(sens=sens,elas=elas,w=w,v=v,lam=lam))
}


#################################################################################################################################
#LTRE


LTRE.f<-function(n.est,p.no,delta){

### Get series of year types and p.est values that define the environment ###
	iter<-iterate.model(p.vec,n.est,store.env=T,DD=F)

p.pert=rep(0,11)
p.pert[p.no]=delta

partial.K=0

lams=rep(NA,n.est)

	for (i in 1:n.est) {
	
		K<-bigmatrix(n.big.matrix,iter$params.yr[i,]);

		Ki=K$Pmatrix+iter$p.est.yr[i]*K$Bmatrix

		lams[i]=Re(eigen(Ki)$values[1])
		
		K.ups=bigmatrix(n.big.matrix,(1+p.pert)*iter$params.yr[i,]);
		K.downs=bigmatrix(n.big.matrix,(1-p.pert)*iter$params.yr[i,]);

		K.up=K.ups$Pmatrix+iter$p.est.yr[i]*K.ups$Bmatrix
		K.down=K.downs$Pmatrix+iter$p.est.yr[i]*K.downs$Bmatrix


		partial.K=partial.K+(K.up-K.down)/(2*delta*iter$params.yr[i,p.no])

	}

	
	
	partial.K<-partial.K/n.est;

return(list(partial.K=partial.K,mean.K=iter$mean.K,var.lam=var(lams)))
}




#################################################################################################################################
##run model and plot some graphs

setwd("~/temp/ellner/stoch 2"); 

n.big.matrix=100
n.est<-5000
n.runin<-500

iter<-iterate.model(p.vec,n.est,store.env=F,DD=T)

#calculate stochastic growth rate

mean(iter$Rt[200:n.est],na.rm=T)

plot(iter$N.year,type="b",pch=19)

mean(iter$N.year[n.runin:n.est])

p.size.dist.fm<-iter$size.dist/sum(iter$size.dist)

mean.size<-sum(p.size.dist.fm*exp(iter$meshpts))
mean.size

p.size.dist.fl.fm<-iter$size.dist.fl/sum(iter$size.dist.fl)

mean.size.fl<-sum(p.size.dist.fl.fm*exp(iter$meshpts))
mean.size.fl

diff<-iter$meshpts[2]-iter$meshpts[1]

quartz()

par(mfrow=c(1,2),bty="l",pty="s")

hist(lst,breaks=10,main="",col="grey",ylim=c(0,420),xlab="Size")
lines(iter$meshpts,p.size.dist.fm*599/sum(p.size.dist.fm*diff))

#text(locator(1),"a)")

hist(size.fl,breaks=10,main="",col="grey",xlab="Flowering size")
lines(iter$meshpts,p.size.dist.fl.fm*21.4/sum(p.size.dist.fl.fm*diff))

#text(locator(1),"b)")

##################################################################################################################################
##Elasticity analysis 

#p.vec[3]=-14.15

stoc.pert=stoc.pert.analysis(n.est)
sum(stoc.pert$elas.s)
sum(stoc.pert$elas.s.sx)
sum(stoc.pert$elas.s.sx.mean)
sum(stoc.pert$elas.s.pf)
sum(stoc.pert$elas.s.pf.mean)
sum(stoc.pert$elas.s.gyx)
sum(stoc.pert$elas.s.gyx.mean)
sum(stoc.pert$elas.s.fn)
sum(stoc.pert$elas.s.fn.mean)


#the sx elasticities should sum to 1, and equal the column sums of the elasticity matrix

max(abs(apply(stoc.pert$elas.s,2,sum)-stoc.pert$elas.s.sx))

#seems OK

meshpts<-stoc.pert$meshpts; h=meshpts[2]-meshpts[1]; 
hinv=1/h; hi2=hinv*hinv; 

elas.s.sd=stoc.pert$elas.s.sd

elas.s.sd.diff=stoc.pert$elas.s-stoc.pert$elas.s.mean

plot(elas.s.sd,elas.s.sd.diff)

abline(0,1)

elas.s.sx.sd=stoc.pert$elas.s.sx-stoc.pert$elas.s.sx.mean

elas.s.pf.sd=stoc.pert$elas.s.pf-stoc.pert$elas.s.pf.mean

elas.s.fn.sd=stoc.pert$elas.s.fn-stoc.pert$elas.s.fn.mean

elas.s.gyx.sd=stoc.pert$elas.s.gyx-stoc.pert$elas.s.gyx.mean

sens.s.sd=stoc.pert$lam.s*stoc.pert$elas.s.sd/stoc.pert$kernel.sd

sens.s.v2=0.5*sens.s.sd/stoc.pert$kernel.sd

plot(sens.s.v2,stoc.pert$sens.s.v2.0)

sen.s.sx.mean=stoc.pert$elas.s.sx.mean*stoc.pert$lam.s/stoc.pert$mean.sx

sen.s.sx.sd=elas.s.sx.sd*stoc.pert$lam.s/sqrt(stoc.pert$var.sx)

sen.s.pf.mean=stoc.pert$elas.s.pf.mean*stoc.pert$lam.s/stoc.pert$mean.pf

sen.s.pf.sd=elas.s.pf.sd*stoc.pert$lam.s/sqrt(stoc.pert$var.pf)

sen.s.fn.mean=stoc.pert$elas.s.fn.mean*stoc.pert$lam.s/stoc.pert$mean.fn

sen.s.sx.var=sen.s.sx.sd/(2*sqrt(stoc.pert$var.sx))

sen.s.pf.var=sen.s.pf.sd/(2*sqrt(stoc.pert$var.pf))

prop.sx=stoc.pert$elas.s.sx.mean/(stoc.pert$elas.s.sx.mean+abs(elas.s.sx.sd))

prop.pf=stoc.pert$elas.s.pf.mean/(stoc.pert$elas.s.pf.mean+abs(elas.s.pf.sd))

prop.fn=stoc.pert$elas.s.fn.mean/(stoc.pert$elas.s.fn.mean+abs(elas.s.fn.sd))

prop.gyx=stoc.pert$elas.s.gyx.mean/(stoc.pert$elas.s.gyx.mean+abs(elas.s.gyx.sd))

arith.size=exp(stoc.pert$meshpts)

quartz()
par(mfrow=c(3,3),bty="l",pty="s",pch=19)

xlab.text="Size (mm)"

plot(arith.size,stoc.pert$elas.s.sx,type="l",xlab=xlab.text,ylab=expression(italic(e[S*","*s(x)])),log="x")

text(locator(1),"a)")

plot(arith.size,stoc.pert$elas.s.sx.mean,type="l",xlab=xlab.text,ylab=expression(italic(e[S*","*s(x)]^mu)),log="x")

text(locator(1),"b)")

plot(arith.size,elas.s.sx.sd,type="l",xlab=xlab.text,ylab=expression(italic(e[S*","*s(x)]^sigma)),log="x")

text(locator(1),"c)")

#plot(stoc.pert$meshpts,prop.sx,type="l",main="Prop mean",ylim=c(0,1),xlab=xlab.text,ylab=expression(italic(e[S*","*s(x)])))

plot(arith.size,stoc.pert$elas.s.pf,type="l",xlab=xlab.text,ylab=expression(italic(e[S*","*p[f]*(x)])),log="x")

text(locator(1),"d)")

plot(arith.size,stoc.pert$elas.s.pf.mean,type="l",xlab=xlab.text,ylab=expression(italic(e[S*","*p[f]*(x)]^mu)),log="x")

text(locator(1),"e)")

plot(arith.size,elas.s.pf.sd,type="l",xlab=xlab.text,ylab=expression(italic(e[S*","*p[f]*(x)]^sigma)),log="x")

text(locator(1),"f)")

#plot(stoc.elas$meshpts,prop.pf,type="l",main="Prop mean",ylim=c(0,1),xlab=xlab.text)

plot(arith.size,stoc.pert$elas.s.fn,type="l",xlab=xlab.text,ylab=expression(italic(e[S*","*f[n]*(x)])),log="x")

text(locator(1),"g)")

plot(arith.size,stoc.pert$elas.s.fn.mean,type="l",xlab=xlab.text,ylab=expression(italic(e[S*","*f[n]*(x)]^mu)),log="x")

text(locator(1),"h)")

#plot(stoc.pert$meshpts,elas.fn.sd,type="l",xlab=xlab.text)

plot(arith.size,stoc.pert$sens.s.fn.v2.0,type="l",xlab=xlab.text,ylab=expression(italic(s[S*","*f[n]*(x)]^{sigma^2*","*0})),log="x")

text(locator(1),"i)")

#plot(stoc.pert$meshpts,prop.fn,type="l",main="Prop mean",ylim=c(0,1),xlab=xlab.text)

#savePlot("ElasticityImagesFunctions.pdf",type="pdf"); 
#savePlot("ElasticityImagesFunctions.bmp",type="bmp"); 
#savePlot("ElasticityImagesFunctions.png",type="png"); 

quartz()
par(mfrow=c(1,3),bty="l",pty="s",pch=19)

cor.test(stoc.pert$elas.s.sx.mean,elas.s.sx.sd)
cor.test(stoc.pert$elas.s.pf.mean,elas.s.pf.sd)
cor.test(stoc.pert$elas.s.fn.mean,stoc.pert$sens.s.fn.v2.0)

plot(stoc.pert$elas.s.sx.mean,elas.s.sx.sd)
abline(0,-1)
abline(1,-1)
plot(stoc.pert$elas.s.pf.mean,elas.s.pf.sd)
abline(0,-1)
#abline(1,-1)
plot(stoc.pert$elas.s.fn.mean,stoc.pert$sens.s.fn.v2.0)

quartz()
par(mfrow=c(1,3),bty="l",pty="s",pch=19)

cor.test(sen.s.sx.mean,sen.s.sx.var)
cor.test(sen.s.pf.mean,sen.s.pf.var)
cor.test(sen.s.fn.mean,stoc.pert$sens.s.fn.v2.0)


plot(sen.s.sx.mean,sen.s.sx.var)
plot(sen.s.pf.mean,sen.s.pf.var)
plot(sen.s.fn.mean,stoc.pert$sens.s.fn.v2.0)

quartz()
par(mfrow=c(3,2),bty="l",pty="s",pch=19,cex.lab=1.2)

xlab.text="Size (mm) "

plot(arith.size,sen.s.sx.mean,type="l",xlab=xlab.text,ylab=expression(italic(s[S*","*s(x)]^mu)),log="x")

text(locator(1),"a)")

plot(arith.size,sen.s.sx.var,type="l",xlab=xlab.text,ylab=expression(italic(s[S*","*s(x)]^sigma^2)),log="x")

text(locator(1),"b)")

#plot(stoc.pert$meshpts,prop.sx,type="l",main="Prop mean",ylim=c(0,1),xlab=xlab.text,ylab=expression(italic(e[S*","*s(x)])))

plot(arith.size,sen.s.pf.mean,type="l",xlab=xlab.text,ylab=expression(italic(s[S*","*p[f]*(x)]^mu)),log="x")

#abline(v=-p.vec[3]/p.vec[4])

#abline(h=0)

text(locator(1),"c)")

plot(arith.size,sen.s.pf.var,type="l",xlab=xlab.text,ylab=expression(italic(s[S*","*p[f]*(x)]^sigma^2)),log="x")

text(locator(1),"d)")

#plot(stoc.elas$meshpts,prop.pf,type="l",main="Prop mean",ylim=c(0,1),xlab=xlab.text)

plot(arith.size,sen.s.fn.mean,type="l",xlab=xlab.text,ylab=expression(italic(s[S*","*f[n]*(x)]^mu)),log="x")

text(locator(1),"e)")

#plot(stoc.pert$meshpts,elas.fn.sd,type="l",xlab=xlab.text)

plot(arith.size,stoc.pert$sens.s.fn.v2.0,type="l",xlab=xlab.text,ylab=expression(italic(s[S*","*f[n]*(x)]^{sigma^2*","*0})),log="x")

text(locator(1),"f)")

#savePlot("SensitivityImagesFunctions.pdf",type="pdf"); 
#savePlot("SensitivityImagesFunctions.bmp",type="bmp"); 
#savePlot("SensitivityImagesFunctions.png",type="png"); 

matrix.image.plot=function(A,x=NULL,y=NULL,xlab=NULL,ylab=NULL,
colors=rainbow(60,start=0.66667,end=0,s=0.9,v=.8),
add.contour=FALSE,...) {
	if(is.null(x)) x=1:ncol(A);
	if(is.null(y)) y=1:nrow(A); 
	if(is.null(xlab)) xlab="column"; if(is.null(ylab)) ylab="row"; 
	nx=length(x); ny=length(y); 
	xlimits=c(min(x)-0.5*(x[2]-x[1]), max(x)+0.5*(x[nx]-x[nx-1]));
 	ylimits=c(min(y)-0.5*(y[2]-y[1]), max(y)+0.5*(y[ny]-y[ny-1]));
	image(x,y,t(A),col=colors,xlim=xlimits,ylim=rev(ylimits),
	xlab=xlab,ylab=ylab) 
	if(add.contour) {
		contour(x,y,t(A),add=TRUE,labcex=0.8,...);
	}
}

xlab.text="Size year t (mm)"

ylab.text="Size year t+1 (mm)"


quartz()
par(mfrow=c(2,2),pty="s",pch=19,cex.lab=1.2,cex.axis=1.2,xaxt="n",yaxt="n")

matrix.image.plot(x=stoc.pert$meshpts,y=stoc.pert$meshpts,hi2*(stoc.pert$elas.s.gyx),xlab=xlab.text,ylab=ylab.text,col=gray(100:400/400),add.contour=TRUE,nlevels=5)
#contour(stoc.pert$meshpts,stoc.pert$meshpts,hi2*t(stoc.pert$elas.s.gyx),add=T,nlevels=4,labcex = 0.8,vfont = c("sans serif", "plain"))
par(xaxt="s",yaxt="s");
axis(side=1,at=log(c(5,10,50)),labels=c("5","10","50"))
axis(side=2,at=log(c(5,10,50)),labels=c("5","10","50"))


#text(2,2,"a)",col="white",cex=1.5)

par(pty="s",pch=19,cex.lab=1.2,cex.axis=1.2,xaxt="n",yaxt="n")

matrix.image.plot(x=stoc.pert$meshpts,y=stoc.pert$meshpts,hi2*(stoc.pert$elas.s.gyx.mean),xlab=xlab.text,ylab=ylab.text,col=gray(100:400/400),add.contour=TRUE,nlevels=5)
#contour(stoc.pert$meshpts,stoc.pert$meshpts,hi2*t(stoc.pert$elas.s.gyx.mean),add=T,nlevels=4,labcex = 0.8,vfont = c("sans serif", "plain"))
par(xaxt="s",yaxt="s");
axis(side=1,at=log(c(5,10,50)),labels=c("5","10","50"))
axis(side=2,at=log(c(5,10,50)),labels=c("5","10","50"))

#text(2,2,"b)",col="white",cex=1.5)

par(pty="s",pch=19,cex.lab=1.2,cex.axis=1.2,xaxt="n",yaxt="n")

#elas.s.gyx.sd[elas.s.gyx.sd>-0.001]=0

matrix.image.plot(x=stoc.pert$meshpts,y=stoc.pert$meshpts,hi2*(elas.s.gyx.sd),xlab=xlab.text,ylab=ylab.text,col=gray(100:400/400),add.contour=TRUE,levels=c(-0.5,-1,-1.5,-2,-2.5,-3))
#contour(stoc.pert$meshpts,stoc.pert$meshpts,hi2*t(elas.s.gyx.sd),add=T,nlevels=4,labcex = 0.8,vfont = c("sans serif", "plain"))
par(xaxt="s",yaxt="s");
axis(side=1,at=log(c(5,10,50)),labels=c("5","10","50"))
axis(side=2,at=log(c(5,10,50)),labels=c("5","10","50"))

#text(2,2,"c)",cex=1.5)

par(pty="s",pch=19,cex.lab=1.2,cex.axis=1.2,xaxt="n",yaxt="n")

sens.s.gyx=stoc.pert$elas.s.gyx.mean/stoc.pert$mean.gyx

matrix.image.plot(x=stoc.pert$meshpts,y=stoc.pert$meshpts,hi2*(sens.s.gyx),xlab=xlab.text,ylab=ylab.text,col=gray(100:400/400),add.contour=TRUE,nlevels=5)
#contour(stoc.pert$meshpts,stoc.pert$meshpts,hi2*t(sens.s.gyx),add=T,nlevels=4,labcex = 0.8,vfont = c("sans serif", "plain"))
par(xaxt="s",yaxt="s");
axis(side=1,at=log(c(5,10,50)),labels=c("5","10","50"))
axis(side=2,at=log(c(5,10,50)),labels=c("5","10","50"))

#text(2,2,"d)",cex=1.5,col="white")




#image(stoc.pert$meshpts,stoc.pert$meshpts,t(prop.gyx),main="Stochastic - sd",xlab=xlab.text,labcex = 0.8,ylab=xlab.text)
#contour(stoc.pert$meshpts,stoc.pert$meshpts,t(prop.gyx),add=T)

#savePlot("ElasticityImagesFunctionsgxy.pdf",type="pdf"); 
#savePlot("ElasticityImagesFunctionsgxy.bmp",type="bmp"); 
#savePlot("ElasticityImagesFunctionsgxy.png",type="png"); 


##############################################
quartz()
par(mfrow=c(1,2),col="grey",pty="s")
hist(stoc.pert$store.sens.s.fn.v2.0.19,col="grey",main="")
hist(log10(stoc.pert$store.sens.s.fn.v2.0.19),col="grey",main="")



################ SPE: plot kernel std. dev 
meshpts<-stoc.pert$meshpts
image(meshpts,meshpts,t(stoc.pert$kernel.sd),xlab="Longest leaf length (log scale) t", 
    ylab="Longest leaf length (log scale) t+1",col=grey(100:300/300))
title(main="Kernel SD")
contour(meshpts,meshpts,t(stoc.pert$kernel.sd),method = "flattest",add=TRUE,
vfont = c("sans serif", "plain"),labcex = 0.8,levels=c(0.05,seq(0,1,by=0.2))) 


############## SPE: modify scaling of the elasticity contour plots 
 
mean.K<-det.pert.analysis(stoc.pert$mean.K); 
sum(mean.K$elas)

quartz()

par(mfrow=c(2,2),pty="s",cex.lab=1.2,cex.axis=1.2,xaxt="n",yaxt="n")
meshpts<-stoc.pert$meshpts; h=meshpts[2]-meshpts[1]; 
hinv=1/h; hi2=hinv*hinv; 

matrix.image.plot(x=meshpts,y=meshpts,hi2*(mean.K$elas),xlab=xlab.text,ylab=ylab.text,col=grey(100:400/400),add.contour=TRUE,nlevels=5)
#title(main="Mean K elasticities")
#contour(meshpts,meshpts,hi2*t(mean.K$elas),method = "flattest",add=TRUE,vfont = c("sans serif", "plain"),labcex = 0.8,nlevels=5) 
par(xaxt="s",yaxt="s");
axis(side=1,at=log(c(5,10,50)),labels=c("5","10","50"))
axis(side=2,at=log(c(5,10,50)),labels=c("5","10","50"))
#text(2,2,"a)",col="white",cex=1.5)

par(pty="s",pch=19,cex.lab=1.2,cex.axis=1.2,xaxt="n",yaxt="n")

matrix.image.plot(x=meshpts,y=meshpts,hi2*(stoc.pert$elas.s),xlab=xlab.text,ylab=ylab.text,col=grey(100:400/400),add.contour=TRUE,nlevels=5)
#title(main="Stochastic elasticities")
#contour(meshpts,meshpts,hi2*t(stoc.pert$elas.s),method = "flattest",add=TRUE,vfont = c("sans serif", "plain"),labcex = 0.8,nlevels=5) 
par(xaxt="s",yaxt="s");
axis(side=1,at=log(c(5,10,50)),labels=c("5","10","50"))
axis(side=2,at=log(c(5,10,50)),labels=c("5","10","50"))
#text(2,2,"b)",col="white",cex=1.5)

par(pty="s",pch=19,cex.lab=1.2,cex.axis=1.2,xaxt="n",yaxt="n")

matrix.image.plot(x=meshpts,y=meshpts,hi2*(stoc.pert$elas.s.mean),xlab=xlab.text,ylab=ylab.text,col=grey(100:400/400),add.contour=TRUE,nlevels=5)
#title(main="Stochastic elasticities, mean")
#contour(meshpts,meshpts,hi2*t(stoc.pert$elas.s.mean),method = "flattest",add=TRUE,vfont = c("sans serif", "plain"),labcex = 0.8,nlevels=5) 
par(xaxt="s",yaxt="s");
axis(side=1,at=log(c(5,10,50)),labels=c("5","10","50"))
axis(side=2,at=log(c(5,10,50)),labels=c("5","10","50"))
#text(2,2,"c)",col="white",cex=1.5)

par(pty="s",pch=19,cex.lab=1.2,cex.axis=1.2,xaxt="n",yaxt="n")

matrix.image.plot(x=meshpts,y=meshpts,hi2*(stoc.pert$elas.s.sd),xlab=xlab.text,ylab=ylab.text,col=grey(80:300/300),add.contour=TRUE,nlevels=5)
#title(main="Stochastic elasticities, sd")
#contour(meshpts,meshpts,hi2*t(stoc.pert$elas.s.sd),method = "flattest",add=TRUE,vfont = c("sans serif", "plain"),labcex = 0.8,nlevels=5) 
par(xaxt="s",yaxt="s");
axis(side=1,at=log(c(5,10,50)),labels=c("5","10","50"))
axis(side=2,at=log(c(5,10,50)),labels=c("5","10","50"))
#text(2,2,"d)",cex=1.5)
#abline(v=1.45); abline(h=5.17); 

#savePlot("ElasticityImages.pdf",type="pdf"); 
#savePlot("ElasticityImages.bmp",type="bmp"); 
#savePlot("ElasticityImages.png",type="png"); 



cor(as.numeric(mean.K$elas),as.numeric(stoc.pert$elas.s))

cor(as.numeric(mean.K$elas),as.numeric(stoc.pert$elas.s.mean))

cor(as.numeric(mean.K$elas),as.numeric(elas.s.sd))

cor(as.numeric(stoc.pert$elas.s),as.numeric(stoc.pert$elas.s.mean))

cor(as.numeric(stoc.pert$elas.s),as.numeric(elas.s.sd))

cor(as.numeric(stoc.pert$elas.s.mean),as.numeric(elas.s.sd))

plot(stoc.pert$elas.s.mean,elas.s.sd)
abline(0,-1)
abline(1,-1)


quartz()

par(mfrow=c(1,3),pty="s")

image(meshpts,meshpts,t(log(mean.K$sens)),xlab="Longest leaf length (log scale) t", ylab="Longest leaf length (log scale) t+1")

contour(meshpts,meshpts,t(log(mean.K$sens)),method = "flattest",add=TRUE,vfont = c("sans serif", "plain"),labcex = 0.8) 

#text(locator(1),"a)")

image(meshpts,meshpts,t(log(stoc.pert$sens.s)),xlab="Longest leaf length (log scale) t", ylab="Longest leaf length (log scale) t+1")

contour(meshpts,meshpts,t(log(stoc.pert$sens.s)),method = "flattest",add=TRUE,vfont = c("sans serif", "plain"),labcex = 0.8) 

image(meshpts,meshpts,log(abs(t(stoc.pert$sens.s.v2.0))),xlab="Longest leaf length (log scale) t", ylab="Longest leaf length (log scale) t+1")

contour(meshpts,meshpts,log(abs(t(stoc.pert$sens.s.v2.0))),method = "flattest",add=TRUE,vfont = c("sans serif", "plain"),labcex = 0.8) 

##################################################################################################################################
##Elasticity analysis to lower level parameters



n.big.matrix=100
n.est<-10000
n.runin<-500

pert.p=stoc.pert.analysis.lower.p(n.est,11,0.05)


sds=sd.vec

sds[5]=sqrt(VarCovar.grow.rec[1,1])
sds[10]=sqrt(VarCovar.grow.rec[2,2])

elas.s.p=pert.p$elas.s.p

round(elas.s.p,2)

elas.s.p.mean=pert.p$elas.s.p.mean

round(elas.s.p.mean,2)

elas.s.p.sd=pert.p$elas.s.p-pert.p$elas.s.p.mean

round(elas.s.p.sd,2)

sen.s.p.mean=pert.p$sens.s.p

round(sen.s.p.mean,2)

sen.s.p.sd=(pert.p$lam.stoch/sds)*elas.s.p.sd

sen.s.p.v2=(1/(2*sds))*sen.s.p.sd

sen.s.p.v2.0=pert.p$sens.s.p.v2.0

sen.s.p.v2.all=sen.s.p.v2

sen.s.p.v2.all[sds==0]=sen.s.p.v2.0[sds==0]

round(sen.s.p.v2.all,2)



elas.sens=data.frame(elas.s.p,elas.s.p.mean,abs(elas.s.p.sd),sen.s.p.mean,abs(sen.s.p.v2.all))

#pairs(elas.sens)

#sen.s.p.mean=(pert.p$lam.stoch/p.vec)*pert.p$elas.s.p.mean

fit=lm(log(abs(sen.s.p.v2.all))~log(sen.s.p.mean))

summary(fit)

fit=lm(log(abs(sen.s.p.v2.all))~log(sen.s.p.mean),subset=sen.s.p.mean>0.3)

summary(fit)

fit2=lm(log(abs(sen.s.p.v2.all))~offset(log(1/2)+2*log(sen.s.p.mean))-1,subset=sen.s.p.mean>0.3)

anova(fit,fit2)


quartz()

par(mfrow=c(1,2),bty="l",pch=19,pty="s")

plot(sen.s.p.mean,sen.s.p.v2.all,xlab=expression("Sensitivity to mean "*partialdiff*lambda [S]/partialdiff*mu*(theta[i])),
ylab=expression("Sensitivity to variance "~~partialdiff*lambda [S]/partialdiff*sigma^2*(theta[i])),type="n")

points(sen.s.p.mean[sds>0],sen.s.p.v2.all[sds>0],pch=19)
points(sen.s.p.mean[sds==0],sen.s.p.v2.all[sds==0],pch=21)

text(locator(1),"a)")

#abline(log(1/2),2)

#s.mean=seq(min(sen.s.p.mean),max(sen.s.p.mean),length=100)
#s.var=-s.mean-0.5*s.mean*s.mean

#points(s.mean,abs(s.var),type="l")

#quartz()

#par(bty="l",pch=19,pty="s")

#plot(pert.p$sens.s.p.v2.0,pert.p$sens.s.p.v2.0.old)

#abline(0,1)

#quartz()

#par(bty="l",pch=19,pty="s")

sens.minus.non.lin=sen.s.p.v2.all-pert.p$t2.s.p

alt.var.est=2*sen.s.p.v2.all/pert.p$lam.stoch-pert.p$t2.s.p+sen.s.p.mean*sen.s.p.mean/(pert.p$lam.stoch*pert.p$lam.stoch)

approx.result=-(2*sen.s.p.v2.all/pert.p$lam.stoch-pert.p$t2.s.p+pert.p$var.vt.partial.k.wt)

plot(sen.s.p.mean,approx.result,xlab=expression("Sensitivity to mean "*partialdiff*lambda [S]/partialdiff*mu*(theta[i])),
ylab=expression(-"("*partialdiff*lambda [S]/partialdiff*sigma^2*(theta[i])*"-NLA-Var("*partialdiff*lambda/partialdiff*theta[i]*"))"),log="xy",type="n")

points(sen.s.p.mean[sds>0],approx.result[sds>0],pch=19)

points(sen.s.p.mean[sds==0],approx.result[sds==0],pch=21)

text(locator(1),"b)")

abline(log(1/(pert.p$lam.stoch*pert.p$lam.stoch)),2)

#abline(0,1)

fit=lm(log(abs(sens.minus.non.lin))~log(sen.s.p.mean))

summary(fit)

fit=lm(log(abs(sens.minus.non.lin))~log(sen.s.p.mean),subset=sen.s.p.mean>0.3)

summary(fit)

abline(fit)


plot(pert.p$var.vt.partial.k.wt,abs(alt.var.est),xlab=expression("Sensitivity to mean "*partialdiff*lambda [S]/partialdiff*mu*(theta[i])),
ylab=expression("|Sensitivity to variance - NLA term |"),log="xy",type="n")

points(pert.p$var.vt.partial.k.wt[sds>0],abs(alt.var.est[sds>0]),pch=19)


points(pert.p$var.vt.partial.k.wt[sds==0],abs(alt.var.est[sds==0]),pch=21)

text(locator(1),"c)")

abline(0,1)


#abline(log(1/2),2)

fit=lm(log(abs(sens.minus.non.lin))~log(sen.s.p.mean))

summary(fit)

fit=lm(log(abs(sens.minus.non.lin))~log(sen.s.p.mean),subset=sen.s.p.mean>0.3)

summary(fit)

abline(fit)



#savePlot("SensMeanVar.pdf",type="pdf"); 
#savePlot("SensMeanVar.bmp",type="bmp"); 
#savePlot("SensMeanVar.png",type="png"); 

abline(0,1)

R0.calc<-function(n,params){


	M<-bigmatrix(n,params);

	N<-solve(diag(n)-M$Pmatrix);
	R<- M$Bmatrix %*% N
	
	ave.R0<-Re(eigen(R)$values[1]);

	lam<-Re(eigen(M$matrix)$values[1]);
	return(list(lam=lam,ave.R0=ave.R0,T=log(ave.R0)/log(lam)))
}

p.est<-1;

R0.pest1<-R0.calc(n.big.matrix,p.vec)

p.est<-1/R0.pest1$ave.R0

K<-bigmatrix(n.big.matrix,p.vec);

det.sen.elas=det.pert.analysis(K$Pmatrix+p.est*K$Bmatrix)

vKw=det.sen.elas$v %*% ((K$Pmatrix+p.est*K$Bmatrix) %*% det.sen.elas$w)


delta=0.05
partial.K=array(0,dim=c(11,n.big.matrix,n.big.matrix))
partial.K2=array(0,dim=c(11,n.big.matrix,n.big.matrix))

for(i in 1:11){

	p.pert=rep(0,11)
	p.pert[i]=delta
			
	K.ups=bigmatrix(n.big.matrix,(1+p.pert)*store.p.vec);
	K.downs=bigmatrix(n.big.matrix,(1-p.pert)*store.p.vec);

	K.0=K$Pmatrix+p.est*K$Bmatrix
	K.up=K.ups$Pmatrix+p.est*K.ups$Bmatrix
	K.down=K.downs$Pmatrix+p.est*K.downs$Bmatrix

	partial.K[i,,]=(K.up-K.down)/(2*delta*p.vec[i])
		
	partial.K2[i,,]=(((K.up-K.0)/(delta*p.vec[i]))-
				   ((K.0-K.down)/(delta*p.vec[i])))/
				   (delta*p.vec[i])
}

v.partial.k.w=rep(0,11)
v.partial.k2.w=rep(0,11)

for(i in 1:11){ 

	v.partial.k.w[i]=det.sen.elas$v %*% (partial.K[i,,] %*% det.sen.elas$w)
	v.partial.k2.w[i]=det.sen.elas$v %*% (partial.K2[i,,] %*% det.sen.elas$w)

}

t.non.lin=v.partial.k2.w/vKw
t.var=(v.partial.k.w*v.partial.k.w)/(vKw*vKw)

quartz()

par(mfrow=c(1,3),bty="l",pch=19,pty="s")

plot(t.non.lin,pert.p$t2.s.p)

abline(0,1)

plot(t.var,pert.p$t1.s.p)

abline(0,1)

sens.sv=0.5*(t.non.lin-t.var)

plot(sens.sv,sen.s.p.v2.all)

abline(0,1)

quartz()

par(mfrow=c(1,3),bty="l",pch=19,pty="s")

plot(t.non.lin,sen.s.p.v2.all)

abline(0,1)

plot(t.var,sen.s.p.v2.all)

abline(0,1)

plot(t.var,abs(sens.minus.non.lin))

abline(0,1)

det.sen.elas$v %*% det.sen.elas$w

########################################################################################################################################
#LTRE code

n.big.matrix=30
n.est<-3000
n.runin<-500

sens=rep(NA,11)

sens.log.lam=rep(NA,11)

var.lam=rep(NA,11)

for(i in 1:11){

	get.stuff=LTRE.f(n.est,i,0.01)

	lam.mean.k=Re(eigen(get.stuff$mean.K)$values[1])

	w<-Re(eigen(get.stuff$mean.K)$vectors[,1]); 
	v<-Re(eigen(t(get.stuff$mean.K))$vectors[,1]);
	vw<-sum(v*w);

	sens[i]= (v %*% (get.stuff$partial.K %*% w))/vw

	sens.log.lam[i]= (1/lam.mean.k)*(v %*% (get.stuff$partial.K %*% w))/vw

	var.lam[i]=get.stuff$var.lam

	cat(i,"\n")

}

contribs=rep(NA,12)

contribs.log.lam=rep(NA,12)


for(i in 1:11){
	contribs[i]=sens[i]*sens[i]*sd.vec[i]*sd.vec[i]

	if(i==5) contribs[i]=sens[i]*sens[i]*VarCovar.grow.rec[1,1]

	if(i==10) contribs[i]=sens[i]*sens[i]*VarCovar.grow.rec[2,2]

	contribs.log.lam[i]=sens.log.lam[i]*sens.log.lam[i]*sd.vec[i]*sd.vec[i]

	if(i==5) contribs.log.lam[i]=sens.log.lam[i]*sens.log.lam[i]*VarCovar.grow.rec[1,1]

	if(i==10) contribs.log.lam[i]=sens.log.lam[i]*sens.log.lam[i]*VarCovar.grow.rec[2,2]

}

contribs[12]=2*sens[5]*sens[10]*VarCovar.grow.rec[1,2]

contribs.log.lam[12]=2*sens.log.lam[5]*sens.log.lam[10]*VarCovar.grow.rec[1,2]





for(i in 1:12){
	
	if (i==12) cat("CoVar(",p.vec.names[5],p.vec.names[10],") contribution",contribs[12],"\n") else
	cat("Var(",p.vec.names[i],") contribution",contribs[i],"\n")

}

sum(contribs)

mean(var.lam)

##########################################################################################
#test 

#generate parameter set

n.big.matrix=50
n.est<-2000
n.runin<-500

iter<-iterate.model(p.vec,n.est,store.env=T,DD=F)

lams=rep(NA,n.est)

for (i in 1:n.est) {
	
		K<-bigmatrix(n.big.matrix,iter$params.yr[i,]);

		Ki=K$Pmatrix+iter$p.est.yr[i]*K$Bmatrix

		lams[i]=Re(eigen(Ki)$values[1])

	}

params=iter$params.yr

################ Linear analysis on lam


fit.lm=lm(lams~params[,1]+params[,2]+params[,3]+params[,5]+params[,6]+params[,10])
summary(fit.lm); par(mfrow=c(2,2)); plot(fit.lm);

# get the linear terms separately 
terms=predict(fit.lm,type="terms"); 
colnames(terms)<-c("p1","p2","p3","p5","p6","p10"); 

## Do the variance decomposition, including Cov(p5,p10) and the errors  
vlam.terms=c(apply(terms,2,var), 2*cov(terms[,4],terms[,6]), mean(fit.lm$residuals^2)); 
names(vlam.terms)<-c(colnames(terms),"C5,10","Error");  

sum(vlam.terms)/var(lams); # should be close to 1
vlam.terms=vlam.terms; 
round(vlam.terms,digits=2); 

################ Linear models adding interaction terms one at a time.  
residual.var=mean((fit.lm$residuals)^2); 
variable.indices=c(1,2,3,5,6,10); 

### add pairwise terms one at a time, check for how much additional variance they explain 
for(i in 1:5) {
for(j in (i+1):6) {
i1=variable.indices[i]; j1=variable.indices[j]; 
fit.gam.ij=lm(lams~params[,1]+params[,2]+params[,3]+params[,5]+params[,6]+params[,10]+params[,i1]:params[,j1])
residual.vij=mean((fit.gam.ij$residuals)^2);  
cat(i1,j1,(residual.var-residual.vij)/residual.var,"\n"); 
}}

#bigest contribution ~17% so go with independent terms model


lin.var.est=round(vlam.terms,digits=2)


################ Linear analysis on log lam


fit.lm=lm(log(lams)~params[,1]+params[,2]+params[,3]+params[,5]+params[,6]+params[,10])
summary(fit.lm); par(mfrow=c(2,2)); plot(fit.lm);

# get the linear terms separately 
terms=predict(fit.lm,type="terms"); 
colnames(terms)<-c("p1","p2","p3","p5","p6","p10"); 

## Do the variance decomposition, including Cov(p5,p10) and the errors  
vlam.terms=c(apply(terms,2,var), 2*cov(terms[,4],terms[,6]), mean(fit.lm$residuals^2)); 
names(vlam.terms)<-c(colnames(terms),"C5,10","Error");  

sum(vlam.terms)/var(log(lams)); # should be close to 1
vlam.terms=vlam.terms; 
round(vlam.terms,digits=2); 

################ Linear models adding interaction terms one at a time.  
residual.var=mean((fit.lm$residuals)^2); 
variable.indices=c(1,2,3,5,6,10); 

### add pairwise terms one at a time, check for how much additional variance they explain 
for(i in 1:5) {
for(j in (i+1):6) {
i1=variable.indices[i]; j1=variable.indices[j]; 
fit.gam.ij=lm(log(lams)~params[,1]+params[,2]+params[,3]+params[,5]+params[,6]+params[,10]+params[,i1]:params[,j1])
residual.vij=mean((fit.gam.ij$residuals)^2);  
cat(i1,j1,(residual.var-residual.vij)/residual.var,"\n"); 
}}

#Conclusion: interactions between params 1 and 2 explains 30% of the remaining variance 

###### Linear model with independent and interaction for parameters 1 and 2 
fit12 = lm(log(lams)~params[,1]*params[,2]+params[,3]+params[,5]+params[,6]+params[,10])

terms=predict(fit12,type="terms");
terms=terms[,c(1,2,7,3,4,5,6)] #re-order terms as lm and gam do things differently
colnames(terms)<-c("p1","p2","p1p2","p3","p5","p6","p10"); 

### make the interaction term uncorrelated with the separate terms for 1 and 2 
fitb=lm(terms[,3]~terms[,1]+terms[,2]); 
terms[,3]=terms[,3]-terms[,1:2]%*%fitb$coef[2:3]; 
terms[,1]=terms[,1]*(1+fitb$coef[2]); 
terms[,2]=terms[,2]*(1+fitb$coef[3]); 

## Do the variance decomposition, including Cov(p5,p10) and the errors  
vlam12.terms=c(apply(terms,2,var), 2*cov(terms[,4],terms[,6]), mean(fit12$residuals^2)); 
names(vlam12.terms)<-c(colnames(terms),"C5,10","Error");  

sum(vlam12.terms)/var(log(lams)); 
round(cor(terms),digits=2); #should be near zero except on the diagonal and 5,10 correlations 

vlam12.terms=vlam12.terms; 
round(vlam12.terms,digits=2); 

lin.log.var.est=round(vlam12.terms,digits=2)

################ GAM analysis on log lam


fit=gam(log(lams)~s(params[,1])+s(params[,2])+s(params[,3])+s(params[,5])+s(params[,6])+s(params[,10]))
summary(fit); gam.check(fit);

# get the spline terms separately 
terms=predict(fit,type="terms"); 
colnames(terms)<-c("p1","p2","p3","p5","p6","p10"); 

## Do the variance decomposition, including Cov(p5,p10) and the errors  
vlam.terms=c(apply(terms,2,var), 2*cov(terms[,4],terms[,6]), mean(fit$residuals^2)); 
names(vlam.terms)<-c(colnames(terms),"C5,10","Error");  

sum(vlam.terms)/var(log(lams)); # should be close to 1
vlam.terms=vlam.terms; 
round(vlam.terms,digits=2); 


################ GAMs adding bivariate terms one at a time.  
residual.var=mean((fit$residuals)^2); 
variable.indices=c(1,2,3,5,6,10); 

### add pairwise terms one at a time, check for how much additional variance they explain 
for(i in 1:5) {
for(j in (i+1):6) {
i1=variable.indices[i]; j1=variable.indices[j]; 
fit.gam.ij=gam(log(lams)~s(params[,1])+s(params[,2])+s(params[,3])+s(params[,5])+s(params[,6])+s(params[,10])+s(params[,i1],params[,j1]))
residual.vij=mean((fit.gam.ij$residuals)^2);  
cat(i1,j1,(residual.var-residual.vij)/residual.var,"\n"); 
}}
# Conclusion: interactions between params 1 and 2 explains 70% of the remaining variance 

###### GAM analysis with independent and joint splines for parameters 1 and 2 
fit12 = gam(log(lams)~s(params[,1])+s(params[,2])+s(params[,1],params[,2])+s(params[,3])+s(params[,5])+s(params[,6])+s(params[,10]))
gam.check(fit12); 

terms=predict(fit12,type="terms"); 
colnames(terms)<-c("p1","p2","p1p2","p3","p5","p6","p10"); 

### make the bivariate(1,2) spline term uncorrelated with the separate terms for 1 and 2 
fitb=lm(terms[,3]~terms[,1]+terms[,2]); 
terms[,3]=terms[,3]-terms[,1:2]%*%fitb$coef[2:3]; 
terms[,1]=terms[,1]*(1+fitb$coef[2]); 
terms[,2]=terms[,2]*(1+fitb$coef[3]); 

## Do the variance decomposition, including Cov(p5,p10) and the errors  
vlam12.terms=c(apply(terms,2,var), 2*cov(terms[,5],terms[,7]), mean(fit12$residuals^2)); 
names(vlam12.terms)<-c(colnames(terms),"C5,10","Error");  

sum(vlam12.terms)/var(log(lams)); 
round(cor(terms),digits=2); #should be near zero except on the diagonal and 5,10 correlations 

vlam12.terms=vlam12.terms; 
round(vlam12.terms,digits=2); 

gam.var.est=round(vlam12.terms,digits=2)

lin.var.est

gam.var.est




########################################################################################################################################
#ESS stuff, some functions then code for analysis


invader.gr<-function(params.yr,p.est.yr,index=NULL,I.params=NULL){

	nt<-rep(1/n.big.matrix,n.big.matrix)
	Rt<-rep(NA,n.est)

	for (year.t in 1:n.est){

		params.t<-params.yr[year.t,]

		if (length(index)>0) params.t[index]<-I.params+params.t[index]
		
		year.K<-bigmatrix(n.big.matrix,params.t)

		p.est<-p.est.yr[year.t]

		nt1<-(year.K$Pmatrix+p.est*year.K$Bmatrix) %*% nt
		
		sum.nt1<-sum(nt1)
		
		Rt[year.t]<-log(sum.nt1)

		nt<-nt1/sum.nt1
		
		#cat(Rt[year.t],"\n")

	}

return(mean(Rt[n.runin:n.est],na.rm=T))
}

I.gr.general<-function(I.params,index,R.env) {

	rate<-invader.gr(R.env$params.yr,R.env$p.est.yr,index,I.params)

	cat(I.params,"\tin function... rate =",rate,"\n");
	return(-rate);

}

### ESS function using one of a variety of optimisation methods ###

find.ESS<-function(p.vec,index,optim.type,scaling=p.vec[index]) {

	converge<-0;
	n.iterations<-0;
	old.params<-rep(0,length(index));
	params<-p.vec[index];

	repeat {
		### Generate resident environment ###
		R.env<-iterate.model(p.vec,n.est,store.env=T,DD=T)
		for (i in index){
			R.env$params.yr[,i]<-R.env$params.yr[,i]-p.vec[i]
		}
		### Find fittest strategy in current resident environment using... ###
		I.best.fitness<-switch(optim.type,
			### 1) ...quasi-Newton method ###
			optim(params,I.gr.general,method="BFGS",
				  control=list(fnscale=2,parscale=scaling,maxit=20),
				  index=index,R.env=R.env),
			### 2) ...Nelder-Mead simplex method ###
			optim(params,I.gr.general,method="Nelder-Mead",
				  control=list(fnscale=2,parscale=scaling,alpha=1.0,beta=0.5,gamma=2.0),
				  index=index,R.env=R.env),
			### 3) ...simulated annealing method ###
			optim(params,I.gr.general,method="SANN",
				  control=list(fnscale=2,parscale=scaling,temp=10),
				  index=index,R.env=R.env),
			### 4) ...one parameter optimiser ###
			optimise(f=I.gr.general,lower=2*params,upper=0.5*params,tol=0.01,
				  index=index,R.env=R.env)
		)
		### Set best invader strategy as new resident ###
		if (optim.type==4) params=I.best.fitness$minimum else params<-I.best.fitness$par;
		
		p.vec[index]<-params;
		cat(params,"\t",I.best.fitness$counts,"\t",I.best.fitness$convergence,"\n\n");
		cat(I.best.fitness$value,"\n");	
		### Exit loop when ESS is found or convergence is too slow ###
		if (sum(signif(old.params,3)==signif(params,3))==length(index)) converge<-1;
		#if(abs(I.best.fitness$objective)<0.001) converge<-1;
		n.iterations<-n.iterations+1;
		old.params<-params;
		if (converge==1 || n.iterations==10) {
			cat("Number of invasions =",n.iterations-1,", ESS parameters =",params,"\n");
			break;
		}
	}

	return(params);

}

##################################################################################################################################
##ESS analysis

n.big.matrix=50
n.est<-5000
n.runin<-500

#intercept and slope of flowering function

p.vec=store.p.vec

p.vec[3]=-400

p.vec[4]=100

index=c(3,4)

#p.vec.ESS.2p<-find.ESS(p.vec,index,4)

p.vec.ESS.2p=c(-200*3.61,200)

p.vec[index]<-p.vec.ESS.2p

iter<-iterate.model(p.vec,n.est,store.env=T,DD=T)

#calculate stochastic growth rate

mean(iter$Rt[n.runin:n.est],na.rm=T)

plot(iter$N.year,type="b",pch=19)

mean(iter$N.year[n.runin:n.est])

p.size.dist.ESS.2p<-iter$size.dist/sum(iter$size.dist)

mean.size<-sum(p.size.dist.ESS.2p*exp(iter$meshpts))
mean.size

p.size.dist.fl.ESS.2p<-iter$size.dist.fl/sum(iter$size.dist.fl)

mean.size.fl<-sum(p.size.dist.fl.ESS.2p*exp(iter$meshpts))
mean.size.fl

diff<-iter$meshpts[2]-iter$meshpts[1]

quartz()

par(mfrow=c(1,2),bty="l",pty="s")

hist(lst,breaks=10,main="",col="grey",ylim=c(0,420),xlab="Size")
lines(iter$meshpts,p.size.dist*599/sum(p.size.dist*diff))
#lines(iter$meshpts,p.size.dist.cf*599/sum(p.size.dist.cf*diff),lty=2)
lines(iter$meshpts,p.size.dist.ESS.2p*599/sum(p.size.dist.ESS.2p*diff),lty=3)

#text(locator(1),"a)")

hist(size.fl,breaks=10,main="",col="grey",xlab="Flowering size",ylim=c(0,75))
lines(iter$meshpts,p.size.dist.fl*21.4/sum(p.size.dist.fl*diff))
#lines(iter$meshpts,p.size.dist.fl.cf*21.4/sum(p.size.dist.fl.cf*diff),lty=2)
lines(iter$meshpts,p.size.dist.fl.ESS.2p*21.4/sum(p.size.dist.fl.ESS.2p*diff),lty=3)


p.vec=store.p.vec

#now just intercept

index=c(3)

#p.vec.ESS<-find.ESS(p.vec,index,4)

p.vec.ESS=-14.15

p.vec[index]<-p.vec.ESS

iter<-iterate.model(p.vec,n.est,store.env=T,DD=T)

#calculate stochastic growth rate

mean(iter$Rt[n.runin:n.est],na.rm=T)

plot(iter$N.year,type="b",pch=19)

mean(iter$N.year[n.runin:n.est])

p.size.dist.ESS<-iter$size.dist/sum(iter$size.dist)

mean.size<-sum(p.size.dist.ESS*exp(iter$meshpts))
mean.size

p.size.dist.fl.ESS<-iter$size.dist.fl/sum(iter$size.dist.fl)

mean.size.fl<-sum(p.size.dist.fl.ESS*exp(iter$meshpts))
mean.size.fl

diff<-iter$meshpts[2]-iter$meshpts[1]

quartz()

par(mfrow=c(1,2),bty="l",pty="s")

hist(lst,breaks=10,main="",col="grey",ylim=c(0,420),xlab="Size")
lines(iter$meshpts,p.size.dist*599/sum(p.size.dist*diff))
lines(iter$meshpts,p.size.dist.ESS*599/sum(p.size.dist.ESS*diff))
lines(iter$meshpts,p.size.dist.ESS.2p*599/sum(p.size.dist.ESS.2p*diff),lty=3)

text(locator(1),"a)")

hist(size.fl,breaks=10,main="",col="grey",xlab="Flowering size",ylim=c(0,40))
lines(iter$meshpts,p.size.dist.fl*21.4/sum(p.size.dist.fl*diff))
lines(iter$meshpts,p.size.dist.fl.ESS*21.4/sum(p.size.dist.fl.ESS*diff))
lines(iter$meshpts,p.size.dist.fl.ESS.2p*21.4/sum(p.size.dist.fl.ESS.2p*diff),lty=3)

text(locator(1),"b)")

#savePlot("ESSImages.pdf",type="pdf"); 
#savePlot("ESSImages.bmp",type="bmp"); 
#savePlot("ESSImages.png",type="png"); 


#########################################################################################
#fitness landscapes - assumes ESS environment stored in iter

store.iter=iter

flow.inter<-seq(-30,-1,1)
n.test<-length(flow.inter)
Ls<-rep(NA,n.test)
mean.size.fl.I<-rep(NA,n.test)
var.size.fl.I<-rep(NA,n.test)

p.size.dist.fl.tmp<-matrix(NA,nrow=n.big.matrix,ncol=n.test)


iter$params.yr[,index]<-iter$params.yr[,index]-p.vec.ESS

for (i in 1:n.test){
	Ls[i]<-invader.gr(iter$params.yr,iter$p.est.yr,index=c(3),I.params=flow.inter[i])
	p.vec[3]<-flow.inter[i]
	i.iter<-iterate.model(p.vec,n.est,store.env=F,DD=T)
	p.size.dist.fl.tmp[,i]<-i.iter$size.dist.fl/sum(i.iter$size.dist.fl)
	mean.size.fl.I[i]<-sum(p.size.dist.fl.tmp[,i]*exp(i.iter$meshpts))
	var.size.fl.I[i]<-sum(p.size.dist.fl.tmp[,i]*exp(2*i.iter$meshpts))-mean.size.fl.I[i]^2

	cat(i,"  ",flow.inter[i],"  ",exp(Ls[i]),"  ",mean.size.fl.I[i],"  ",var.size.fl.I[i],"\n")
}

range(mean.size.fl.I[exp(Ls)>(max(exp(Ls))-0.01)])
range(mean.size.fl.I[exp(Ls)>(max(exp(Ls))-0.001)])

quartz()

par(mfrow=c(1,2),bty="l",pty="s",pch=19)

plot(flow.inter,Ls,type="n",xlab=expression("Intercept of flowering function " * italic(beta * scriptstyle(0))),
	ylab=expression("Fitness log("*lambda[s]*")"),ylim=c(-0.8,0.05))


SE<-fit.flow.r$coef.sd[1]

x.corners=as.numeric(c(fit.flow.r$coefficients[1]-2*SE,fit.flow.r$coefficients[1]-2*SE,fit.flow.r$coefficients[1]+2*SE,fit.flow.r$coefficients[1]+2*SE))

polygon(x.corners,c(-0.8,0,0,-0.8),col="grey90",border=0)

points(flow.inter,Ls,type="b")

abline(h=max(Ls))
#abline(h=max(Ls)-0.01,col="red")

abline(v=fit.flow.r$coefficients[1])

text(locator(1),"a)")

plot(mean.size.fl.I,Ls,type="n",xlab="Mean flowering size mm",
	ylab=expression("Fitness log("*lambda[s]*")"),ylim=c(-0.5,0.05))


SE<-sd(exp(size.fl))/sqrt(length(size.fl))

polygon(c(mean(exp(size.fl))-2*SE,mean(exp(size.fl))-2*SE,mean(exp(size.fl))+2*SE,mean(exp(size.fl))+2*SE),c(-0.8,0,0,-0.8),col="grey90",border=0)

points(mean.size.fl.I,Ls,type="b")

abline(h=max(Ls))
#abline(h=max(Ls)-0.01,col="red")

abline(v=mean(exp(size.fl)))

text(locator(1),"b)")

p.vec[3]<-p.vec.ESS

#savePlot("Fitness landscapes.pdf",type="pdf"); 
#savePlot("Fitness landscapes.bmp",type="bmp"); 
#savePlot("Fitness landscapes.png",type="png"); 


##############################################################################################
#Fitness landscape for varying amounts of variability in flowering intercepts

iter=store.iter

alpha=seq(0.0,5,length=21)
Ls.var=rep(NA,21)

actual.var=iter$params.yr[,index]-p.vec.ESS

for(i in 1:length(alpha)){
	iter$params.yr[,index]=p.vec.ESS+actual.var*alpha[i]
	Ls.var[i]<-invader.gr(iter$params.yr,iter$p.est.yr)
	cat(alpha[i],"  ",exp(Ls.var[i]),"\n")
}


##############################################################################################
#Fitness landscape for varying flowering slopes

iter=store.iter

bs=seq(3.5,4.5,length=21)

Ls.bs=rep(NA,length(bs))

iter$params.yr[,4]=0.0

for(i in 1:length(bs)){
	Ls.bs[i]<-invader.gr(iter$params.yr,iter$p.est.yr,index=c(4),I.params=bs[i])
	cat(bs[i],"  ",exp(Ls.bs[i]),"\n")
}


quartz()

par(mfrow=c(1,2),pch=19,bty="l",pty="s")

plot(var(actual.var)*alpha^2,exp(Ls.var),type="b",xlab=expression("Variance in flowering intercepts " * italic(beta * scriptstyle(0))),
	ylab=expression("Fitness "*lambda[s]))

abline(h=1)

abline(v=1.029^2)

#text(locator(1),"a)")


plot(bs,exp(Ls.bs),type="b",xlab=expression("Slope of flowering function " * italic(beta * scriptstyle(s))),ylim=c(0.95,1.002),
	ylab=expression("Fitness "*lambda[s]))

abline(h=1)

abline(v=3.861)

#text(locator(1),"b)")

#save(det.elas,stoc.elas,iter,p.size.dist.fl,p.size.dist.fl.ESS,p.size.dist.fl.ESS.2p,p.size.dist,p.size.dist.ESS,p.size.dist.ESS.2p,
#p.vec.ESS,p.vec.ESS.2p,flow.inter,Ls,mean.size.fl.I,file="c:\\temp\\ellner\\Carlina 100.Rdata")

############################################################################################
#effect of varying sd's

p.vec=store.p.vec

sd.vec=store.sd.vec

sd.index=1

n.ESSs=10

sd.ESSs=rep(NA,n.ESSs)

sds=seq(0,3*sd.vec[sd.index],length=n.ESSs)

index=c(3)

#now just intercept

for(i in 1:n.ESSs){

	sd.vec[sd.index]=sds[i]

	sd.ESSs[i]<-find.ESS(p.vec,index,4)

	cat(i,"  ",sds[i],"   ",sd.ESSs[i],"\n")

}

plot(sds,sd.ESSs)

###############################################################################################
#effect of varying covariance between 

p.vec=store.p.vec

sd.vec=store.sd.vec

limits=sqrt(VarCovar.grow.rec[1,1]*VarCovar.grow.rec[2,2])

n.ESSs=10

sd.ESSs=rep(NA,n.ESSs)

covs=seq(-1*limits,1*limits,length=n.ESSs)

index=c(3)

for(i in 1:n.ESSs){

	VarCovar.grow.rec[1,2]=VarCovar.grow.rec[2,1]=covs[i]

	sd.ESSs[i]<-find.ESS(p.vec,index,4)

	cat(i,"  ",covs[i],"\n")

}

plot(covs/limits,sd.ESSs)



############################################################################################
#effect of varying sd's

p.vec=store.p.vec

sd.vec=store.sd.vec

sd.index=1

n.ESSs=10

sd.ESSs=rep(NA,n.ESSs)

sds=seq(0,3*sd.vec[sd.index],length=n.ESSs)

index=c(3)

#now just intercept

for(i in 1:n.ESSs){

	sd.vec[sd.index]=sds[i]

	sd.lambda.s<-find.ESS(p.vec,index,4)

	cat(i,"  ",sds[i],"   ",sd.ESSs[i],"\n")

}

plot(sds,sd.ESSs)

































