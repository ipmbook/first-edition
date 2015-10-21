# rm(list=ls(all=TRUE))


############################################# 
# Kernel functions 
#############################################
fyx<-function(y,x,params) {
	return(params["f"]*exp(2*x)*dbeta(y+3,3,3));	
}

pyx<-function(y,x,params) {
	meanlog=log(exp(x)+params["g.rate"]); 
	return(params["stay"]*dnorm(y,meanlog,params["sigma.g"]))
}	

############################################################
# Find stable dist'n and lambda by iteration: deterministic
############################################################
ssd=function(Kmat,tol=1e-8,verbose=FALSE) {
	Kmat; qmax=1; lam=1; iter=1; 
	Nt=rep(1,nrow(Kmat));
	while (qmax>tol){
		Nt1=Kmat%*%Nt;
		qmax=sum(abs(Nt1-lam*Nt)); 
		lam=sum(Nt1); 
		if(verbose) cat(iter,sum(Nt),"  ",sum(Nt1),"  ",lam, qmax,"\n");
		iter=iter+1; 
	#set Nt equal to NT1 and normalise so Nt sums to 1
		Nt = Nt1/lam;
	}
	return(list(lambda=lam,stable.dist=Nt))
}

######################################################### 
# Compute spread rate(s): Deterministic model 
#########################################################
GoverS=function(s,P,F) {
	H = P + (1/(1-(s*L)^2))*F;
	g=log(ssd(H)$lambda);
	return(g/s); 
}


####################################################################
# 
Compute spread rate(s): Stochastic model, two year types 

####################################################################

GamOverS2=function(s,theta1,theta2,years=NULL,nyears=1000,nskip=250) {
	
if(is.null(years)) years=1+as.numeric(runif(nyears)<0.5); 
	
nyears=length(years); 
	
F=P=as.list(1:2); 
	
F[[1]]=h*outer(y,y,fyx,params=theta1);
	
P[[1]]=h*outer(y,y,pyx,params=theta1);
	
F[[2]]=h*outer(y,y,fyx,params=theta2);
	
P[[2]]=h*outer(y,y,pyx,params=theta2);
	
Kbar=(F[[1]]+F[[2]]+P[[1]]+P[[2]])/2;
	
n0=ssd(Kbar)$stable.dist; 
	Nt=n0/sum(n0); 
lam=numeric(nyears);
	
L1=theta1["L"]; L2=theta2["L"];
	
H=as.list(1:2);
	
H[[1]] = P[[1]] + (1/(1-(s*L1)^2))*F[[1]];
	
H[[2]] = P[[2]] + (1/(1-(s*L2)^2))*F[[2]];
	
for(j in 1:nyears){
		
yt=years[j];
		
Nt1=H[[yt]]%*%Nt;
		
lam[j]=sum(Nt1); 
		
Nt = Nt1/lam[j];
	
}
	
g=mean(log(lam[-(1:nskip)]))
	
return(g/s)
}	




####################################################################
# 
Compute spread rate(s): Stochastic model, uniform distribution

####################################################################

GamOverSU=function(s,theta1,theta2,years=NULL,nyears=1000,nskip=250) {

	cat(s,"\n")
	if(is.null(years)) years=runif(nyears)

	nyears=length(years);
 
	F=P=as.list(1:2);
  
	F[[1]]=h*outer(y,y,fyx,params=theta1);
  
	P[[1]]=h*outer(y,y,pyx,params=theta1);
  
	F[[2]]=h*outer(y,y,fyx,params=theta2);
  
	P[[2]]=h*outer(y,y,pyx,params=theta2);
  
	Kbar=(F[[1]]+F[[2]]+P[[1]]+P[[2]])/2;
  
	n0=ssd(Kbar)$stable.dist;
   
	Nt=n0/sum(n0); lam=numeric(nyears);
	
    L1=theta1["L"]; L2=theta2["L"];
	
    Ft=Pt=Ht=matrix(0,length(y),length(y))
   
		for(j in 1:nyears){
   
		yt=years[j];
   
		thetat=yt*theta1+(1-yt)*theta2

   		Lt=yt*L1+(1-yt)*L2;
   		
Ht=outer(y,y,pyx,params=thetat) + (1/(1-(s*Lt)^2))*outer(y,y,fyx,params=thetat)

   		Nt1=h*(Ht%*%Nt);

   		lam[j]=sum(Nt1);
   		Nt = Nt1/lam[j];


   		}

   	g=mean(log(lam[-(1:nskip)]))
   	
return(g/s)

   	}	



####################################################################
# Compute lambda_S: Stochastic model, two year types 
####################################################################
lambda.S=function(theta1,theta2,years=NULL,nyears=1000,nskip=250) {
	if(is.null(years)) years=1+as.numeric(runif(nyears)<0.5); 
	nyears=length(years); 
	F=P=as.list(1:2); 
	F[[1]]=h*outer(y,y,fyx,params=theta1);
	P[[1]]=h*outer(y,y,pyx,params=theta1);
	F[[2]]=h*outer(y,y,fyx,params=theta2);
	P[[2]]=h*outer(y,y,pyx,params=theta2);
	Kbar=(F[[1]]+F[[2]]+P[[1]]+P[[2]])/2;
	n0=ssd(Kbar)$stable.dist; 
	Nt=n0/sum(n0); lam=numeric(nyears);
	K=as.list(1:2);
	K[[1]] = P[[1]] + F[[1]];
	K[[2]] = P[[2]] + F[[2]];
	for(j in 1:nyears){
		yt=years[j];
		Nt1=K[[yt]]%*%Nt;
		lam[j]=sum(Nt1); 
		Nt = Nt1/lam[j];
	}
	lam.S=exp(mean(log(lam[-(1:nskip)])))
	return(lam.S)
}	

################################################################## 
# Kernel perturbation analysis of two-year-types model
# Assumes model iteration by midpoint rule
##################################################################
kernel.pert.analysis<-function(params1,params2,s,years,y){
	n.est=length(years);  h=y[2]-y[1];
	F1=h*outer(y,y,fyx,params=params1);
	P1=h*outer(y,y,pyx,params=params1);
	F2=h*outer(y,y,fyx,params=params2);
	P2=h*outer(y,y,pyx,params=params2);
	L1=params1["L"]; L2=params2["L"];
	H=as.list(1:2);
	H[[1]] = P1 + (1/(1-(s*L1)^2))*F1;
	H[[2]] = P2 + (1/(1-(s*L2)^2))*F2;
	M1 = (1/(1-(s*L1)^2));
	M2 = (1/(1-(s*L2)^2));
	mean.K = 0.5*(H[[1]] + H[[2]])
	mean.F = 0.5*(F1 + F2);
	mean.P = 0.5*(P1 + P2); 

### Get wt and Rt time series ###
	wt<-matrix(1/n.big.matrix, nrow=n.est+1, ncol=n.big.matrix);
	Rt<-rep(NA, n.est);
	for (i in 1:n.est) {
		Kt<-H[[years[i]]];
		wt[i+1,]<-Kt%*%wt[i,]
		Rt[i]<-sum(wt[i+1,]); 
		wt[i+1,]<-wt[i+1,]/Rt[i];
		if(i%%250==0) cat("wt and Rt ",i,"\n")
	}

### Get vt time series ###
	vt<-matrix(1/n.big.matrix, nrow=n.est+1, ncol=n.big.matrix);
	for (i in (n.est+1):2) {
		Kt<-H[[years[i-1]]];
		vt[i-1,]<-vt[i,]%*%Kt;
		vt[i-1,]<-vt[i-1,]/sum(vt[i-1,]);
		if(i%%250==0) cat("vt  ",i,"\n")
	}

############################################################## 
# Generate sensitivity and elasticity matrices for lambda_S with
# respect to kernel entries K(y,x), their mean and their sd 
##############################################################
	sens<-matrix(0,nrow=n.big.matrix,ncol=n.big.matrix);
	elas<-elas.mean<-elas.sd<-kernel.sd<-sens;
	
	for (i in 1:n.est) {
		#standard calculations needed for the various formulae
		Kt<-H[[years[i]]];
		vt1.wt=vt[i+1,]%*%t(wt[i,])
		Rt.vt1.wt1=(Rt[i]*t(vt[i+1,])%*%wt[i+1,])[1]; # this = <v_t+1,K_t w_t> 

		#calculation of the standard sensitivities and elasticities
		sens<-sens + vt1.wt/Rt.vt1.wt1;
		elas<-elas + Kt*(vt1.wt/Rt.vt1.wt1);
		elas.mean<-elas.mean + mean.K*(vt1.wt/Rt.vt1.wt1);
        elas.sd = elas.sd + (Kt-mean.K)*(vt1.wt/Rt.vt1.wt1); 
		kernel.sd=kernel.sd+(Kt-mean.K)^2	
		if(i%%250==0) cat("Elasticities and sensitivities ",i,"\n")
	}

	lam.stoch=exp(mean(log(Rt[251:n.est]),na.rm=T))
	sens<-lam.stoch*sens/n.est;
	elas<-elas/n.est;
	elas.mean=elas.mean/n.est;
	elas.sd = elas.sd/n.est; 
	kernel.sd=sqrt(kernel.sd/(n.est-1)); 

return(list(sens=sens, elas=elas, elas.mean=elas.mean, elas.sd=elas.sd, 
kernel.sd=kernel.sd, lam.S=lam.stoch, Rt=Rt, wt=wt, vt=vt))
}

################################################################## 
# Expectations for (Pbar,Fbar) perturbation analysis of two-year-types model
# Assumes model iteration by midpoint rule
##################################################################
entry.pert.analysis<-function(params1,params2,s,years,y,n.burnin=250){
	n.est=length(years);  h=y[2]-y[1];
	F1=h*outer(y,y,fyx,params=params1);
	P1=h*outer(y,y,pyx,params=params1);
	F2=h*outer(y,y,fyx,params=params2);
	P2=h*outer(y,y,pyx,params=params2);
	L1=params1["L"]; L2=params2["L"];
	H<-as.list(1:2);
	H[[1]] = P1 + (1/(1-(s*L1)^2))*F1;
	H[[2]] = P2 + (1/(1-(s*L2)^2))*F2;
	M1 = (1/(1-(s*L1)^2));
	M2 = (1/(1-(s*L2)^2));
	M=c(M1,M2); 
	mean.K = 0.5*(H[[1]] + H[[2]])
	mean.F = 0.5*(F1 + F2);
	mean.P = 0.5*(P1 + P2); 

### Get wt and Rt time series ###
	wt<-matrix(1/n.big.matrix, nrow=n.est+1, ncol=n.big.matrix);
	Rt<-rep(NA, n.est);
	for (i in 1:n.est) {
		Kt<-H[[years[i]]];
		wt[i+1,]<-Kt%*%wt[i,]
		Rt[i]<-sum(wt[i+1,]); 
		wt[i+1,]<-wt[i+1,]/Rt[i];
		if(i%%250==0) cat("wt and Rt ",i,"\n")
	}

### Get vt time series ###
	vt<-matrix(1/n.big.matrix, nrow=n.est+1, ncol=n.big.matrix);
	for (i in (n.est+1):2) {
		Kt<-H[[years[i-1]]];
		vt[i-1,]<-vt[i,]%*%Kt;
		vt[i-1,]<-vt[i-1,]/sum(vt[i-1,]);
		if(i%%250==0) cat("vt  ",i,"\n")
	}

############################################################## 
# Generate sensitivity and elasticity matrices for lambda_S with
# respect to kernel entries K(y,x), their mean and their sd 
##############################################################
	Evw<-EvwM<-matrix(0,nrow=n.big.matrix,ncol=n.big.matrix);

	j=0;
	for (i in n.burnin:(n.est-n.burnin)) {
		#standard calculations needed for the formulas
		Kt<-H[[years[i]]];
		vt1.wt=vt[i+1,]%*%t(wt[i,])
		Rt.vt1.wt1=(Rt[i]*t(vt[i+1,])%*%wt[i+1,])[1]; # this = <v_t+1,K_t w_t> 

		#expectations needed for P and F sensitivity analysis
		Evw<-Evw + vt1.wt/Rt.vt1.wt1;
		EvwM<-EvwM + M[years[i]]*(vt1.wt/Rt.vt1.wt1);
    	if(i%%250==0) cat("Elasticities and sensitivities ",i,"\n")
		j=j+1; 
	}

	lam.stoch=exp(mean(log(Rt[n.burnin:n.est]),na.rm=T))
	Evw<-Evw/j;
	EvwM<-EvwM/j;
	
return(list(Evw=Evw, EvwM=EvwM, lam.S=lam.stoch, gam.S=log(lam.stoch),Rt=Rt, wt=wt, vt=vt,mean.F=mean.F,mean.P=mean.P))
}



