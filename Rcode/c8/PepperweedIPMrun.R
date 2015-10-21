# rm(list=ls(all=TRUE))

setwd("/Users/stephenellner/Repos/ipm_book/Rcode/c8")
source("PepperweedIPMFunctions.R") 

setwd("/Users/stephenellner/Repos/ipm_book/c8/figures")

#############################################
# Iteration matrix setup 
#############################################
Minsize = -3.5; Maxsize=4.5; 
n.big.matrix = 150;

# boundary points b and mesh points y
h = (Maxsize-Minsize)/n.big.matrix;
b = Minsize+c(0:n.big.matrix)*h;
y = 0.5*(b[1:n.big.matrix]+b[2:(n.big.matrix+1)]);

################################################ 
# Deterministic model first 
################################################
L=15; g.rate=0.5; sigma.g=0.1; stay=0.94; f=.04; 
params=c(L,g.rate,sigma.g,f,stay); 
names(params)<-c("L","g.rate","sigma.g","f","stay"); 

### Growth rate and ssd of total population 
Fmat=h*outer(y,y,fyx,params=params);
Pmat=h*outer(y,y,pyx,params=params); 
out=ssd(Pmat+Fmat); lambda=out$lambda; 
w=out$stable.dist; w=w/(h*sum(w));
plot(y,w,type="o"); 
out$lambda;
B=h*sum(Fmat%*%w); B; # check that mean per-patch fecundity ~ 0.3

### Asymptotic spread rate 
cstar=optimize(GoverS,lower=1e-6,upper=(1-1e-6)/L,F=Fmat,P=Pmat)$objective; 
cstar; 

################################################################## 
# Stochastic model with equally common wet/dry years 
# Wet is good for dispersal, bad for growth & fecundity
##################################################################
params=c(L,g.rate,sigma.g,f,stay)
names(params)<-c("L","g.rate","sigma.g","f","stay"); 

## changes representing a big difference between wet and dry 
dL = L; dg=g.rate; ds=(1-stay); df=f;  

########## Case 1: only change demography, not dispersal distance		
nyears=10000; years=1+as.numeric(runif(nyears)<0.5); 
V=seq(0,0.9,length=11);
cstar=numeric(length(V));
for(j in 1:length(V)) {
	params1=params + V[j]*c(0,-dg,0,-df,-ds);  
	params2=params + V[j]*c(0,dg,0,df,ds) 
	cstar[j]=optimize(GamOverS2,lower=0.0001,upper=0.9999/params["L"],years=years,
				theta1=params1,theta2=params2)$objective;
	cat(j,cstar[j],"\n"); 
}
		
########### Case 2: only change dispersal distance, not demography 		
cstar2=numeric(length(V));
for(j in 1:length(V)) {
	params1=params + V[j]*c(dL,0,0,0,0);  
	params2=params + V[j]*c(-dL,0,0,0,0) 
	cstar2[j]=optimize(GamOverS2,lower=0.0001,upper=0.9999/params1["L"],years=years,
				theta1=params1,theta2=params2)$objective;
	cat(j,cstar2[j],"\n"); 
}

########### Case 3: change both, positive correlation
cstar3=numeric(length(V));
for(j in 1:length(V)) {
	params2=params + V[j]*c(-dL,-dg,0,-df,-ds);   
	params1=params + V[j]*c(dL,dg,0,df,ds); 
	cstar3[j]=optimize(GamOverS2,lower=0.0001,upper=0.9999/params1["L"],years=years,
				theta1=params1,theta2=params2)$objective;
	cat(j,cstar3[j],"\n"); 
}

########### Case 3b: change both, negative correlation
cstar3b=numeric(length(V));
for(j in 1:length(V)) {
	params1=params + V[j]*c(dL,-dg,0,-df,-ds);   
	params2=params + V[j]*c(-dL,dg,0, df,ds); 
	cstar3b[j]=optimize(GamOverS2,lower=0.0001,upper=0.9999/params1["L"],years=years,
				theta1=params1,theta2=params2)$objective;
	cat(j,cstar3[j],"\n"); 
}


dev.new(width=6,height=5); 
par(cex.axis=1.35,cex.lab=1.35,mgp=c(2.5,1,0),mar=c(5,4.25,2,1));
matplot(V,cbind(cstar,cstar2,cstar3,cstar3b),xlab="Fractional range of variation, V",ylab="Spread rate (m/yr)", 
	lty=c(1,1,2,3), lwd=2,col=c("black","gray60","gray60","gray60"),type="l",pch=16,bty="l"); 
legend(-0.06,9,legend=c("Vary demography","Vary dispersal","Vary both, C = 1", "Vary both, C = -1"),
col=c("black","gray60","gray60","gray60"),lty=c(1,1,2,3),lwd=2,bty="n",cex=1); 


#dev.copy2eps(file="PepperweedSpreadRate.eps"); 
#dev.copy2pdf(file="PepperweedSpreadRate.pdf"); 
