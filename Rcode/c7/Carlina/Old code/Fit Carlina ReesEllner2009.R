rm(list=ls(all=TRUE))
library(nlme)
library(MASS)
library(rjags)
library(coda)

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c7/Carlina",sep="")); 


dataf<-read.csv("Carlina survive or flower.fre",header=T)

attach(dataf)


s.dataf<-dataf[flow==1,]
s.dataf$flow<-0
names(s.dataf)[5]<-"die"

detach(dataf)

dataf<-read.csv("Carlina survive or died.txt",header=T)

attach(dataf)

dataf<-rbind(dataf,s.dataf)
attach(dataf)

rec.size<-z[age==1]
rec.yeart<-yeart[age==1]

fit.rec.lm<-lm(rec.size~factor(rec.yeart))
anova(fit.rec.lm)
summary(fit.rec.lm)


fit.rec.lm<-lm(rec.size~factor(rec.yeart)-1)


par(mfrow=c(1,1),bty="l",pch=19)
qqnorm(fit.rec.lm$coef)
qqline(fit.rec.lm$coef)

fit.rec.r<-lme(rec.size~1,data=dataf,random=~1|rec.yeart)


detach(dataf)

###################################################################################################################################
##growth analysis

dataf<-data.frame(read.csv("crlnagrw.txt",header=T));
attach(dataf)



fit.grow.i<-lm(lst1~factor(dataf$yeart)*lst-1)

summary(lm(fit.grow.i$coef[1:16]~fit.grow.i$coef[17:32]))

plot(fit.grow.i$coef[1:16],fit.grow.i$coef[17:32])

cor.test(fit.grow.i$coef[1:16],fit.grow.i$coef[17:32])

lst.c<-lst-mean(lst)


fit.grow<-lm(lst1~lst.c+factor(dataf$yeart))

fit.grow.i<-lm(lst1~lst.c*factor(dataf$yeart))

anova(fit.grow,fit.grow.i)

fit.grow<-lm(lst1~factor(dataf$yeart)*lst.c-1)

summary(lm(fit.grow$coef[1:16]~fit.grow$coef[17:32]))

plot(fit.grow$coef[1:16],fit.grow$coef[17:32])

cor.test(fit.grow$coef[1:16],fit.grow$coef[17:32])


par(mfrow=c(1,2),bty="l",pty="s",pch=19)

qqnorm(fit.grow$coef[1:16])
qqline(fit.grow$coef[1:16])
qqnorm(fit.grow$coef[17:32])
qqline(fit.grow$coef[17:32])

g.dataf<-groupedData(lst1~lst.c|yeart,data=data.frame(dataf,lst.c))

fit.grow.r<-lme(lst1~lst.c,data=g.dataf,random=~1|yeart)
fit.grow.r.s<-lme(lst1~lst.c,data=g.dataf,random=~lst.c|yeart)
fit.grow.r.s.diag<-lme(lst1~lst.c,data=g.dataf,random=pdDiag(~lst.c))
fit.grow.r.s.diag.v<-lme(lst1~lst.c,data=g.dataf,random=pdDiag(~lst.c),weight=varExp(form=~fitted(.)))

anova(fit.grow.r, fit.grow.r.s, fit.grow.r.s.diag, fit.grow.r.s.diag.v)

summary(fit.grow.r.s.diag)
intervals(fit.grow.r.s.diag)


par(mfrow=c(1,2),bty="l",pty="s",pch=19)

qqnorm(ranef(fit.grow.r.s.diag)[,1])
qqline(ranef(fit.grow.r.s.diag)[,1])
qqnorm(ranef(fit.grow.r.s.diag)[,2])
qqline(ranef(fit.grow.r.s.diag)[,2])


#####################################################################################################################################
#fit bivariate model to growth and recruitment

rec.yeart=rec.yeart-1
rec.size=rec.size[rec.yeart>0]
rec.yeart=rec.yeart[rec.yeart>0]

lst=lst[yeart<16]
lst1=lst1[yeart<16]
yeart=yeart[yeart<16]

mean.lst=mean(lst)

lst=lst-mean(lst)

reclm<-lm(rec.size~factor(rec.yeart)-1)
grlm<-lm(lst1~factor(yeart)*lst-1)

plot(reclm$coef[1:15],grlm$coef[1:15])
cor.test(reclm$coef[1:15],grlm$coef[1:15])

N1=length(lst)
N2=N1+length(rec.size)

#####################################################################################################################################
#fit it with BUGS Cam formulation

n.chains=3;

rho=rnorm(1,0,0.1)
a=rnorm(1,1,0.1)
b=rnorm(1,1,0.1)

X=c(lst,rep(NA,length(rec.size))); Y=c(lst1,rec.size); ngroup=15; group=c(yeart,rec.yeart);

#data=list("N1","N2","X","Y","ngroup","group"); 

data = list(N1=N1,N2=N2,X=X,Y=Y,ngroup=ngroup,group=group)

inits=function() {
	list(mu.A=rnorm(1,0,1),mu.Ar=rnorm(1,0,1),B=rnorm(1,0,1),
	     prec.B=rnorm(1,1,0.1),prec.e=rnorm(1,1,0.1),prec.er=rnorm(1,1,0.1),a=a,b=b,rho=rho,
	     x=rnorm(ngroup,0.0,1),y=rnorm(ngroup,0.0,1),B.g=rnorm(ngroup)
	    )
}

parameters=c("mu.A","mu.Ar","B","sd.B","rho","sd.e","sd.er","sd.A","sd.Ar");

glmm.sim=jags.model(data=data,inits=inits,
file="model grow rec Car.txt",n.chains=n.chains)

samps <- coda.samples(glmm.sim,parameters,n.iter=20000,thin=10)
plot(samps[,1:4])
plot(samps[,5:8])
plot(samps[,9])
autocorr.plot(samps)
gelman.plot(samps)

samps <- coda.samples(glmm.sim,parameters,n.iter=200000,thin=200)

burn.in <- 50000
summary(window(samps, start=burn.in))



