## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Evolutionary calculations for the Carlina stochastic IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
library(lme4)
library(parallel)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
setwd("~/Repos/ipm_book/Rcode/c9")
source("Carlina Demog Funs.R");

source("../utilities/Standard Graphical Pars.R");

# Set simulation parameters

dataf<-read.csv("Carlina survive or flower.fre",header=T)

#fit model ignoring year effects, used in ESS analysis

fit.flow<-glm(flow~lst,family=binomial,data=dataf)

fit.flow.r<-glmer(flow~lst+(1|yeart),data=dataf,family=binomial)

s.dataf<-dataf[dataf$flow==1,]

size.fl<-with(s.dataf,lst[flow==1])

s.dataf$flow<-0
names(s.dataf)[5]<-"die"

dataf<-read.csv("Carlina survive or died.txt",header=T)

dataf<-rbind(dataf,s.dataf)

minsize<-1.4481; 
maxsize<-5.170528;


#iterate model

iterate_model<-function(params,params.sd,n.est,store.env=F) {

	if(store.env==T){
		params.yr<-matrix(NA,ncol=12,nrow=n.est)
	} else {
		params.yr<-NULL
	}
 
		
	nt<-rep(1/nBigMatrix,nBigMatrix)
	Rt<-rep(NA,n.est)
	size.dist<-rep(0,nBigMatrix)
	size.dist.fl<-rep(0,nBigMatrix)
	N.year<-rep(NA,n.est)


	for (year.t in 1:n.est){
		if(year.t%%1000==0) cat("iterate: ", year.t,"\n");

		#params.t<-rnorm(11,params,sd.vec)

		params.t=params + qnorm(runif(12,0.001,0.999))*params.sd

		#sample from multivariate normal distribution for growth intercepts and yearly mean recruit size
		#params.t[c(5,10)]<-mvrnorm(1,mean.grow.rec,VarCovar.grow.rec)

		params.t[c("grow.int","rcsz.int")] <- params[c("grow.int","rcsz.int")]+matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 
		
		year.K<-mk_K(nBigMatrix,params.t,minsize,maxsize)
		
		h <- diff(year.K$meshpts)[1]

		#calculate total seed production

		bt<-sum(year.K$F %*% nt)*h

		params.t["p.r"]<-nrec[sample(1:16,1)]/bt

		
		if (store.env==T){
			params.yr[year.t,] <- params.t
		}

		nt1<-(year.K$P+params.t["p.r"]*year.K$F) %*% nt

				
		sum.nt1<-sum(nt1)*h
		sum.nt<-sum(nt)*h

		Rt[year.t]<-log(sum.nt1/sum.nt)

		N.year[year.t]<-sum.nt

		dist.fl.year<-nt*s_z(year.K$meshpts,params.t)*p_bz(year.K$meshpts,params.t)

		size.dist<-size.dist+nt

		size.dist.fl<-size.dist.fl+dist.fl.year

		nt<-nt1

		
	}

	

	if (store.env==T){
			colnames(params.yr) <- names(params.t)
		}
	
	return(list(size.dist=size.dist,size.dist.fl=size.dist.fl,params.yr=params.yr,
		    N.year=N.year,Rt=Rt,meshpts=year.K$meshpts))
}


##run model and plot some graphs
nBigMatrix=100
n.est<-5000
n.runin<-500

iter<-iterate_model(m.par.true ,m.par.sd.true,n.est,store.env=F)

#calculate stochastic growth rate

mean(iter$Rt[n.runin:n.est],na.rm=T)

plot(iter$N.year,type="b",pch=19)

mean(iter$N.year[n.runin:n.est])

p.size.dist<-iter$size.dist/sum(iter$size.dist)

mean.size<-sum(p.size.dist*exp(iter$meshpts))
mean.size

p.size.dist.fl<-iter$size.dist.fl/sum(iter$size.dist.fl)

mean.size.fl<-sum(p.size.dist.fl*exp(iter$meshpts))
mean.size.fl

mean(exp(size.fl),na.rm=TRUE)

diff<-iter$meshpts[2]-iter$meshpts[1]

quartz()

par(mfrow=c(1,2),bty="l",pty="s")

hist(dataf$lst,breaks=10,main="",col="grey",ylim=c(0,420),xlab="Size")
lines(iter$meshpts,p.size.dist*599/sum(p.size.dist*diff))

#text(locator(1),"a)")

hist(size.fl,breaks=10,main="",col="grey",xlab="Flowering size")
lines(iter$meshpts,p.size.dist.fl*21.4/sum(p.size.dist.fl*diff))

#text(locator(1),"b)")




########################################################################################################################################
#ESS stuff, some functions then code for analysis


invader_gr<-function(params.yr){

	nt<-rep(1/nBigMatrix,nBigMatrix)
	Rt<-rep(NA,n.est)

	for (year.t in 1:n.est){

		params.t<-params.yr[year.t,]

		year.K<-mk_K(nBigMatrix,params.t,minsize,maxsize)

		nt1<-(year.K$P+params.t["p.r"]*year.K$F) %*% nt
		
		sum.nt1<-sum(nt1)
		
		Rt[year.t]<-log(sum.nt1)

		nt<-nt1/sum.nt1
		
		#cat(Rt[year.t],"\n")

	}

return(mean(Rt[n.runin:n.est],na.rm=T))
}

##################################################################################################################################
##ESS analysis
stoch_lambda_beta <- function (beta.0,params){
	params[,"flow.int"] <- params[,"flow.int"]+beta.0
	lambda.s <- invader_gr(params)
	return(lambda.s)
}

find_ESS <- function (params,tol) {
repeat{
#calculate env  for resident - current strategy
	iter<-iterate_model(params ,m.par.sd.true,n.est,store.env=T)

#find best invader
	opt.beta <- optimize(stoch_lambda_beta,lower=-5,upper=5,params=iter$params.yr,maximum=TRUE)

#make best invader the resident
	params["flow.int"] <- params["flow.int"] + opt.beta$maximum

#if best invader had lambda.s=0 to some tolerance we're done
	if(abs(opt.beta$objective)<tol) break()
	
	cat(opt.beta$objective,"   ",opt.beta$maximum,"\n")
	
}
return(params["flow.int"])
}

ESS.iterative.invasion <- find_ESS(m.par.true,0.001)

cat("ESS flowering intercept ",ESS.iterative.invasion,"\n")
ESS.params <- m.par.true
ESS.params["flow.int"] <- ESS.iterative.invasion

iter<-iterate_model(ESS.params ,m.par.sd.true,n.est,store.env=T)
p.size.dist.fl.ESS<-iter$size.dist.fl/sum(iter$size.dist.fl)
p.size.dist.ESS<-iter$size.dist/sum(iter$size.dist)

ESS.size.fl<-sum(p.size.dist.fl.ESS*exp(iter$meshpts))
ESS.size.fl

mean(exp(size.fl),na.rm=TRUE)

diff<-iter$meshpts[2]-iter$meshpts[1]

quartz()
par(mfrow=c(1,2),bty="l",pty="s")
hist(dataf$lst,breaks=7,main="",col="grey",ylim=c(0,420),xlab="Size")
lines(iter$meshpts,p.size.dist*599/sum(p.size.dist*diff))
lines(iter$meshpts,p.size.dist.ESS*599/sum(p.size.dist.ESS*diff),col="red")


#text(locator(1),"a)")

tmp <- hist(size.fl,breaks=seq(2.75,5,length.out=8),main="",col="grey",xlab="Flowering size",ylim=c(0,40))
lines(iter$meshpts,p.size.dist.fl*34.39286/sum(p.size.dist.fl*diff))
lines(iter$meshpts,p.size.dist.fl.ESS*34.39286/sum(p.size.dist.fl.ESS*diff),col="red")

#text(locator(1),"b)")

#savePlot("ESSImages.pdf",type="pdf"); 
#savePlot("ESSImages.bmp",type="bmp"); 
#savePlot("ESSImages.png",type="png"); 


#########################################################################################
#fitness landscapes - assumes ESS environment stored in iter

flow.inter<-seq(-30,-1,1)
flow.inter.diffs <- flow.inter - ESS.params["flow.int"]
n.test<-length(flow.inter)
Ls<-rep(NA,n.test)

# # for (i in 1:n.test){
	
	
	# Ls[i]<-stoch_lambda_beta(flow.inter.diffs[i],iter$params.yr)
	# # p.vec[3]<-flow.inter[i]
	# # i.iter<-iterate.model(p.vec,n.est,store.env=F,DD=T)
	# # p.size.dist.fl.tmp[,i]<-i.iter$size.dist.fl/sum(i.iter$size.dist.fl)
	# # mean.size.fl.I[i]<-sum(p.size.dist.fl.tmp[,i]*exp(i.iter$meshpts))
	# #var.size.fl.I[i]<-sum(p.size.dist.fl.tmp[,i]*exp(2*i.iter$meshpts))-mean.size.fl.I[i]^2

	# cat(i,"  ",flow.inter[i],"  ",exp(Ls[i]),"\n")
# }

#Ls <- simplify2array(mclapply(1:n.test,function(i) stoch_lambda_beta(flow.inter.diffs[i],iter$params.yr),mc.cores=4))

Ls <- simplify2array(mclapply(flow.inter.diffs, stoch_lambda_beta,params=iter$params.yr,mc.cores=24))

# times.mc <- rep(NA,24)

# for(i in 1:12){
	
	# tmp <- system.time(Ls <- simplify2array(mclapply(1:n.test,function(i) stoch_lambda_beta(flow.inter.diffs[i],iter$params.yr),mc.cores=2*i)))
	
	# times.mc[i] <- tmp["elapsed"]
	# cat(i,"\n")
	
# }

# plot(1:12,times.mc[1:12],type="b",pch=19)


# system.time(Ls <- simplify2array(mclapply(1:n.test,function(i) stoch_lambda_beta(flow.inter.diffs[i],iter$params.yr),mc.cores=24)))
# system.time(Ls <- simplify2array(mclapply(flow.inter.diffs, stoch_lambda_beta,params=iter$params.yr,mc.cores=12)))
# system.time(for (i in 1:n.test){ Ls[i]<-stoch_lambda_beta(flow.inter.diffs[i],iter$params.yr)})


quartz()

par(mfrow=c(1,1),bty="l",pty="s",pch=19)

plot(flow.inter,Ls,type="n",xlab=expression("Intercept of flowering function " * italic(beta * scriptstyle(0))),
	ylab=expression("Fitness log("*lambda[s]*")"),ylim=c(-0.8,0.05))


SE <- summary(fit.flow.r)$coef[1,2]
inter <- summary(fit.flow.r)$coef[1,1]
x.corners=as.numeric(c(inter-2*SE,inter-2*SE,inter+2*SE,inter+2*SE))

polygon(x.corners,c(-0.8,0,0,-0.8),col="grey90",border=0)
lines(c(inter,inter),c(-0.8,0),col="turquoise",lwd=2)

points(flow.inter,Ls,type="b")

abline(h=max(Ls))
#abline(h=max(Ls)-0.01,col="red")

dev.copy2eps(file="~/Repos/ipm_book/c9/figures/CarlinaFitnessLandscape.eps")

#text(locator(1),"a)")
#Plot pips

nBigMatrix=100
n.est<-5000
n.runin<-300


n.betas<-200
beta.s <- seq(-30,-1,length=n.betas)
params <- m.par.true

lambdas.pip <- matrix(NA,ncol=n.betas,nrow=n.betas)

params.res <- m.par.true

#Constructing the Pip takes a while as we're calculating 40,000 stochastic growth rates
#each using a 100x100 matrix and 5000 iterates, let's use one we made earlier...

# # for(i in 1:n.betas){
	
	# params.res["flow.int"] <- beta.s[i]
	
	# iter.R<-iterate_model(params.res ,m.par.sd.true, n.est, store.env=T)
	
	# beta.diffs <- beta.s - beta.s[i]

    # lambdas.pip[i,] <- simplify2array(mclapply(beta.diffs,stoch_lambda_beta,params=iter.R$params.yr,mc.cores=12))
    
    # #lambdas.pip[i,] <- sapply(beta.diffs,stoch_lambda_beta,params=iter.R$params.yr)

# cat(i,"\n")

# }

#save(beta.s,lambdas.pip,file="Pip data.Rdata")
load("Pip data.Rdata")

lambdas.pip <- lambdas.pip[beta.s < -5,beta.s < -5]

beta.s <- beta.s[beta.s < -5]

set_graph_pars("panel2"); 
image(beta.s,beta.s,lambdas.pip<0, xlab=expression("Resident strategy " * italic(beta) [0]),
ylab=expression("Invader strategy " * italic(beta) [0]))

abline(v=ESS.iterative.invasion)

add_panel_label("a")

mip <- lambdas.pip * t(lambdas.pip)

diag(mip) <- -1

image(beta.s,beta.s,mip<0, xlab=expression("Resident strategy " * italic(beta) [0]),
ylab=expression("Invader strategy " * italic(beta) [0]))

#abline(v=ESS.iterative.invasion)
abline(0,1)

add_panel_label("b")

#dev.copy2eps(file="~/Repos/ipm_book/c9/figures/CarlinaPip.eps")

#library(mgcv)

# Ls <- as.vector(lambdas.pip)

# n.betas <- length(beta.s)

# beta.R <- rep(beta.s,rep(n.betas,n.betas))

# beta.I <- rep(beta.s,n.betas)

# fit.surf <- gam(Ls ~ s(beta.R) + s(beta.I) + s(beta.I,beta.R) )
# summary(fit.surf) 

# fitted.fits <- matrix(predict(fit.surf),nrow=n.betas)


n.betas <- length(beta.s)
smoothed.fit <- matrix(apply(lambdas.pip,1,smooth),nrow=n.betas, byrow=TRUE)
smoothed.fit <- matrix(apply(smoothed.fit,2,smooth),nrow=n.betas)


set_graph_pars("panel2"); 
image(beta.s,beta.s,smoothed.fit<0, xlab=expression("Resident strategy " * italic(beta) [0]),
ylab=expression("Invader strategy " * italic(beta) [0]))

abline(v=ESS.iterative.invasion,lwd=1)
text(-23,-15,expression(italic(lambda[S])*" > 1"),cex=1.1)
text(-20,-25,expression(italic(lambda[S])*" < 1"),cex=1.1)
text(-10,-15,expression(italic(lambda[S])*" > 1"),cex=1.1)
text(-11,-7,expression(italic(lambda[S])*" < 1"),cex=1.1)


add_panel_label("a")
mip <- smoothed.fit * t(smoothed.fit)
diag(mip) <- -1
image(beta.s,beta.s,mip<0, xlab=expression("Flowering strategy " * italic(beta) [0]),
ylab=expression("Flowering strategy " * italic(beta) [0]))

#abline(v=ESS.iterative.invasion)
abline(0,1)
add_panel_label("b")

dev.copy2eps(file="~/Repos/ipm_book/c9/figures/CarlinaPip.eps")





















