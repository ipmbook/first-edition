## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Simulate the Carlina stochastic density independent IBM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
library(lme4)
library(parallel)
set.seed(523886743)

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c2",sep="")); 

source("../utilities/Standard Graphical Pars.R");

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c7/Carlina",sep="")); 

source("Carlina Demog Funs DI.R") 

#######################################################################################################
#  Simulate DI Carlina IBM
#######################################################################################################
# Set simulation parameters
init.pop.size <- 100
init.mean.z <- 3
init.sd.z <- 0.5
n.yrs <- 50
n.reps <- 30

# pop.sizes.sim <-matrix(NA, nrow=n.yrs,ncol=n.reps)

# for(reps in 1:n.reps){
 # source("Carlina DI IBM.R") 
 # pop.sizes.sim[,reps] <- pop.size.t
 # store.min.z.t        <- if (reps==1) min.z.t else c(store.min.z.t,min.z.t)  
 # store.max.z.t        <- if (reps==1) max.z.t else c(store.max.z.t,max.z.t)
 # store.mean.z.t       <- if (reps==1) mean.z.t else c(store.mean.z.t,mean.z.t)  
 # store.mean.fl.z.t    <- if (reps==1) mean.fl.z.t else c(store.mean.fl.z.t,mean.fl.z.t)  
# }

# save(sim.data,pop.sizes.sim,store.min.z.t,store.max.z.t,store.mean.z.t,store.mean.fl.z.t,file="CarlinaIBMsim.Rdata")

load(file="CarlinaIBMsim.Rdata")

#dev.copy2eps(file="~/Repos/ipm_book/c7/figures/CarlinaStocDIIBM.eps")

#####################################################################
#Run the IPM with the true parameters estimates
#####################################################################

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

	for (year.t in 1:n.est){
		if(year.t%%1000==0) cat("iterate: ", year.t,"\n");

		#params.t<-rnorm(11,params,sd.vec)

		params.t=params + qnorm(runif(12,0.001,0.999))*params.sd

		#sample from multivariate normal distribution for growth intercepts and yearly mean recruit size
		#params.t[c(5,10)]<-mvrnorm(1,mean.grow.rec,VarCovar.grow.rec)

		params.t[c("grow.int","rcsz.int")] <- params[c("grow.int","rcsz.int")]+matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 
		
		year.K<-mk_K(nBigMatrix,params.t,minsize,maxsize)
	
		
		if (store.env==T){
			params.yr[year.t,] <- params.t
		}

		nt1<-year.K$K %*% nt
	
		sum.nt1<-sum(nt1)
		
		Rt[year.t]<-log(sum.nt1)

		dist.fl.year <- nt * s_z(year.K$meshpts,params.t) * p_bz(year.K$meshpts,params.t)

		size.dist <- size.dist + nt

		size.dist.fl<-size.dist.fl+dist.fl.year

		nt <- nt1 / sum(nt1)
		
	}

    if (store.env==T){
			colnames(params.yr) <- names(params.t)
		}
	
	size.dist <- size.dist / sum(size.dist)
	size.dist.fl <- size.dist.fl / sum(size.dist.fl)
	
	return(list(size.dist=size.dist,size.dist.fl=size.dist.fl,params.yr=params.yr,
		    Rt=Rt,meshpts=year.K$meshpts))
}


##run model and plot some graphs


nBigMatrix <- 100
n.est <- 10500
n.runin <- 500
minsize <- 1.5
maxsize <- 5

iter <- iterate_model(m.par.true ,m.par.sd.true, n.est, store.env=F)

dev.new(height=4,width=8);
set_graph_pars("panel2"); 
matplot(pop.sizes.sim,type="l",log="y",xlab="Time (years)",ylab="Population size",yaxt="n")
axis(2,at=c(1,10,1000,100000),labels=c("1","10","1000","100000"))
abline(h=2000000,lty=3,col="red")
add_panel_label("a")
pop.sizes.sim[pop.sizes.sim==0] <- NA
lambda.t <- as.vector(log(pop.sizes.sim[11:n.yrs,]/pop.sizes.sim[10:(n.yrs-1),]))
hist(lambda.t,main="",xlab = expression(r[t] *"=log(N(t+1)/N(t))"),col="grey",freq=FALSE,ylim=c(0,0.55),breaks=16)
add_panel_label("b")

#calculate stochastic growth rate

stoch.lambda <- mean(iter$Rt[n.runin:n.est],na.rm=T)
SE.stoch.lambda <- sd(iter$Rt[n.runin:n.est],na.rm=T) / sqrt(length(iter$Rt[n.runin:n.est]))
stoch.lambda.IBM <- mean(lambda.t,na.rm=T)
SE.lambda.IBM <- sd(lambda.t,na.rm=TRUE)/sqrt(sum(!is.na(lambda.t)))
cat("Log stochastic growth rate IPM =",stoch.lambda," from IBM =",stoch.lambda.IBM,"\n")
cat("Approximate 95% CI for log stochastic growth rate from IBM ",stoch.lambda.IBM - 2*SE.lambda.IBM, " to ", stoch.lambda.IBM + 2*SE.lambda.IBM,"\n")
cat("Approximate 95% CI for log stochastic growth rate from IPM ",stoch.lambda - 2*SE.stoch.lambda, " to ", stoch.lambda + 2*SE.stoch.lambda,"\n")
lines(c(stoch.lambda,stoch.lambda),c(0,0.55),col="turquoise",lwd=2)

points(density(iter$Rt[n.runin:n.est]),col="red",type="l")
#dev.copy2eps(file="~/Repos/ipm_book/c7/figures/CarlinaStocDIIBM.eps")

#Mean plant size from IBM (log scale)
mean(store.mean.z.t )
#from IPM (log scale)
sum(iter$size.dist * iter$meshpts)

#Mean flowering size from IBM
mean(store.mean.fl.z.t,na.rm=TRUE)
#from IPM
sum(iter$size.dist.fl * exp(iter$meshpts))

#########################################################
#Lambda S from finite runs
#Generating sets of 20 years parameters and then calculating lambda S


iterate_model.FR<-function(params,params.sd,n.years,n.est) {

#Construct the yearly kernels

	K.year.i <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		
		params.t=params + qnorm(runif(12,0.001,0.999))*params.sd

		params.t[c("grow.int","rcsz.int")] <- params[c("grow.int","rcsz.int")]+matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 

		year.K<-mk_K(nBigMatrix,params.t,minsize,maxsize)
		K.year.i[i,,] <- year.K$K
	}
 
#initialize variables	

	nt<-rep(1/nBigMatrix,nBigMatrix)
	Rt<-rep(NA,n.est)

#Iterate model

	for (year.t in 1:n.est){
		if(year.t%%5000==0) cat("iterate: ", year.t,"\n");

		#Select year at random
		
		year.i <- sample(1:n.years,1)
		
		#iterate model with year-specific kernel
		nt1<-K.year.i[year.i,,] %*% nt
	
		sum.nt1<-sum(nt1)
		
		#Calculate log growth rate we normalise so sum.nt=1 
		Rt[year.t]<-log(sum.nt1)

		nt <- nt1 /sum.nt1
		
		#cat(Rt[year.t],"  ",sum.nt,"\n")
		
				
	}
	
	return(mean(Rt[n.runin:n.est],na.rm=T))
}


nBigMatrix <- 100
n.est <- 10500
n.runin <- 500
minsize <- 1.5
maxsize <- 5


lambdaS <- rep(NA,500)

# # for(i in 1:500){
	# lambdaS[i] <- iterate_model.FR(m.par.true,m.par.sd.true,20,n.est)
	# cat(i,"\n")
# }

lambdaS <- simplify2array(mclapply(1:500, function(i) iterate_model.FR(m.par.true,m.par.sd.true,20,n.est),mc.cores=6))
dev.new(height=6,width=8);
set_graph_pars("panel1"); 
hist(lambdaS,breaks=15,freq=FALSE,col="grey",main="",xlab=expression("log("*lambda[S]*")"))
lines(c(stoch.lambda,stoch.lambda),c(0,2.5),col="turquoise",lwd=3)
sum(lambdaS<0)/length(lambdaS)

#dev.copy2eps(file="~/Repos/ipm_book/c7/figures/LambdaSfinite.eps")


#Simple minded approach - split Rt time series into 20 years segments and calculate log(lambdaS) for each

no.runs <- length(iter$Rt[(n.runin+1):n.est])/20

lambda.test <- data.frame(Rt=iter$Rt[(n.runin+1):n.est], run = gl(no.runs,20))

lambdas <- summaryBy(Rt ~ run, data=lambda.test)

with(lambdas,hist(Rt.mean,breaks=15,freq=FALSE,col="grey",main="",xlab=expression("log("*lambda[t]*")")))

#looks remarkably similar range -0.6 to 0.6, mean value about 0.12

#dev.copy2eps(file="~/Repos/ipm_book/c7/figures/LambdaSfinite.eps")







