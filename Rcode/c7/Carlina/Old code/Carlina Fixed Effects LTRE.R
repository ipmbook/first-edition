## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run LTRE analysis using large parameter set from fitted mixed models
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
library(lme4)
library(MCMCglmm)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c2",sep="")); 

source("../utilities/Standard Graphical Pars.R");

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c7/Carlina",sep="")); 


source("Carlina Demog Funs DI.R") 

#####################################################################
#LTRE analysis
#####################################################################

nBigMatrix <- 100
n.est <- 20000
n.runin <- 500
minsize <- 1.5
maxsize <- 5
n.years <-20

#function adds yearly lambda to param matrix, also removed parameters that
#don't vary with time, and calculates d lambda / d parameters numerically
#returns data.frame with parameters and lambda as columns for analysis
#plus other stuff


LTRE.calcs<-function(params,delta=.Machine$double.eps^(1/3)){

	K.year.i     <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	lambda.i  <- rep(NA,n.years)
	n.params <- dim(params)[1]
	
	for(i in 1:n.years){
		year.K         <- mk_K(nBigMatrix,params[,i],minsize,maxsize)
		lambda.i[i] <- Re(eigen(year.K$K,only.values=TRUE)$values[1])
		K.year.i[i,,] <- year.K$K
		#cat(lambda.i[i],"\n")
	}

	h <- year.K$h; 
	meshpts <- year.K$meshpts
	
#Calculate mean kernel, and mean of each params

	mean.kernel   <- apply(K.year.i,2:3,mean)
	mean.params <- apply(params,1,mean)
	
#Calculate sensitivities to mean kernel for each parameter

parm.to.pert         <- 1:n.params
d.lambda.d.theta <- rep(NA,n.params)

for(i in 1:n.params){
	p.pert                 <-  (i==parm.to.pert)*delta
	K.up                    <- mk_K(nBigMatrix,mean.params+p.pert,minsize,maxsize)$K
	lambda.up         <- Re(eigen(K.up,only.values=TRUE)$values[1])
	
	
	K.down               <-mk_K(nBigMatrix,mean.params-p.pert,minsize,maxsize)$K
	lambda.down    <- Re(eigen(K.down,only.values=TRUE)$values[1])
	
	d.lambda.d.theta[i] = 	(lambda.up - lambda.down) / (2*delta)
	}
	
	names(d.lambda.d.theta) <- rownames(params)
	
	which.vary<- which(apply(params,1,sd)>0);	
	LTRE.data <- rbind(params[which.vary,],lambda.i)
	LTRE.data.df<- data.frame(t(LTRE.data)); 
	d.lambda.d.theta <- d.lambda.d.theta[which.vary]
	 return(list(meshpts=meshpts, h=h,mean.kernel=mean.kernel,LTRE.data=LTRE.data.df,
	 d.lambda.d.theta=d.lambda.d.theta))

}

load(file="Big param matrix.Rdata")

n.years <-5000
LTRE.nums <- LTRE.calcs(store.params)

ps.and.l <- LTRE.nums$LTRE.data
d.lambda.d.theta <- LTRE.nums$d.lambda.d.theta

sens.mat <- outer(d.lambda.d.theta,d.lambda.d.theta,FUN="*")

var.cov <- cov(ps.and.l[,1:6])

var.terms.1.order <-sens.mat*var.cov
var(ps.and.l[,"lambda.i"])
sum(var.terms.1.order)

fit <- lm(lambda.i~surv.int+surv.z+flow.int+grow.int+grow.z+rcsz.int, data=ps.and.l)

slopes <- as.numeric(coef(fit)[2:7])

plot(d.lambda.d.theta,slopes)
abline(0,1)

sens.mat <- outer(slopes,slopes,FUN="*")
var.terms.lm <- sens.mat*var.cov
sum(var.terms.lm) + (summary(fit)$sigma)^2

var(ps.and.l[,"lambda.i"])


terms <- predict(fit,type="terms")

cov(terms)

sum(cov(terms)) + (summary(fit)$sigma)^2

fit <- gam(ps.and.l[,"lambda.i"] ~ s(ps.and.l[,"surv.int"]) + 
                                                            s(ps.and.l[,"surv.z"])    +
                                                            s(ps.and.l[,"flow.int"]) +
                                                            s(ps.and.l[,"grow.int"])+
                                                            s(ps.and.l[,"grow.z"])   +
                                                            s(ps.and.l[,"rcsz.int"]) , data=ps.and.l)

terms <- predict(fit,type="terms")

var.terms.gam <- cov(terms)

rownames(var.terms.gam) <- rownames(var.terms.1.order)
colnames(var.terms.gam)  <- colnames(var.terms.1.order)

var.terms.gam
sum(var.terms.gam)+summary(fit)$scale

var(ps.and.l[,"lambda.i"])




#save(params.plus.lambda,file="Big params plus lambda.Rdata")



