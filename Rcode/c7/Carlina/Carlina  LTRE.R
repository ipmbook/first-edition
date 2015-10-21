## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run LTRE analysis using large parameter set from fitted mixed models
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))


require(mgcv)
require(randomForest)

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

#function adds yearly lambda to param matrix, also removes parameters that
#don't vary with time, and calculates d lambda / d parameters numerically
#returns data.frame with parameters and lambda as columns for analysis
#plus other stuff


LTRE.calcs<-function(params,delta=.Machine$double.eps^(1/3)){

	K.year.i     <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	lambda.i  <- rep(NA,n.years)
	n.params <- dim(params)[1]
	
#Calculate lambda for each yearly kernel
	
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
	mean.params   <- apply(params,1,mean)
	
#Calculate sensitivities to mean kernel for each parameter using finite differences

parm.to.pert         <- 1:n.params
d.lambda.d.theta <- rep(NA,n.params)

for(i in 1:n.params){
	p.pert                 <-  (i==parm.to.pert)*delta
	K.up                   <- mk_K(nBigMatrix,mean.params+p.pert,minsize,maxsize)$K
	lambda.up              <- Re(eigen(K.up,only.values=TRUE)$values[1])
	
	
	K.down                 <-  mk_K(nBigMatrix,mean.params-p.pert,minsize,maxsize)$K
	lambda.down            <-  Re(eigen(K.down,only.values=TRUE)$values[1])
	
	d.lambda.d.theta[i] = 	(lambda.up - lambda.down) / (2*delta)
	}
	
	names(d.lambda.d.theta) <- rownames(params)
	
#remove parameters that don't vary

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

#extract parameters/lambda and sensitivities to formulae look neater

ps.and.l <- LTRE.nums$LTRE.data
d.lambda.d.theta <- LTRE.nums$d.lambda.d.theta

#calculate sensitivity matrix
sens.mat <- outer(d.lambda.d.theta,d.lambda.d.theta,FUN="*")

#calculate sensitivity matrix - cols 1 to 6 contain parameters
var.cov <- cov(ps.and.l[,1:6])

#Calculate 1st order approx
var.terms.1.order <-sens.mat*var.cov

#sum terms
var(ps.and.l[,"lambda.i"])
sum(var.terms.1.order)

#error
(var(ps.and.l[,"lambda.i"])-sum(var.terms.1.order))/var(ps.and.l[,"lambda.i"])*100

#calculate Cont terms and sort
order.terms.1.order <- sort(apply(var.terms.1.order,1,sum)/sum(var.terms.1.order))
order.terms.1.order 

#fit lm model for variation in lambda

fit <- lm(lambda.i~surv.int+surv.z+flow.int+grow.int+grow.z+rcsz.int, data=ps.and.l)

#extract slopes
slopes <- as.numeric(coef(fit)[2:7])

#check slopes ~ sensitivites
plot(d.lambda.d.theta,slopes)
abline(0,1)

#cacluate sensitivity matrix using slopes
sens.mat <- outer(slopes,slopes,FUN="*")

#calculate weighted var-covar matrix
var.terms.lm <- sens.mat*var.cov

#check total variance equal
sum(var.terms.lm) + (summary(fit)$sigma)^2
var(ps.and.l[,"lambda.i"])

#calculate Cont terms and sort
order.terms.lm <- sort(apply(var.terms.lm,1,sum)/sum(var.terms.lm))
order.terms.lm

#do calculation with predict
terms <- predict(fit,type="terms")

cov(terms)

sum(cov(terms)) + (summary(fit)$sigma)^2

#fit gam model
fit <- gam(ps.and.l[,"lambda.i"] ~ s(ps.and.l[,"surv.int"]) + 
                                                            s(ps.and.l[,"surv.z"])    +
                                                            s(ps.and.l[,"flow.int"]) +
                                                            s(ps.and.l[,"grow.int"])+
                                                            s(ps.and.l[,"grow.z"])   +
                                                            s(ps.and.l[,"rcsz.int"]) , data=ps.and.l)

#extract terms using predict
terms <- predict(fit,type="terms")

#weighted var-cov matrix
var.terms.gam <- cov(terms)

rownames(var.terms.gam) <- rownames(var.terms.1.order)
colnames(var.terms.gam)  <- colnames(var.terms.1.order)

var.terms.gam
sum(var.terms.gam)+summary(fit)$scale

var(ps.and.l[,"lambda.i"])

#calculate Cont terms and sort
order.terms.gam <- sort(apply(var.terms.gam,1,sum)/sum(var.terms.gam))
order.terms.gam

###############################################
#Random Forests
###############################################

# 'Tune' the random forest parameter mtry
year.p <- subset(ps.and.l,select=-lambda.i)
lambda.t <- ps.and.l$lambda.i
out<- tuneRF(x=year.p, y=lambda.t,
ntreeTry=500,stepfactor=1,improve=0.02); 
plot(out); ## tells us to use mtry=4 

# Fit the model, computing importance measures
#x=subset(ps.and.l,select=-lambda.i)
lambdaRF <- randomForest(x=year.p,y=lambda.t,ntree=500,importance=TRUE,mtry=4)

graphics.off(); 
dev.new(width=8,height=5); set_graph_pars("panel2");
# verify that the number of trees is enough
# the plot of R^2 vs. number of trees used should saturate 
plot(lambdaRF$rsq,type="l",xlab="Number of trees",ylab="Model r^2"); add_panel_label("a")

# Plot the importance results 
varImpPlot(lambdaRF,type=1,scale=FALSE,main="")
add_panel_label("b")

#sorted Cont terms
order.terms.1.order 
order.terms.lm
order.terms.gam





























