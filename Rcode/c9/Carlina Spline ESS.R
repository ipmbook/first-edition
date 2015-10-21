###################################################################################################
#  This program seeks an ESS function-valued flowering strategy in the stochastic monocarp
#  monocarp model with Carlina parameters and density-dependent recruitment, by finding a resident
#  for which the gradient of the invader fitness is zero (1st order condition for ESS).    
###################################################################################################

rm(list=ls(all=TRUE))
## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c9",sep="")); 

source("../utilities/Standard Graphical Pars.R");
source("Monocarp Demog SplineFuns.R");  
source("Carlina Demog Funs.R"); # Carlina functions and parameters 
require(fda); 

###################################################################################################
#  Function to iterate monocarp stochastic IPM with size variation and a function-valued flowering
#  strategy specified by a B-spline basis object B and spline coefficients cj. Argument params.yr  
#  is a sequence of year-specific parameters except for recruitment probability. It is returned
#  unchanged except that p.r is set to the sequence of actual recruitment probabilities resulting
#  from the resident flowering strategy.     
###################################################################################################
iterate_spline_resident <- function(params.yr,B,cj,m,L,U) {
    h <- (U-L)/m; meshpts <- L+ ((1:m)-0.5)*h
    nt <- matrix(1,m,1); 
    n.iter <- nrow(params.yr); 
    Nt <- numeric(n.iter); Nt[1] <- sum(h*nt); 
	
 	for(gen in 2:n.iter){	
        m.par.year <- params.yr[gen-1,]
        m.par.year["p.r"] <- 1; # So that F%*%n = total seeds, not seeds*p.r 
              
        # P <- h * (outer(meshpts, meshpts, P_z1z_spline, m.par = m.par.year, B = B, cj = cj))
        # Equivalent, but much faster (m spline evaluations instead of m^2): 
        G <- outer(meshpts, meshpts, G_z1z, m.par = m.par.year);
        d <- (1 - p_bz_spline(meshpts, B,cj)) * s_z(meshpts, m.par.year); 
        P <- h*(G%*%diag(as.vector(d)))
    
    
		F <- h* F_z1z_spline_vec(meshpts, meshpts, m.par=m.par.year, B=B, cj=cj)
        
        seeds=F%*%nt; 
        total.seeds=sum(h*seeds); 
        p.r.yr <- Recr.year[gen]/total.seeds; 
        
        nt <- seeds*p.r.yr + P%*%nt; 
        params.yr[gen-1,"p.r"] <- p.r.yr;
        Nt[gen] <- sum(h*nt); 
 	
        if(gen%%100==0) cat(gen,p.r.yr,sum(h*nt),"\n")
}

return(list(nt=nt,Nt=Nt,params.yr=params.yr))
}
###################################################################################################
#  Function to compute long-term growth rate of an invader with spline coefficient cj, 
#  experiencing year-specific parameters and recruitment probabilities in params.yr  
#  It returns final population structure, year-specific growth rates, and estimated log(lambda_{S,I}).   
###################################################################################################
spline_invader_gr <- function(params.yr,B,cj,m,L,U) {
    n.iter=nrow(params.yr); 
    h <- (U-L)/m; meshpts <- L+ ((1:m)-0.5)*h
    nt <- matrix(1/m,m,1); Rt <- numeric(n.iter);  
 
	for(gen in 2:n.iter){
        m.par.year <- params.yr[gen-1,]; 
        # P <- h * (outer(meshpts, meshpts, P_z1z_spline, m.par = m.par.year, B = B, cj = cj))
        # Equivalent, but much faster (m spline evaluations instead of m^2): 
        G <- outer(meshpts, meshpts, G_z1z, m.par = m.par.year);
        d <- (1 - p_bz_spline(meshpts, B,cj)) * s_z(meshpts, m.par.year); 
        P <- h*(G%*%diag(as.vector(d)))
        
		F <- h* F_z1z_spline_vec(meshpts, meshpts, m.par=m.par.year, B=B, cj=cj)
		
        nt <- (P+F)%*%nt; 
        Rt[gen] <- sum(nt); nt<-nt/Rt[gen]; 
        if(gen%%300==0) cat(gen,Rt[gen],"\n")
}
		lambda.hat = mean(log(Rt[101:n.iter]))
return(list(nt=nt,Rt=Rt,lambda.hat=lambda.hat))
}

#################################################################### 
# Objective function: penalized lambda of invader 
####################################################################
lambdaS.inv=function(cj,params.yr,B,m,L,U,max.slope,X) {
    slope <- diff(X%*%cj)/dz; 
    mslope <- max(slope); 
    pen=ifelse(mslope > max.slope, 1000*((mslope-max.slope)^2), 0) 
    lam<-spline_invader_gr(params.yr,B,cj,m,L,U)$lambda.hat
    
    return(lam-pen)
}    

#################################################################### 
# 1 + rms(gradient) of objective function 
####################################################################
gradnorm.inv <- function(cj,params.yr,B,m,L,U,max.slope,X) {
	    res<- iterate_spline_resident(params.yr,B,cj,m,L,U)
		cat("resident","\n"); 
		params.yr=as.matrix(res$params.yr)
		npar=length(cj); mat=diag(1+npar); out=numeric(1+npar); 
		for(i in 1:(1+npar)){
			out[i]=lambdaS.inv(cj+0.01*mat[1:npar,i],params.yr,B,m,L,U,max.slope,X)
			cat("gradient", i, "\n")
		}
        grad=100*(out[1:npar]-out[1+npar])
		return(1+sqrt(sum(grad^2))); 
}        

################################################################################
## Generate yearly parameter values ONCE AND FOR ALL
## This is crucial, so that strategy fitness is only a function of parameters
## and two identical calls to the objective function give the same value
###############################################################################
n.iter=600;   
params.yr <- matrix(NA,ncol=length(m.par.true),nrow=n.iter)
params.yr <- data.frame(params.yr); 
names(params.yr) <- names(m.par.true); 
for(gen in 1:n.iter){#Iterate the model
    m.par.year <- m.par.true + qnorm(runif(12,0.001,0.999))*m.par.sd.true
    m.par.year[c("grow.int","rcsz.int")] <- m.par.true[c("grow.int","rcsz.int")] + 
           matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 
    params.yr[gen,] <- m.par.year; 
}
params.yr[,"p.r"] <- 1; # So that F%*%n = total seeds, not seeds*p.r
params.yr <- as.matrix(params.yr); 

## Generate sequence of recruit numbers  
nrec<-c(20,42,12,17,8,19,58,45,44,2,56,25,75,92,94,6,4,34,104); # Carlina values
Recr.year <- sample(nrec,5000,replace=TRUE) # more than we'll need


#################################################################################
#  Set things up for optimizing the strategy  
#################################################################################
## create spline basis for flowering strategy 
L <- (-1); U <- 8; ## WIDE size range
B <- create.bspline.basis(rangeval=c(L,U),nbasis=6, norder=4)
zvals <- seq(L,U,length=250); 
X <- eval.basis(zvals,B); dz=(U-L)/250; # for computing slope 

## initialize spline coefficients to approximate estimated true strategy (not an ESS!) 
cj <- runif(ncol(X)); 
yvals <- m.par.true["flow.int"] + zvals*m.par.true["flow.z"]; 
mse <- function(cj) { sqrt(mean((yvals - X%*%cj)^2)) }
out <- optim(par=cj,f=mse,control=list(trace=0,maxit=2500)); 
for(j in 1:10) {out <- optim(par=out$par,f=mse,control=list(trace=0,maxit=2500)); cat(out$value,"\n") } 
cj <- out$par;

## Test: does the model run? 
res <- iterate_spline_resident(params.yr,B,cj,m=90,L,U);   
lambdaS.inv(cj,params.yr=res$params.yr,B,m=60,L,U,max.slope=1.05*m.par.true["flow.z"],X); 
gradnorm.inv(cj,params.yr,B,m=60,L,U,max.slope=1.05*m.par.true["flow.z"],X); 

#################################################################### 
# Minimize norm of invader gradient 
####################################################################
cjStep=as.list(1:10); cjStep[[1]]=cj; 
for(j in 2:10) {
    out=optim(par=cjStep[[j-1]],f=gradnorm.inv,control=list(maxit=100,trace=4),
        params.yr=params.yr,B=B,m=90,L=L,U=U,max.slope=1.05*m.par.true["flow.z"],X=X)
    cjStep[[j]]=out$par
    u <- X%*%cjStep[[1]]; p <- exp(u)/(1+exp(u));
    for(k in 2:j) {u<- X%*%cjStep[[k]]; p<- cbind(p,exp(u)/(1+exp(u)));}   
    matplot(zvals,p,col=c("black",rep("grey40",20)),type="l",lwd=c(2,rep(1,20)), xlab="Size z",
    lty=1,ylab="Flowering probability",xlim=c(1.5,6)); 
    points(zvals,p[,k],col="red",type="l",lwd=2); 
    u1 <- -14.3 + zvals*m.par.true["flow.z"]; p1=exp(u1)/(1+exp(u1))
    points(zvals,p1,type="p",col="blue"); 
}

save.image("Carlina Spline ESS2.Rdata")  
#################################################################### 
# Plot the results 
####################################################################
graphics.off(); 
dev.new(height=5.5,width=9,mar=c(4,4,1,1)); set_graph_pars("panel1"); 
par(cex.axis=1.2,cex.lab=1.4,mgp=c(2.2,1,0),yaxs="i"); 
cjt=cjStep[[1]]; u <- X%*%cjt; p <- exp(u)/(1+exp(u));
for(k in 2:10) {
        cjt=cjStep[[k]]; u <- X%*%cjt; p <- cbind(p,exp(u)/(1+exp(u)));
}    
matplot(zvals,p[,c(1,10)],col="black",type="l",lwd=2,lty=c(2,1), xlab="Size z",
ylab="Flowering probability",xlim=c(1.5,6),ylim=c(0,1)); 

u1 <- -14.3388 + zvals*m.par.true["flow.z"]; p1=exp(u1)/(1+exp(u1))
points(zvals,p1,type="p",col="black",cex=0.8,lwd=2,pch=19); 

legend("topleft",legend=c("Initial strategy","Function-valued ESS","Logistic ESS"),
lwd=2,lty=c(2,1,NA),pch=c(NA,NA,16),bty="n",inset=0,cex=1.2);  

dev.copy2eps(file="../../c9/figures/CarlinaSplineESS.eps"); 
