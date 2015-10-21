###################################################################################################
#  This program simulates a sequence of invasions by alternative function-valued flowering
#  strategies, in the stochastic monocarp model with Carlina parameters and density-dependent 
#  recruitment.    
###################################################################################################

rm(list=ls(all=TRUE))
## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 

source("../utilities/Standard Graphical Pars.R");
source("Monocarp Demog SplineFuns.R");  
source("Carlina Demog Funs.R"); # Carlina functions and parameters 
require(fda); 



###################################################################################################
#  Function to iterate monocarp stochastic IPM with size variation and a function-valued flowering
#  strategy specified by a B-spline basis object B and spline coefficients cj 
#  It stores and returns the sequence of year-specific parameters and population sizes    
###################################################################################################
iterate_spline_resident <- function(params,n.iter,B,cj,m,L,U) {
    
    nrec<-c(20,42,12,17,8,19,58,45,44,2,56,25,75,92,94,6,4,34,104);
    
    h <- (U-L)/m; meshpts <- L+ ((1:m)-0.5)*h
    nt <- matrix(1,m,1); 
    Nt <- numeric(n.iter); Nt[1] <- sum(h*nt); 
    params.yr <- matrix(NA,ncol=length(params),nrow=n.iter)
    params.yr <- data.frame(params.yr); 
    names(params.yr) <- names(params); 
    
	for(gen in 2:n.iter){	#Iterate the model
	    #generate yearly parameters
        m.par.year <- params + qnorm(runif(12,0.001,0.999))*m.par.sd.true
        m.par.year[c("grow.int","rcsz.int")] <- params[c("grow.int","rcsz.int")] + 
                matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 
        
        m.par.year["p.r"] <- 1; # So that F%*%n = total seeds, not seeds*p.r 
        params.yr[gen-1,] <- m.par.year; 
        
        Recr.year <- sample(nrec,1)
        P <- h * (outer(meshpts, meshpts, P_z1z_spline, m.par = m.par.year, B = B, cj = cj))
		F <- h* F_z1z_spline_vec(meshpts, meshpts, m.par=m.par.year, B=B, cj=cj)
        
        seeds=F%*%nt; 
        total.seeds=sum(h*seeds); 
        p.r.yr <- Recr.year/total.seeds; 
        
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
#  It returns final population structure, year-specific growth rates, and estimated log(lambda_I).   
###################################################################################################
spline_invader_gr <- function(params.yr,B,cj,m,L,U) {
    n.iter=nrow(params.yr); 
    h <- (U-L)/m; meshpts <- L+ ((1:m)-0.5)*h
    nt <- matrix(1/m,m,1); Rt <- numeric(n.iter);  
 
	for(gen in 2:n.iter){
        m.par.year <- params.yr[gen-1,]; 
        P <- h * (outer(meshpts, meshpts, P_z1z_spline, m.par = m.par.year, B = B, cj = cj))
		F <- h* F_z1z_spline_vec(meshpts, meshpts, m.par=m.par.year, B=B, cj=cj)
		
        nt <- (P+F)%*%nt; 
        Rt[gen] <- sum(nt); nt<-nt/Rt[gen]; 
        if(gen%%100==0) cat(gen,Rt[gen],"\n")
}
		lambda.hat = mean(log(Rt[101:n.iter]))
return(list(nt=nt,Rt=Rt,lambda.hat=lambda.hat))
}

#################################################################### 
# Objective function for maximizing invader lambda 
####################################################################
objfun.inv=function(cj,params.yr,B,m,L,U,max.slope,X) {
    slope <- diff(X%*%cj)/dz; 
    mslope <- max(slope); 
    pen=ifelse(mslope > max.slope, 1000*((mslope-max.slope)^2), 0) 
    lam<-spline_invader_gr(params.yr=params.yr,B=B,cj=cj,m=m,L=L,U=U)$lambda.hat
    return(-lam+pen)
}    

d.objfun.inv <- function(cj,params.yr,B,m,L,U,max.slope,X) {
		npar=length(cj); 
	    mat=diag(1+npar); 
        gfun=function(i) {
            return(objfun.inv(cj+0.01*mat[1:npar,i],params.yr,B,m,L,U,max.slope,X))
        } 
        out=numeric(npar+1); for(i in 1:(npar+1)) {out[i]=gfun(i); cat("gradient",i,"\n"); } 
        cat(100*(out[1:npar]-out[1+npar]),"\n"); 
        return(100*(out[1:npar]-out[1+npar])); 
}  

#################################################################### 
# Simulate a sequence of invasions 
####################################################################

## create spline basis 
L <- (-2); U <- 8; ## WIDE size range
B <- create.bspline.basis(rangeval=c(L,U),nbasis=8, norder=4)
X <- eval.basis(seq(L,U,length=250),B); dz=(U-L)/250; # for computing slope 

## initialize parameters to approximate estimated true flowering strategy (not an ESS!) 
cj <- runif(ncol(X)); 
zvals <- seq(L,U,length=250); 
yvals <- m.par.true["flow.int"] + zvals*m.par.true["flow.z"]; 
mse <- function(cj) { sqrt(mean((yvals - X%*%cj)^2)) }
out <- optim(par=cj,f=mse,control=list(trace=0,maxit=2500)); 
for(j in 1:20) { out <- optim(par=out$par,f=mse,control=list(trace=0,maxit=2500)); cat(out$value,"\n") } 
cj <- out$par;

## Test: does the model run? 
res <- iterate_spline_resident(params=m.par.true,n.iter=600,B=B,cj=cj,m=60,L=L,U=U);   
objfun.inv(cj,params.yr=as.matrix(res$params.yr),B,m=60,L,U,max.slope=1.05*m.par.true["flow.z"],X=X); 
d.objfun.inv(cj,params.yr=as.matrix(res$params.yr),B,m=60,L,U,max.slope=1.05*m.par.true["flow.z"],X) 

### simulate sequential invasion
cjStep=as.list(1:10); cjStep[[1]]=cj;  
for(j in 2:10) { 
    res <- iterate_spline_resident(params=m.par.true,n.iter=600,B=B,cj=cjStep[[j-1]],m=60,L=L,U=U);  
	
    out=optim(par=cjStep[[j-1]],method="BFGS",f=objfun.inv,gr=d.objfun.inv,control=list(maxit=40+5*j,trace=4,REPORT=1),
        params.yr=as.matrix(res$params.yr),B=B,m=60,L=L,U=U,max.slope=1.05*m.par.true["flow.z"],X=X)
    
	cjStep[[j]]=out$par;
	save.image(file="CarlinaFuncValESSsim=BFGS.Rdata"); 
    cat("Invasion step ", j, "completed")
    cj=cjStep[[1]]; u <- X%*%cj; p <- exp(u)/(1+exp(u));
    plot(zvals,p,type="l",lty=1,col=1); 
    for(k in 2:j) {
        cj=cjStep[[k]]; u <- X%*%cj; p <- exp(u)/(1+exp(u));
        points(zvals,p,type="l",lty=1,col=k); 
    }    
}    

#################################################################### 
# Plot the results 
####################################################################
graphics.off(); 
dev.new(); set_graph_pars("panel1"); 
cj=cjStep[[1]]; u <- X%*%cj; p <- exp(u)/(1+exp(u));
for(k in 2:12) {
        cj=cjStep[[k]]; u <- X%*%cj; p <- cbind(p,exp(u)/(1+exp(u)));
}    
matplot(zvals,p,col=c("black",rep("grey40",20)),type="l",lwd=c(2,rep(1,20)), xlab="Size z",
lty=1,ylab="Flowering probability",xlim=c(1.5,6)); 

points(zvals,p[,11],col="red",type="l",lwd=2); 
points(zvals,p[,12],col="blue",type="l",lwd=2); 
    
u1 <- -14.3 + zvals*m.par.true["flow.z"]; p1=exp(u1)/(1+exp(u1))
points(zvals,p1,type="p",col="blue"); 

dev.copy2eps(file="../../c10/figures/CarlinaSplineESS.eps"); 