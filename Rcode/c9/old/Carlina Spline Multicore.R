rm(list=ls(all=TRUE));
set.seed(53241986)
## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep="")); 

source("../utilities/Standard Graphical Pars.R");
source("Carlina Demog Funs.R"); # Carlina functions and parameters 
source("Monocarp Demog SplineFuns.R");  

require(parallel); 
cl=makeCluster(6,useXDR=FALSE);
clusterEvalQ(cl,require(fda)); 
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
        F <- h * (outer(meshpts, meshpts, F_z1z_spline, m.par = m.par.year, B = B, cj = cj))
        
        seeds=F%*%nt; 
        total.seeds=sum(h*seeds); 
        p.r.yr <- Recr.year/total.seeds; 
        
        nt <- seeds*p.r.yr + P%*%nt; 
        params.yr[gen-1,"p.r"] <- p.r.yr;
        Nt[gen] <- sum(h*nt); 
 	
        if(gen%%10==0) cat(gen,p.r.yr,sum(h*nt),"\n")
}

return(list(nt=nt,Nt=Nt,params.yr=params.yr))
}

spline_invader_gr <- function(params.yr,B,cj,m,L,U) {
    n.iter=nrow(params.yr); 
    h <- (U-L)/m; meshpts <- L+ ((1:m)-0.5)*h
    nt <- matrix(1/m,m,1); Rt <- numeric(n.iter);  
 
	for(gen in 2:n.iter){
        m.par.year <- params.yr[gen-1,]; 
        P <- h * (outer(meshpts, meshpts, P_z1z_spline, m.par = m.par.year, B = B, cj = cj))
        F <- h * (outer(meshpts, meshpts, F_z1z_spline, m.par = m.par.year, B = B, cj = cj))
        nt <- (P+F)%*%nt; 
        Rt[gen] <- sum(nt); nt<-nt/Rt[gen]; 
        if(gen%%10==0) cat(gen,Rt[gen],"\n")
}
		lambda.hat = mean(log(Rt[101:n.iter]))
return(list(nt=nt,Rt=Rt,lambda.hat=lambda.hat))
}


### Objective function for maximizing invader lambda 
objfun.inv=function(cj,params.yr,B,m,L,U,max.slope,X) {
    slope <- diff(X%*%cj)/dz; 
	pen=ifelse(max(slope) > max.slope, 1e12*(max(slope)-max.slope)^2, 0) 
    lam<-spline_invader_gr(params.yr=params.yr,B=B,cj=cj,m=100,L=L,U=U)$lambda.hat
    return(-lam+pen)
} 

d.objfun.inv <- function(cj,...) {
		npar=length(cj)
	    mat=diag(1+npar)
        gfun=function(i) {
            return(objfun.inv(cj+0.01*mat[1:npar,i],params.yr,B,m,L,U,max.slope,X))
        }    
        out=unlist(clusterApply(cl,x=1:(1+npar),fun=gfun))
        return(100*(out[1:npar]-out[1+npar])); 
}        

clusterExport(cl,objects())

## create spline basis 
L <- (-3); U <- 8; ## WIDE size range
B <- create.bspline.basis(rangeval=c(L,U),nbasis=8, norder=4)
X <- eval.basis(seq(L,U,length=250),B); dz=(U-L)/250; # for computing slope 

## initialize parameters to approximate estimated true flowering strategy (not an ESS!) 
cj <- runif(5); 
zvals <- seq(L,U,length=250); 
yvals <- m.par.true["flow.int"] + zvals*m.par.true["flow.z"]; 
mse <- function(cj) { sqrt(mean((yvals - X%*%cj)^2)) }
out <- optim(par=cj,f=mse,control=list(trace=4,maxit=2500)); 
out <- optim(par=out$par,f=mse,control=list(trace=4,maxit=2500)); 
cj <- out$par;

## Test: does the model run? 
res <- iterate_spline_resident(params=m.par.true,n.iter=600,B=B,cj=cj,m=80,L=L,U=U);   
inv <- spline_invader_gr(params.yr=as.matrix(res$params.yr),B=B,cj=cj,m=80,L=L,U=U) 

params.yr=as.matrix(res$params.yr); m=80; mat=diag(6); max.slope=4;  
out=d.objfun.inv(cj,params.yr,B,m,L,U,max.slope=4,X)

### simulate sequential invasion
cjStep=list(10); cjStep[[1]]=cj;  
for(j in 2:10) { 
    res <- iterate_spline_resident(params=m.par.true,n.iter=600,B=B,cj=cjStep[[j-1]],m=80,L=L,U=U); 
    clusterExport(cl,varlist=objects()); 
    out=optim(par=cjStep[[j-1]],f=objfun.inv,gr=d.objfun.inv,method="BFGS",control=list(maxit=5,trace=4),params.yr=as.matrix(res$params.yr),
            B=B,m=80,L=L,U=U,max.slope=0.05+m.par.true["flow.z"])
    cjStep[[j]]=out$par
    cat("Invasion step ", j, "completed")
    cj=cjStep[[1]]; u <- X%*%cj; p <- exp(u)/(1+exp(u));
    plot(zvals,p,type="l",lty=1,col=1); 
    for(k in 2:j) {
        cj=cjStep[[k]]; u <- X%*%cj; p <- exp(u)/(1+exp(u));
        points(zvals,p,type="l",lty=1,col=k); 
    }    
}    

# find a good invader against the resident 
out=optim(par=cj,f=objfun.inv,control=list(maxit=100,trace=4),params.yr=as.matrix(res$params.yr),
            B=B,m=100,L=L,U=U,max.slope=)
cj2 <- out$par; 
u1 <- X%*%cj; u2 <- X%*%cj2; 
p1 <- exp(u1)/(1+exp(u1));
p2 <- exp(u2)/(1+exp(u2));
matplot(zvals,cbind(p1,p2),col=c("black","blue"),type="l",lty=1,xlab="z",ylab="Flowering probability");

# let the new invader become the new resident
res <- iterate_spline_resident(params=m.par.true,n.iter=800,B=B,cj=cj2,m=100,L=L,U=U);
# find a good invader against the resident 
out=optim(par=cj2,f=objfun.inv,control=list(maxit=100,trace=4),params.yr=as.matrix(res$params.yr),
            B=B,m=100,L=L,U=U,max.slope=max.slope)
            
cj3 <- out$par; 
u3 <- X%*%cj3; 
p3 <- exp(u3)/(1+exp(u3));
matplot(zvals,cbind(p1,p2,p3),col=c("black","blue","red"),type="l",lty=1,xlab="z",ylab="Flowering probability");           
  