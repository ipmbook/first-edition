## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Find the ESS function-valued flowering strategy and compare to the ESS logistic regression strategy
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))
require(fda); 

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c9",sep="")); 

source("../c2/Monocarp Demog Funs.R");
source("Monocarp Demog SplineFuns.R"); 
source("../utilities/Standard Graphical Pars.R");

# Just to be safe, we re-specify the parameter vector for Oenothera 
m.par.true <- c(## survival
                surv.int  =  -0.65,
                surv.z    =   0.75,
                ## flowering
                flow.int  = -18.00,
                flow.z    =   6.9,
                ## growth
                grow.int  =   0.96,
                grow.z    =   0.59,
                grow.sd   =   0.67,
                ## recruit size
                rcsz.int  =   -.08, 
                rcsz.sd   =   0.76,
                ## seed size
                seed.int  =   1.00,
                seed.z    =   2.20,
                ## recruitment probability
                p.r       =   0.007)  

L <- (-6); U <- 7; ## WIDE size range

####################################################################
# For comparison, find ESS with logist regression flowering 
####################################################################

# Function to compute R0-tilde
Rtilde0 <- function(beta0,parms) {
    pars <- parms; pars["flow.int"] <- beta0; 
    IPM <- mk_K(100,m.par=pars,L=L,U=U) 
    N <- solve(diag(100)-IPM$P);
    R0 <- abs(eigen(IPM$F%*%N)$values[1])
    return(R0)
}    

### confirm that we replicate previous results 
betaVals <- R0vals <- seq(-30,-20,length=30); 
for(j in 1:30) R0vals[j]=Rtilde0(betaVals[j],m.par.true)
plot(betaVals,R0vals);  

## find the ESS beta0 exactly 
out <- optimize(function(b) -Rtilde0(b,m.par.true),lower=-30,upper=-20); 
beta0.ESS <- out$minimum
R0.ESS <- -out$objective; 


abline(v=beta0.ESS); # So far so good. 

####################################################################
# Now find the function-valued ESS 
####################################################################
L <- (-6); U <- 7; ## WIDE size range
B <- create.bspline.basis(rangeval=c(L,U),nbasis=12, norder=4)
X <- eval.basis(seq(L,U,length=250),B); dz=(U-L)/250; # for computing slope 
source("Monocarp Demog SplineFuns.R"); 

# Function to compute R0-tilde with spline flowering 
Rtilde0_spline <- function(cj,parms,B) {
    IPM <- mk_K_spline(100,m.par=parms,L=L,U=U,B,cj) 
    N <- solve(diag(100)-IPM$P);
    R0 <- abs(eigen(IPM$F%*%N)$values[1])
    return(R0)
}  
Rtilde0_spline(runif(12),m.par.true,B);  # test 

objfun=function(cj,parms,B) {
    #check: is max slope too large? 
    slope <- diff(X%*%cj)/dz; 
    if(max(slope) > 0.5+parms["flow.z"]) {
        return(1e12)
    }else{
    R0 <- Rtilde0_spline(cj,parms,B) 
    return(-R0) 
    }    
}   
p0=runif(12); objfun(p0,m.par.true,B); # test

pbfit <- optim(par=p0,fn=objfun,control=list(trace=4,maxit=2500),parms=m.par.true,B=B); 
pbfit <- optim(par=pbfit$par,fn=objfun,control=list(trace=4,maxit=2500),parms=m.par.true,B=B); 
cj.opt <- pbfit$par; R0.opt <- pbfit$value; 

dev.new(height=6,width=8);
set_graph_pars("panel1"); 
zvals=seq(1,5,length=200); 
p_b.opt <- p_bz_spline(zvals,B,cj.opt);
plot(zvals,p_b.opt,type="l",lwd=2,xlab="Size z",ylab="Flowering probability"); 

## compare 
m.par <- m.par.true; m.par["flow.int"] <- beta0.ESS; 
zvals=seq(1,5,length=50); 
p_b.ESS <- p_bz(zvals, m.par)
points(zvals,p_b.ESS,type="p",col="grey30",cex=1.4,lwd=2,pch=19); 

#dev.copy2eps(file="../../c9/figures/OenotheraSplineFlower.eps"); 
