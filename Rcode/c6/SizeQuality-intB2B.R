# R code illustrating implementation of a size x quality IPM 
# using Bin-to-Bin and Cumulative Kernel, compared with midpoint
# rule. 

rm(list=ls(all=TRUE))
require(Matrix); 
setwd("c:/repos/ipm_book/Rcode/c6"); 
source("domEig.R"); source("GaussQuadSub.R"); 

# Integrate a function FUN of 2 variables using specified weights and
# nodes in each variable
quad2D <- function(FUN,wts1,wts2,nodes1,nodes2) { 
    X=expand.grid(nodes1,nodes2); 
    W=expand.grid(wts1,wts2); 
    fval=FUN(X[,1],X[,2]); 
    int=sum(fval*W[,1]*W[,2]); 
	return(int)
}

vec <- function(nmat) matrix(nmat,ncol=1); 

#============================================================================ 
#  Define the kernel 
#============================================================================ 
# Functions that make up the kernels  
g_x <- function(xp,x,q) dnorm(xp,m=1 + q + 0.7*x,sd=0.3)  #Growth   
g_q <- function(qp,q) dnorm(qp,m=0.5*q,sd=0.2)  #Quality dynamics
s_x <- function(x) exp(x-1)/(1+exp(x-1)) # Survival 
c_x <- function(xp,x) dnorm(xp,mean=0.8+x/5,sd=0.1) #Offspring size
c_q <- function(qp) dnorm(qp,mean=0,sd=0.5) #Offspring quality 

# Kernel functions 
p_xq <- function(xp,qp,x,q) s_x(x)*g_x(xp,x,q)*g_q(qp,q)
f_xq <- function(xp,qp,x,q) 0.75*s_x(x)*c_x(xp,x)*c_q(qp) 
k_xq <- function(xp,qp,x,q) p_xq(xp,qp,x,q)+f_xq(xp,qp,x,q)

#============================================================================ 
#  Define the bin-integrated kernel 
#============================================================================ 
intg_x <- function(xp,x,q,hx) pnorm(xp+hx/2,m=1 + q + 0.7*x,sd=0.02) - pnorm(xp-hx/2,m=1 + q + 0.7*x,sd=0.02)
intg_q <- function(qp,q,hq) pnorm(qp+hq/2,m=0.5*q,sd=0.2) - pnorm(qp-hq/2,m=0.5*q,sd=0.2)    
intp_xq <- function(xp,qp,x,q,hx,hq) s_x(x)*intg_x(xp,x,q,hx)*intg_q(qp,q,hq)

intc_x <- function(xp,x,hx) pnorm(xp+hx/2,mean=0.8+x/5,sd=0.1) - pnorm(xp-hx/2,mean=0.8+x/5,sd=0.1)
intc_q <- function(qp,hq) pnorm(qp+hq/2,mean=0,sd=0.5) - pnorm(qp-hq/2,mean=0,sd=0.5)
intf_xq <- function(xp,qp,x,q,hx,hq) 0.75*s_x(x)*intc_x(xp,x,hx)*intc_q(qp,hq) 
intk_xq <- function(xp,qp,x,q,hx,hq) intp_xq(xp,qp,x,q,hx,hq)+intf_xq(xp,qp,x,q,hx,hq)    

################################################################################ 
# Construct iteration matrices by BIN-TO-BIN and midpoint rule 
################################################################################
mvals=c(seq(20,60,by=10),1000); # number of bins for size 

lvals=matrix(0,length(mvals),2);   
for(mj in 1:length(mvals)){
       
# Compute meshpoints: a small matrix so this runs fast 
mx <- mvals[mj]; mq <- mx/2; 
Lx <- (-1); Ux<- 7; Lq <- (-2.5); Uq <- (2.5); 
hx <- (Ux-Lx)/mx; yx <- Lx + hx*((1:mx)-0.5);
hq <- (Uq-Lq)/mq; yq <- Lq + hq*((1:mq)-0.5);

# weights and nodes for GL(3) in initial state (x,q) 
out <- gaussQuadInt(-hx/2,hx/2,order=3)
xWts <-  out$weights; xNodes <- out$nodes;
out <- gaussQuadInt(-hq/2,hq/2,order=3)
qWts <-out$weights; qNodes=out$nodes;

# Function eta to kernel values in their proper place in A 
eta_ij <- function(i,j,mx) {(j-1)*mx+i}

# matrix whose (i,j) entry is eta(ij) 
Eta <- outer(1:mx,1:mq,eta_ij,mx=mx); 

### Iteration matrix by midpoint rule 
A=A2=matrix(0,mx*mq,mx*mq); Kvals=Kvals2=array(0,c(mx,mq,mx,mq));  
for(i in 1:mx){
	for(j in 1:mq){
		for(k in 1:mq){
				kvals=k_xq(yx,yq[k],yx[i],yq[j])
				A[Eta[,k],Eta[i,j]]=kvals
				Kvals[,k,i,j]=kvals
			
	}}
	cat(i,"Midpoint","\n"); 
}

A<-hx*hq*A; out<-domEig(A); lvals[mj,1]=out$lambda; 
cat(mx,out$lambda,"\n"); 

### Iteration matrix by B2B with cumulative kernel 
if(mx < 50) {
  for(i in 1:mx){ 
	for(j in 1:mq){
		for(k in 1:mx){
			for(l in 1:mq){
                FUN=function(x,q,theta) intk_xq(yx[i],yq[j],x,q,hx,hq); 
                kvals2=quad2D(FUN=FUN,xWts,qWts,xNodes+yx[k],qNodes+yq[l])/(hx*hq)
				A2[Eta[i,j],Eta[k,l]]=kvals2
				Kvals2[i,j,k,l]=kvals2
	}}}
	cat(i,"IntB2B","\n"); 
}

out<-domEig(A2); lvals[mj,2]=out$lambda; 
cat(mx,out$lambda,"\n"); 
}

}

# RESULTS
        [,1]     [,2]
#[1,] 1.737108 1.216370
#[2,] 2.292719 1.214917
#[3,] 1.612873 1.214767
#[4,] 1.199876 1.214709
#[5,] 1.220360 1.214670
#[6,] 1.214589 0.000000

