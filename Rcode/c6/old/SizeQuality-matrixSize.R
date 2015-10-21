# R code illustrating implementation of a size x quality IPM 

rm(list=ls(all=TRUE))
require(Matrix); 
setwd("c:/repos/ipm_book/Rcode/c6"); 
source("domEig.R"); source("GaussQuadSub.R"); 

vec <- function(nmat) matrix(nmat,ncol=1); 

#============================================================================ 
#  Define the kernel and iteration matrix
#============================================================================ 
# Functions that make up the kernels  
g_x <- function(xp,x,q) dnorm(xp,m=1 + q + 0.7*x,sd=0.05)  #Growth   
g_q <- function(qp,q) dnorm(qp,m=0.5*q,sd=0.2)  #Quality dynamics
s_x <- function(x) exp(x-1)/(1+exp(x-1)) # Survival 
c_x <- function(xp,x) dnorm(xp,mean=0.8+x/5,sd=0.1) #Offspring size
c_q <- function(qp) dnorm(qp,mean=0,sd=0.5) #Offspring quality 

# Kernel functions 
p_xq <- function(xp,qp,x,q) s_x(x)*g_x(xp,x,q)*g_q(qp,q)
f_xq <- function(xp,qp,x,q) 0.75*s_x(x)*c_x(xp,x)*c_q(qp) 
k_xq <- function(xp,qp,x,q) p_xq(xp,qp,x,q)+f_xq(xp,qp,x,q)      

# Integrate a function FUN of 4 variables with respect to its first 
# two arguments, using specified weights and nodes (matrices with ncol=2), 
# at all values in a vector of its third argument and one value of its fourth
quad4Dp <- function(FUN,xWts,qWts,xNodes,qNodes,xvals,q) {
	order=length(xNodes);  
    X=expand.grid(xNodes,qNodes,xvals);
    W=expand.grid(xWts,qWts,xvals); 
    fvals=FUN(X[,1],X[,2],X[,3],q);
    terms=matrix(fvals*W[,1]*W[,2],nrow=order^2); 
    return(colSums(terms))
}

mvals=c(seq(20,60,by=5),100); lvals=matrix(0,length(mvals),2);  
for(mj in 1:length(mvals)){       
# Compute meshpoints: a small matrix so this runs fast 
mx <- mvals[mj]; mq <- mx/2; 
Lx <- (-1); Ux<- 7; Lq <- (-2.5); Uq <- (2.5); 
hx <- (Ux-Lx)/mx; yx <- Lx + hx*((1:mx)-0.5);
hq <- (Uq-Lq)/mq; yq <- Lq + hq*((1:mq)-0.5);

# Compute the 4D kernel and 2D iteration matrix. With some simple 
# vectorizing it's not too slow. The shortcuts here have 
# been validated against code that uses loops for everything. 

# Function eta to kernel values in their proper place in A 
eta_ij <- function(i,j,mx) {(j-1)*mx+i}

# matrix whose (i,j) entry is eta(ij) 
Eta <- outer(1:mx,1:mq,eta_ij,mx=mx); 

A=matrix(0,mx*mq,mx*mq); Kvals=array(0,c(mx,mq,mx,mq));  
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

out <- gaussQuadInt(-hx/2,hx/2,order=9)
xWts <-  out$weights; xNodes <- out$nodes;
out <- gaussQuadInt(-hq/2,hq/2,order=9)
qWts <-out$weights; qNodes=out$nodes;

if(mx < 70) {
A2 = matrix(0,mx*mq,mx*mq); Kvals2=array(0,c(mx,mq,mx,mq));  
for(i in 1:mx){
	for(j in 1:mq){
			for(l in 1:mq){
				kvals=quad4Dp(FUN=k_xq,xWts,qWts,xNodes+yx[i],qNodes+yq[j],yx,q=yq[l])/(hx*hq)
				A2[Eta[i,j],Eta[,l]]=kvals
				Kvals2[i,j,,l]=kvals
	}}
	cat(i,"Point 2 Bin","\n"); 
}
A2<-hx*hq*A2; out<-domEig(A2); lvals[mj,2]=out$lambda; 
cat(mx,out$lambda,"\n"); 
}

}


