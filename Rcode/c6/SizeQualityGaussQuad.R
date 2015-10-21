# R code illustrating implementation of a size x quality IPM 
# Using Gauss-Legendre quadrature 

rm(list=ls(all=TRUE))
require(Matrix); 
setwd("c:/repos/ipm_book/Rcode/c6"); 
source("domEig.R"); source("GaussQuadSub.R"); 
source("../Utilities/MatrixImage.R"); 

vec <- function(nmat) matrix(nmat,ncol=1); 

#============================================================================ 
#  Define the kernel and iteration matrix
#============================================================================ 
# Survival probability is logistic regression s(x)=exp(x-1)/1+exp(x-1)
# New size = Normal(mean=1 + quality + 0.7*old size,sd=0.3)
# New quality = Normal(mean=.5*old quality,sd=.5)  
# Fecundity = 0.75 per surviving individual
# Offspring size = Normal(mean=0.8+0.2*parent size, sd=0.35)
# Offspring quality=Normal(mean=0,sd=.5)
      
# Functions that make up the kernels  
g_x <- function(xp,x,q) dnorm(xp,m=1 + q + 0.7*x,sd=0.3)  #Growth   
g_q <- function(qp,q) dnorm(qp,m=0.5*q,sd=0.5)  #Quality dynamics
s_x <- function(x) exp(x-1)/(1+exp(x-1)) # Survival 
c_x <- function(xp,x) dnorm(xp,mean=0.8+x/5,sd=0.35) #Offspring size
c_q <- function(qp) dnorm(qp,mean=0,sd=0.5) #Offspring quality 

# Kernel functions 
p_xq <- function(xp,qp,x,q) s_x(x)*g_x(xp,x,q)*g_q(qp,q)
f_xq <- function(xp,qp,x,q) 0.75*s_x(x)*c_x(xp,x)*c_q(qp) 
k_xq <- function(xp,qp,x,q) p_xq(xp,qp,x,q)+f_xq(xp,qp,x,q)      
      
# Compute meshpoints and weights 
mx <- 50; mq <- 25; Lx <- (-1); Ux<- 7; Lq <- (-2.5); Uq <- (2.5); 
out=gaussQuadInt(Lx,Ux,mx); yx=out$nodes; Wx=out$weights; 
out=gaussQuadInt(Lq,Uq,mq); yq=out$nodes; Wq=out$weights; 

# Compute the 4D kernel and 2D iteration matrix. 

# Function eta to kernel values in their proper place in A 
eta_ij <- function(i,j,mx) {(j-1)*mx+i}

# matrix whose (i,j) entry is eta(ij) 
Eta <- outer(1:mx,1:mq,eta_ij,mx=mx); 

A=matrix(0,mx*mq,mx*mq); W=Kvals=array(0,c(mx,mq,mx,mq));  
for(i in 1:mx){
	for(j in 1:mq){
		for(k in 1:mq){
				kvals=k_xq(yx,yq[k],yx[i],yq[j])
				A[Eta[,k],Eta[i,j]]=kvals*Wx[i]*Wq[j]
                W[,k,i,j]=Wx*Wq[k]*Wx[i]*Wq[j] #weight for 4D integration 
				Kvals[,k,i,j]=kvals
	}}
	cat(i,"\n"); 
}

#============================================================================# 
#  Find lambda, w,v by iteration 
#============================================================================# 
out<-domEig(A); out2<-domEig(t(A)); 
lam.stable=out$lambda; lam.stable.t <- out2$lambda; 
w <- Re(matrix(out$w,mx,mq)); 
w <- w/sum(outer(Wx,Wq)*w); 
v <- Re(matrix(out2$w,mx,mq)); v <- v/sum(v); 

# Compute elasticity matrix 
repro.val=matrix(v,mx,mq); stable.state=w; 
v.dot.w=sum(outer(Wx,Wq)*stable.state*repro.val)
sens=outer(repro.val,stable.state)/v.dot.w
elas=sens*Kvals/lam.stable;

# Compute matrix of total(=integrated) elasticities for all transitions (x_i,q_j) -> anywhere 
total.elas=matrix(0,mx,mq); 
for(i in 1:mx) {
  for(j in 1:mq) {
    total.elas[i,j]=sum(elas[,,i,j]*outer(Wx,Wq))
}}

# Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(W*elas)," and it should = 1","\n"); 

#============================================================================# 
#  Plot stable state distribution and state-dependent total elasticity 
#============================================================================# 
source("../Utilities/Standard Graphical Pars.R"); 
graphics.off(); dev.new(); set_graph_pars("panel4"); 

plot(yx,apply(stable.state,1,function(x) sum(x*Wq)),xlab="Size x",ylab="frequency",type="l"); 
# title(main="Stable size distribution",cex.main=0.8); 
add_panel_label("a"); 

matrix.image(stable.state,yq,yx,xlab="Quality q",ylab="Size x",do.legend=FALSE);
title(main="Size-quality distribution",cex.main=0.9); 
add_panel_label("b"); 

matrix.image(repro.val,yq,yx,xlab="Quality q",ylab="Size x",do.legend=FALSE);
title(main="Reproductive value",cex.main=0.9); 
add_panel_label("c"); 

matrix.image(total.elas,yq,yx,xlab="Quality q",ylab="Size x",do.legend=FALSE);
title(main="Total elasticity",cex.main=0.9); 
add_panel_label("d"); 

#dev.copy2eps(file="../../c6/figures/SizeQualityGQ.eps"); 

#============================================================================# 
#  Compute the Fundamental Matrix 
#============================================================================# 
P=F=matrix(0,mx*mq,mx*mq); 
for(i in 1:mx){
	for(j in 1:mq){
		for(k in 1:mq){
				P[Eta[,k],Eta[i,j]]=p_xq(yx,yq[k],yx[i],yq[j])*Wx[i]*Wq[j]
                F[Eta[,k],Eta[i,j]]=f_xq(yx,yq[k],yx[i],yq[j])*Wx[i]*Wq[j]
	}}
	cat(i,"\n"); 
}

N <- solve( diag(rep(1),nrow(P))-P); ## fundamental matrix
Times <- N%*%(F%*%matrix(stable.state,ncol=1)); 

