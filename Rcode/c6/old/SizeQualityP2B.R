# R code illustrating implementation of a size x quality IPM 

rm(list=ls(all=TRUE))
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c6",sep="")); 

require(Matrix); require(statmod); 
source("domEig.R"); source("GaussQuadSub.R"); 
source("../Utilities/MatrixImage.R"); 

# Integrate a function FUN of 4 variables with respect to its first 
# two arguments, using specified weights and nodes, at all values in  
# a vector of its third argument, at a specified value of its fourth. 
quad4Dp <- function(FUN,xWts,qWts,xNodes,qNodes,xvals,q) {
	order=length(xNodes);  
    X=expand.grid(xNodes,qNodes,xvals);
    W=expand.grid(xWts,qWts,xvals); 
    fvals=FUN(X[,1],X[,2],X[,3],q);
    terms=matrix(fvals*W[,1]*W[,2],nrow=order^2); 
    return(colSums(terms))
}

######## IPM Example 

vec <- function(nmat) matrix(nmat,ncol=1); 

#============================================================================ 
#  Define the kernel and iteration matrix
#============================================================================ 
# Functions that make up the kernels  
g_x <- function(xp,x,q) dnorm(xp,m=1 + q + 0.7*x,sd=0.25)  #Growth   
g_q <- function(qp,q) dnorm(qp,m=0.5*q,sd=0.5)  #Quality dynamics
s_x <- function(x) exp(x-1)/(1+exp(x-1)) # Survival 
c_x <- function(xp,x) dnorm(xp,mean=0.8+x/5,sd=0.3) #Offspring size
c_q <- function(qp) dnorm(qp,mean=0,sd=0.25) #Offspring quality 

# Kernel functions 
p_xq <- function(xp,qp,x,q) s_x(x)*g_x(xp,x,q)*g_q(qp,q)
f_xq <- function(xp,qp,x,q) 0.75*s_x(x)*c_x(xp,x)*c_q(qp) 
k_xq <- function(xp,qp,x,q) p_xq(xp,qp,x,q)+f_xq(xp,qp,x,q)      
      
# Compute cell centers for P2B iteration matrix
mx <- 50; mq <- 30; 
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

out <- gaussQuadInt(-hx/2,hx/2,order=7)
xWts <-  out$weights; xNodes <- out$nodes;
out <- gaussQuadInt(-hq/2,hq/2,order=7)
qWts <-out$weights; qNodes=out$nodes;

A= matrix(0,mx*mq,mx*mq); Kvals=array(0,c(mx,mq,mx,mq));  
for(i in 1:mx){
	for(j in 1:mq){
			for(l in 1:mq){
				kvals=quad4Dp(FUN=k_xq,xWts,qWts,xNodes+yx[i],qNodes+yq[j],yx,q=yq[l])/(hx*hq)
				A[Eta[i,j],Eta[,l]]=kvals
				Kvals[i,j,,l]=kvals
	}}
	cat(i,"\n"); 
}
A<-hx*hq*A; 

#============================================================================# 
#  Find lambda, w,v by iteration 
#============================================================================# 
out<-domEig(A); out2<-domEig(t(A)); 
lam.stable=out$lambda; lam.stable.t <- out2$lambda; 
w <- Re(matrix(out$w,mx,mq)); 
w <- w/(hx*hq*sum(w)); 
v <- Re(matrix(out2$w,mx,mq)); v <- v/sum(v); 

# Compute elasticity matrix 
repro.val=matrix(v,mx,mq); stable.state=matrix(w,mx,mq); 
v.dot.w=sum(hx*hq*stable.state*repro.val)
sens=outer(repro.val,stable.state)/v.dot.w
elas=sens*Kvals/lam.stable;

# Compute matrix of total(=integrated) elasticities for all transitions (x_i,q_j) -> anywhere 
total.elas=hx*hq*apply(elas,c(3,4),sum); 

# Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hx*hq*hx*hq*elas)," and it should = 1","\n"); 

#============================================================================# 
#  Plot stable state distribution and state-dependent total elasticity 
#============================================================================# 
source("../Utilities/Standard Graphical Pars.R"); 
graphics.off(); dev.new(); set_graph_pars("panel4"); 

plot(yx,apply(stable.state,1,sum),xlab="Size x",ylab="frequency",type="l"); 
# title(main="Stable size distribution",cex.main=0.8); 
add_panel_label("a"); 

#plot(yq,apply(stable.state,2,sum),xlab="Quality q",ylab="frequency",type="l"); 
# title(main="Stable quality distribution",cex.main=0.8); 
#add_panel_label("b"); 

matrix.image(stable.state,yq,yx,xlab="Quality q",ylab="Size x",do.legend=FALSE);
title(main="Size-quality distribution",cex.main=0.9); 
add_panel_label("b"); 

matrix.image(repro.val,yq,yx,xlab="Quality q",ylab="Size x",do.legend=FALSE);
title(main="Reproductive value",cex.main=0.9); 
add_panel_label("c"); 

matrix.image(total.elas,yq,yx,xlab="Quality q",ylab="Size x",do.legend=FALSE);
title(main="Total elasticity",cex.main=0.9); 
add_panel_label("d"); 

#dev.copy2eps(file="../../c6/figures/SizeQualityEx.eps"); 
