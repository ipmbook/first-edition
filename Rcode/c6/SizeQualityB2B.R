# R code illustrating Bin-to-Bin implementation of a size-quality IPM
# using 5th order Gauss Legendre for initial and final bins  

rm(list=ls(all=TRUE))
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c6",sep="")); 

require(Matrix); require(statmod); 
source("domEig.R"); 
source("GaussQuadSub.R"); 

source("../Utilities/MatrixImage.R"); 

# Integrate a function FUN of 4 variables using specified weights and
# nodes in each dimension (matrices with ncol=4) 
quad4D <- function(FUN,weights,nodes) {
	order=nrow(nodes); int=0;  
    X=expand.grid(data.frame(nodes)); 
    W=expand.grid(data.frame(weights)); 
    fval=FUN(X[,1],X[,2],X[,3],X[,4]); 
    int=sum(fval*W[,1]*W[,2]*W[,3]*W[,4]); 
	return(int)
}

### tests: integrate function on [-1,1]^4
### because these are polynomials, the results should be exactly right 
FUN1 <- function(x,y,z,w) 1 # integral=16 
FUN2 <- function(x,y,z,w) x + y + z + w #integral=0 
FUN3 <- function(x,y,z,w) x^3 + y^3 + z^5 + w^5 # integral=0 
FUN4 <- function(x,y,z,w) x^2 + y^2 + z^2 + w^2 # integral=8 * (2/3)* 4 = 64/3
wts=nodes=matrix(0,7,4)
for(j in 1:4) {
	out=gaussQuadInt(-1,1,7)
	wts[,j]=out$weights;
	nodes[,j]=out$nodes
}
quad4D(FUN4,wts,nodes); # so far so good 

######## IPM Example 

vec <- function(nmat) matrix(nmat,ncol=1); 

#============================================================================ 
#  Define the kernel and iteration matrix
#============================================================================ 
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
      
# Compute cell centers for B2B iteration matrix
mx <- 25; mq <- 25; 
Lx <- (-1); Ux<- 7; Lq <- (-2.5); Uq <- (2.5); 
hx <- (Ux-Lx)/mx; yx <- Lx + hx*((1:mx)-0.5);
hq <- (Uq-Lq)/mq; yq <- Lq + hq*((1:mq)-0.5);

# Compute the 4D kernel and 2D iteration matrix. 

# Function eta to kernel values in their proper place in A 
eta_ij <- function(i,j,mx) {(j-1)*mx+i}

# matrix whose (i,j) entry is eta(ij) 
Eta <- outer(1:mx,1:mq,eta_ij,mx=mx); 

out <- gaussQuadInt(-hx/2,hx/2,order=5);
xWts <-  out$weights; xNodes <- out$nodes;
out <- gaussQuadInt(-hq/2,hq/2,order=5)
qWts <-out$weights; qNodes=out$nodes;
xqWeights=cbind(xWts,qWts,xWts,qWts)
A= matrix(0,mx*mq,mx*mq); Kvals=array(0,c(mx,mq,mx,mq));  

xqWeights=cbind(xWts,qWts,xWts,qWts)
system.time({
for(i in 1:mx){
	for(j in 1:mq){
		for(k in 1:mx){
			for(l in 1:mq){
				xqNodes=cbind(xNodes+yx[i],qNodes+yq[j],xNodes+yx[k],qNodes+yq[l])
				kvals=quad4D(FUN=k_xq,weights=xqWeights,nodes=xqNodes)/(hx*hx*hq*hq)
				A[Eta[i,j],Eta[k,l]]=kvals
				Kvals[i,j,k,l]=kvals
	}}}
	cat(i,"\n"); 
}
A<-hx*hq*A;  
})

#============================================================================# 
#  Repeat, but use the vectorised version. 
# !!! N.B. needs a huge amount of memory (~32GB) so don't run unless you have this
#============================================================================# 

system.time({

xorder <- length(xNodes)
qorder <- length(qNodes)

quad <- expand.grid(xp=yx, qp=yq, x=yx, q=yq, 
                    xpmap=seq.int(xorder), qpmap=seq.int(qorder),
                    xmap =seq.int(xorder), qmap =seq.int(qorder))

quad <- transform(quad, 
                  xpnodes = xNodes[xpmap], qpnodes = qNodes[qpmap],
                  xnodes  = xNodes[ xmap], qnodes  = qNodes[ qmap],
                  weights = xWts[xmap] * qWts[qmap] * xWts[xpmap] * qWts[qpmap])

quad$xmap <- NULL; quad$qmap <- NULL; quad$xpmap <- NULL; quad$qpmap <- NULL 

A.noloop <- with(quad, k_xq(xp+xpnodes, qp+qpnodes, x+xnodes, q+qnodes) * weights)

dim(A.noloop) <- c(mx * mq, mx * mq, xorder^2 * qorder^2)
A.noloop <- apply(A.noloop, c(1,2), sum) / (hx*hq)

})

all(round(A.noloop, 1e-16) == round(A, 1e-16))

#============================================================================# 
# Repeat, but use the vectorised version + sparse grid for quadrature
# This is pretty fast and uses a sensible amount of memory
#============================================================================# 

out <- createIntegrationGrid(type="GQU", dimension=4, k=5, sym = FALSE)
nodes <- out$nodes
nodes[,1] <- hx*(nodes[,1]-1/2)
nodes[,2] <- hq*(nodes[,2]-1/2)
nodes[,3] <- hx*(nodes[,3]-1/2)
nodes[,4] <- hq*(nodes[,4]-1/2)
weights <- out$weights

system.time({
quad <- expand.grid(xp=yx, qp=yq, x=yx, q=yq, map=seq_along(weights))
quad <- transform(quad, 
                  xpnodes = nodes[map,1], qpnodes = nodes[map,2],
                  xnodes  = nodes[map,3], qnodes  = nodes[map,4],
                  weights = out$weights[map])
A.noloop <- with(quad, k_xq(xp+xpnodes, qp+qpnodes, x+xnodes, q+qnodes) * weights)
rm(quad); gc()
dim(A.noloop) <- c(mx * mq, mx * mq, length(weights))
A.noloop <- apply(A.noloop, c(1,2), sum) * hx * hq
})

domEig(apply(A.noloop, c(1,2), sum))[[1]]

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

