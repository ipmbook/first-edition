##########################################################################
### Functions to compute nodes and weights for Gauss-Legendre quadrature
### on one interval, or on a set of subintervals within an interval 
##########################################################################
require(statmod); 

# Gauss-Legendre quadrature on interval (L,U) 
gaussQuadInt <- function(L,U,order=7) {
	# nodes and weights on [-1,1]
	out <- gauss.quad(order); #GL is the default 
    w <- out$weights; x <- out$nodes;  
    weights=0.5*(U-L)*w; 
    nodes=0.5*(U+L) + 0.5*(U-L)*x; 
	return(list(weights=weights,nodes=nodes)); 
}

### For comparison, midpoint rule
midpointInt <- function(L,U,m) {
    h=(U-L)/m; nodes=L+(1:m)*h-h/2;
    weights=rep(h,m); 
    return(list(weights=weights,nodes=nodes)); 
}   

## Gaussian-Legendre quadrature on subintervals of (L,U). 
## User can specify the number of subintervals (intervals)
## or the locations of breaks between subintervals (breaks) 
gaussQuadSub <- function(L,U,order=7,intervals=1,breaks=NULL) {
	# nodes and weights on [-1,1]
	out <- gauss.quad(order); w <- out$weights; x <- out$nodes;  

	# compute subinterval endpoints 
	if(is.null(breaks)){
		h <- (U-L)/intervals;
		b <- L + (1:intervals)*h;
		a <- b-h; 
	} else {
		intervals <- 1+length(breaks);
		a <- c(L,breaks);
		b <- c(breaks,U);
	} 
    weights=as.vector(sapply(1:intervals,function(j) 0.5*(b[j]-a[j])*w))
    nodes=as.vector(sapply(1:intervals,function(j) 0.5*(a[j]+b[j]) + 0.5*(b[j]-a[j])*x))
	return(list(weights=weights,nodes=nodes)); 
}

# Integrate a function FUN of 2 variables using specified 
# weights and nodes in each variable. FUN must accept two 
# vectors as arguments and return a vector of function values
quad2D <- function(FUN,wts1,wts2,nodes1,nodes2) { 
    X=expand.grid(nodes1,nodes2); 
    W=expand.grid(wts1,wts2); 
    fval=FUN(X[,1],X[,2]); 
    int=sum(fval*W[,1]*W[,2]); 
	return(int)
}