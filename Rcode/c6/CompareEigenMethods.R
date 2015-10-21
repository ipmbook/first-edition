matmul <- function(x, A) { (A %*% x) }

domEig=function(A,tol=1e-8) {
	qmax=10*tol; lam=1; x=rep(1,nrow(A));   
	while(qmax>tol) {
		x1=matmul(x,A);
		qmax=sum(abs(x1-lam*x));  
		lam=sum(x1); 
		x=x1/lam; 
	} 
	return(list(lambda=lam,w=x/sum(x)))
}  	

n=25; 
A=matrix(runif(n^4),n^2,n^2); A=A/n^2; 
diag(A) <- 1 + diag(A)+runif(n^2);
system.time(domEig(A)); 

library(igraph); 
matmul <- function(x, A) { (A %*% x) }
system.time(arpack(matmul, extra=A, options=list(n=nrow(A), nev=1)))
# system.time(try(eigen(A)$values[1])); 

A=Matrix(A); matmul <- function(x, A) { (A %*% x)@x }
system.time(domEig(A));
system.time(arpack(matmul, extra=A, options=list(n=nrow(A), nev=1)))