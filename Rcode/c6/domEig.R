domEig=function(A,tol=1e-8) {
	qmax=10*tol; lam=1; x=rep(1,nrow(A));   
	while(qmax>tol) {
		x1=A%*%x;
		qmax=sum(abs(x1-lam*x));  
		lam=sum(x1); 
		x=x1/lam; 
	} 
    # Having found w (within tol), get lambda 
    x1 = A%*%x; lam=sum(x1); x=x1/lam;   
	return(list(lambda=lam,w=x/sum(x)))
}  	
