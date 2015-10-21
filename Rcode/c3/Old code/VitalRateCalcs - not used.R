#mk.K <- function(m,m.par,L,U) 

Mtrue <- mk.K(100,m.par.true,-3,8); 
Pmat <- Mtrue$P
Fmat <- Mtrue$F
Kmat <- Mtrue$K
z <- Mtrue$meshpts
e<-rep(1,dim(Kmat)[1]); 

f2 <- e%*%Fmat%*%Pmat%*%Pmat; 
plot(z,f2); 