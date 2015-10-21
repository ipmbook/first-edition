require(Matrix); require(igraph);

n=50; A=40; set.seed(10661492)
F=matrix(runif(n*n),n,n)/n;
diag(F) <- diag(F)+runif(n,2,3) 
Q=diag(rep(1,A));
M=kronecker(Q,F)
M.s<- Matrix(M,sparse=TRUE)

## one big matrix 
n=runif(nrow(M));
system.time(for(j in 1:100) n <- M%*%n)
system.time(for(j in 1:100) n <- M.s%*%n)

# loop 
x=runif(nrow(F))
system.time(for(j in 1:100){for (i in 1:A) x<-F%*%x})

system.time(domEig(M))
system.time(domEig(M.s))

matmul <- function(x, A) { (A %*% x) }
system.time(arpack(matmul, extra=M, options=list(n=nrow(M), nev=1)))
system.time(arpack(matmul, extra=M.s, options=list(n=nrow(M.s), nev=1)))