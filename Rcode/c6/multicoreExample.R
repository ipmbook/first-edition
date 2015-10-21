#### using mclapply, on a mac 
### make the F[,,a] matrices into a list
f=array(rnorm(250*250*6),c(250,250,6))
flist=as.list(1:6); for(j in 1:6) flist[[j]]=f[,,j]

### make the n[,a] matrices into a list
n=matrix(rnorm(250*6),250,6)
nlist=as.list(1:6); for(j in 1:6) nlist[[j]]=n[,j]

### define a function such that f(i) = the age-i Fn multiplication 
fmult=function(i) flist[[i]]%*%nlist[[i]]
out1=mclapply(1:6,fmult,mc.cores=2)

### PC 
require(parallel); 
cl=makeCluster(2); 
fmult=function(i) f[,,i]%*%n[,i]

f=array(rnorm(200*200*100),c(200,200,100))
n=matrix(rnorm(200*100),200,100);
nt=array(0,c(200,100,500)); 
nt[,,1]=n; 
clusterExport(cl,varlist=c("f","fmult")); 

### multicore
system.time({
for(k in 1:499) {
n=nt[,,k]; clusterExport(cl,"n"); 
out2=parLapply(cl,1:100,fmult); 
for(j in 1:100) nt[,j,k+1]=out2[[j]]; 
}
})

### non-multicore
ntz=nt; 
system.time({
for(k in 1:499) {
for(j in 1:100) ntz[,j,k+1]=f[,,j]%*%ntz[,j,k]
}})

#   user  system elapsed 
# 114.64    1.30  115.97 

