cl=makeCluster(4,useXDR=FALSE);
func=function(i,Nmin,Nmax) {
	return(mean(seq(Nmin,Nmax,by=0.001)^i)^(1/i))
}

out=unlist(parLapply(cl,1:25,func,Nmin=1,Nmax=4))
plot(1:25,out); 
out2=unlist(parLapply(cl,1:25,func,Nmin=2,Nmax=4))
points(1:25,out2); 

## Note how additional arguments to func are passed in larLapply. 
## This is faster than doing clusterExport, but depends on
## func not looking for any variables in its .GlobalEnv 