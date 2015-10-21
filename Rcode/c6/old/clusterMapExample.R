cl=makeCluster(6); 
clusterExport(cl,list=ls());

XQ=expand.grid(1:mx,1:mq); 
out=clusterMap(cl, fun=k_xq, xp=XQ[,1], qp=XQ[,2],MoreArgs=c(x=yx[40],q=yx[20]),SIMPLIFY=TRUE );
