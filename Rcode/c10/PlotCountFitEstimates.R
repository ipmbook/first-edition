## Plotting results from the count data fitting 
rm(list=ls(all=TRUE))

## Working directory must be set here
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c10",sep=""));

X=read.csv(file="CountFitEstimates.csv"); 
names(X)<-c("Set","  a","  b","  sigma","  q","  delta","  nu","  tau","  Kf"); 

X.labels <- c(expression(italic(a)),expression(italic(b)),expression(sigma),
expression(italic(q)),expression(delta),expression(nu),expression(tau),expression(italic(K)[italic(f)])); 

X1=X[X$Set=="Early",-1]; X2=X[X$Set=="Late",-1]

a=0.9; b = 2; sigma=0.8; q=0.85; Delta=0.6; nu=.5; tau=.5; Kf=8; 
parms=matrix(c(a,b,sigma,q,Delta,nu,tau,Kf),1,8); 
parms=data.frame(parms); 
names(parms)<-c("  a","  b","  sigma","  q","  delta","  nu","  tau","  Kf"); 

dev.new(height=5,width=8);
par(cex.axis=1.2,cex.lab=1.3,mar=c(3,4,1,1),mgp=c(2.35,1,0),bty="l"); 
par(tck=0); 
boxplot(log(X1),ylim=c(-1.3,3.1),xlim=c(0,24),at=c(1,4,7,10,13,16,19,22),ylab="Log(estimate)",names=NA);
boxplot(log(X2),ylim=c(-1.3,3.1),xlim=c(0,24),at=c(1,4,7,10,13,16,19,22)+1,add=TRUE,col="grey75",names=NA);
points(c(1,4,7,10,13,16,19,22)-0.75, log(parms), pch=15,col="red") 
axis(1,at=c(1,4,7,10,13,16,19,22),labels=X.labels)

legend("topleft", c("Early Data","Late Data","True value"), fill=c("white","grey75","red"),horiz=FALSE,bty="n",cex=1.3)

dev.copy2eps(file="../../c10/figures/CountFitEstimates.eps"); 