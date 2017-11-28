## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fit and run fixed and mixed effects effects Carlina stochastic IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

set.seed(53241986)

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c2",sep="")); 

source("../utilities/Standard Graphical Pars.R");
source("../utilities/MatrixImage.R");

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c7/Carlina",sep="")); 

load("Yearly parameters.Rdata")

source("Carlina Demog Funs DI.R") 

#####################################################################
#Stochastic perturbation analysis
#####################################################################

nBigMatrix <- 100
n.est <- 20000
n.runin <- 500
minsize <- 1.5
maxsize <- 5
n.years <-20


stoc_pert_analysis<-function(params,n.est,n.runin){
	
	year.i <- sample(1:n.years,n.est+1,replace=TRUE)

	K.year.i <- array(NA,c(n.years,nBigMatrix,nBigMatrix))
	
	for(i in 1:n.years){
		year.K<-mk_K(nBigMatrix,params[,i],minsize,maxsize)
		K.year.i[i,,] <- year.K$K
	}

	h <- year.K$h; 
	
#Calculate mean kernel, v and w

	mean.kernel <- apply(K.year.i,2:3,mean)

	w <- Re(eigen(mean.kernel)$vectors[,1]); 
	v <- Re(eigen(t(mean.kernel))$vectors[,1]);

	# scale eigenvectors <v,w>=1 
	w <- abs(w)/sum(h*abs(w))
	v <- abs(v)
	v <- v/(h*sum(v*w))
    cat(h*sum(v*w)," should = 1","\n")

#Esimate Lambda s
#initialize variables	

	nt<-rep(1/nBigMatrix,nBigMatrix)
	rt.V <- rt.N <- rep(NA,n.est)
	
#Iterate model

	for (year.t in 1:n.est){
		if(year.t%%10000==0) cat("iterate: ", year.t,"\n");

			
		#iterate model with year-specific kernel
		nt1<-K.year.i[year.i[year.t],,] %*% nt
	
		sum.nt1<-sum(nt1)
		
		#Calculate log growth rates  
		
		rt.V[year.t] <- log(sum(nt1*v)/sum(nt*v))
		rt.N[year.t] <- log(sum(nt1)/sum(nt))
		nt <- nt1 / sum.nt1  
	
	}

Ls <- exp(mean(rt.V))
	
	
### Get wt and Rt time series ###
	wt<-matrix(1/nBigMatrix, nrow=n.est+1, ncol=nBigMatrix);
	
	for (i in 1:n.est) {
		
		K             <- K.year.i[year.i[i],,]
		wt[i+1,]  <-K %*% wt[i,]
		wt[i+1,]  <-wt[i+1,]/sum(wt[i+1,]);
		if(i%%10000==0) cat("wt ",i,"\n")

	}

	
### Get vt time series ###
	vt<-matrix(1/nBigMatrix, nrow=n.est+1, ncol=nBigMatrix);
	for (i in (n.est+1):2) {
		
		K           <- K.year.i[year.i[i],,]
		vt[i-1,]  <- vt[i,] %*% K
		vt[i-1,]  <- vt[i-1,]/sum(vt[i-1,]);
		if(i%%10000==0) cat("vt  ",i,"\n")

	}

sens.s <- matrix(0,nrow=nBigMatrix,ncol=nBigMatrix);

elas.s <- sens.s


for (i in n.runin:(n.est-n.runin)) {

		#standard calculations needed for the various formulae
	
		K                 <-   K.year.i[year.i[i],,]
		vt1.wt        <-   vt[i+1,]%*%t(wt[i,])
		vt1.K.wt   <- sum(vt[i+1,] * (K %*% wt[i,]))

		#calculation of the standard sensitivities and elasticities
	
		sens.s<-sens.s+vt1.wt/vt1.K.wt;

		elas.s<-elas.s+K*(vt1.wt/vt1.K.wt);
		
}

 elas.s <- elas.s/(n.est-2*n.runin+1)
 sens.s <- Ls*sens.s/(n.est-2*n.runin+1)
 
 return(list(meshpts=year.K$meshpts,h=h,sens.s=sens.s,elas.s=elas.s,mean.kernel=mean.kernel,Ls=Ls))

}

pert.K <- stoc_pert_analysis(m.par.est,n.est,n.runin)

meshpts <- pert.K$meshpts
h.inv.2 <- 1 / (pert.K$h^2)

Kmean.elas <- pert.K$mean.kernel * pert.K$sens.s / pert.K$Ls

K.sens <- pert.K$sens.s * h.inv.2
K.elas <-  pert.K$elas.s * h.inv.2
K.mean.elas <- Kmean.elas * h.inv.2
K.sd.elas <- K.elas - K.mean.elas

cat("sum elasticities = ",sum(pert.K$elas.s)," should be 1")

## plot these
ikeep <- which(meshpts>1.5 & meshpts<5) # use to extract a region to plot
#postscript("KernSensElas.eps", 
#          width=8, height=8, horizontal=FALSE, paper="special")
## plot the sensitivity kernel and the three elasticity surfaces
set_graph_pars("panel4")
image(meshpts[ikeep], meshpts[ikeep], t(K.sens[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Size (t), z", ylab="Size (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(K.sens[ikeep,ikeep]), add=TRUE)
add_panel_label("a")
image(meshpts[ikeep], meshpts[ikeep], t(K.elas[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Size (t), z", ylab="Size (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(K.elas[ikeep,ikeep]), 
        add=TRUE,levels=c(0.05,0.2,0.4))
add_panel_label("b")
image(meshpts[ikeep], meshpts[ikeep], t(K.mean.elas[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Size (t), z", ylab="Size (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(K.mean.elas[ikeep,ikeep]), 
        add=TRUE,levels=c(0.1,0.4,0.8))
add_panel_label("c")
image(meshpts[ikeep], meshpts[ikeep], t(K.sd.elas[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Size (t), z", ylab="Size (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(K.sd.elas[ikeep,ikeep]), 
        add=TRUE,levels=c(-0.05,-0.2,-0.4))
add_panel_label("d")
# dev.copy2eps(file="~/Repos/ipm_book/c7/figures/CarlinaElasSens.eps")







