setwd("~/Repos/ipm_book/Rcode/c5")
source("../utilities/Standard Graphical Pars.R")
source("Soay DD Demog Funs.R")
L <- 0.5; U <- 3.8; m <- 110; 

#####################################################################
# Iterate density-dependent IPM to find the equilibrium
######################################################################
IPM.dd <- mk_K(Nt=350, m=m, m.par=m.par.true, L=L,U=U)
meshpts <- IPM.dd$meshpts; h<-(U-L)/m;

nbar <- Re(eigen(IPM.dd$K)$vectors[,1]); 
nbar = 350*nbar/sum(h*nbar); 
for(k in 1:50) {
	IPM.dd <- mk_K(Nt=sum(h*nbar), m=m, m.par=m.par.true, L=L,U=U)
	nbar <- IPM.dd$K%*%nbar
}
Nbar <- sum(h*nbar)

### Compute Jacobian matrix at equilibrium and find its dominant eigenvalue
J <- mk_J(nbar=nbar, m=m, m.par=m.par.true, L=L, U=U, eps=0.01);
lam_J <- eigen(J)$values[1]; cat(lam_J,"\n")

######################################################################
# Redefine recruitment probability to depend on 
# steepness parameter rho
######################################################################

pr_z <- function(Nt, m.par) {
    linear.p <- m.par["recr.int"] + m.par["recr.Nt"] * Nt - rho*(Nt-Nbar)
    1/(1+exp(-linear.p))
}

# Sanity check: when rho=0 it should all be the same 
rho <- 0; 
J <- mk_J(nbar=nbar, m=m, m.par=m.par.true, L=L, U=U, eps=0.01);
lam_J <- eigen(J)$values; cat(lam_J[1],"\n") # so far so good 


######################################################################
# Calculate and plot 3 largest eigenvalues as rho increases
# It turns out all three are real 
######################################################################
rhovals <- seq(0,0.1,length=21); 
Revals <- Imvals <- matrix(0,21,5);
Revals[1,] <- Re(lam_J[1:5]); Imvals[1,] <- Im(lam_J[1:5])
for(j in 1:21) {
	rho <- rhovals[j];
	J <- mk_J(nbar=nbar, m=m, m.par=m.par.true, L=L, U=U, eps=0.01);
    lam_J <- eigen(J)$values; cat(j, rhovals[j],lam_J[1], abs(lam_J[1]),"\n") 
    Revals[j,] <- Re(lam_J[1:5]); Imvals[j,] <- Im(lam_J[1:5])
} 
Revals; Imvals; # three dominant eigenvalues are all real 

# plot the three dominant eigenvalues, then add symbols to
# show which one is dominant 

set_graph_pars("panel4")
lam3 <- t(apply(Revals[,1:3],1,sort)); 
matplot(rhovals,lam3,type="l",lty=1,col="black",xlab="Steepness parameter rho",
ylab="Eigenvalues");
abline(h=-1,lty=2)
for(j in 1:21) {
	lamj <- lam3[j,];
	imax <- which(abs(lamj)==max(abs(lamj)))
	points(rhovals[j],lamj[imax])
} 
add_panel_label("a")

rho <- 0.04; 
nreps<-5; N0<-seq(150,400,length=nreps)
Ntvals<-matrix(NA,21,nreps)
for(j in 1:nreps) {
	IPM.dd <- mk_K(Nt=N0[j], m=m, m.par=m.par.true, L=L,U=U)
	nt <- Re(eigen(IPM.dd$K)$vectors[,1]);
	nt <- N0[j]*nt/sum(h*nt); 
	Ntvals[1,j] <- sum(h*nt)
	for(k in 2:21) {
	   IPM.dd <- mk_K(Nt=sum(h*nt), m=m, m.par=m.par.true, L=L,U=U)
	   nt <- IPM.dd$K%*%nt
	   Ntvals[k,j] <- sum(h*nt)
	}
}
matplot(0:20,Ntvals,type="o",pch=1,col="black",lty=1,xlab="Year t",ylab="Total number of females")
text(16,400,"rho=0.04")
add_panel_label("b")

rho <- 0.08; 
Ntvals<-matrix(NA,21,nreps)
for(j in 1:nreps) {
	IPM.dd <- mk_K(Nt=N0[j], m=m, m.par=m.par.true, L=L,U=U)
	nt <- Re(eigen(IPM.dd$K)$vectors[,1]);
	nt <- N0[j]*nt/sum(h*nt); 
	Ntvals[1,j] <- sum(h*nt)
	for(k in 2:21) {
	   IPM.dd <- mk_K(Nt=sum(h*nt), m=m, m.par=m.par.true, L=L,U=U)
	   nt <- IPM.dd$K%*%nt
	   Ntvals[k,j] <- sum(h*nt)
	}
}
matplot(0:20,Ntvals,type="o",pch=1,col="black",lty=1,xlab="Year t",ylab="Total number of females")
text(16,400,"rho=0.08")
add_panel_label("c")

if(FALSE) { # run this only once, after that load the results  
rhovals <- seq(0.03,0.25,length=500); 
Ntvals<-matrix(NA,500,1000);
m <- 50 
for(j in 1:500) {
	rho <- rhovals[j]; 
    IPM.dd <- mk_K(Nt=300, m=m, m.par=m.par.true, L=L,U=U)
    nt <- Re(eigen(IPM.dd$K)$vectors[,1]);
    nt <- 300*nt/sum(h*nt); 
    Ntvals[j,1] <- sum(h*nt)
    for(k in 2:1000) {
	   IPM.dd <- mk_K(Nt=sum(h*nt), m=m, m.par=m.par.true, L=L,U=U)
	   nt <- IPM.dd$K%*%nt
	   Ntvals[j,k] <- sum(h*nt)
	}
	cat(j,"\n")
}
save(list=c("Ntvals","rhovals"),file="Ntvals.Rdata")
} # end if(TRUE/FALSE)

load(file="Ntvals.Rdata");

matplot(rhovals,Ntvals[,901:1000],type="p",pch=".",col="black",
xlab="Steepness parameter rho", ylab="Total population")
add_panel_label("d")

dev.copy2eps(file="../../c5/figures/SoayBif.eps")

graphics.off(); 
dev.new(width=8,height=4)
set_graph_pars("panel2")
rho=0.2; m <- 110;  
Nt<-numeric(1000); ntvals<-matrix(NA,1000,110); 

IPM.dd <- mk_K(Nt=300, m=m, m.par=m.par.true, L=L,U=U)
nt <- Re(eigen(IPM.dd$K)$vectors[,1])
nt <- 300*nt/sum(h*nt) 
Nt[1] <- sum(h*nt)
ntvals[1,] <- nt; 
for(k in 2:1000) {
	IPM.dd <- mk_K(Nt=sum(h*nt), m=m, m.par=m.par.true, L=L,U=U)
	nt <- IPM.dd$K%*%nt
	Nt[k] <- sum(h*nt)
	ntvals[k,] <- nt; 
	cat(k,"\n")
}
plot(Nt[500:999],Nt[501:1000],cex=0.5,xlab="Total #females N(t)", ylab="Total #females N(t+1)")
add_panel_label("a")
Nt1<-Nt[501:1000]; nt1<-ntvals[501:1000,]

imin <- which(Nt1 < min(Nt1) + 1)
imax <- which(Nt1 > max(Nt1)-1)
imean <- which( abs(Nt1-mean(Nt1))<1)

for(j in 1:500) nt1[j,] <- nt1[j,]/sum(h*nt1[j,]) # scale to N=1
meshpts <- IPM.dd$meshpts
par(yaxs="i")
matplot(meshpts,t(nt1[imin,]), type="l",xlim=c(1.8,U),col="black",xlab="Size z", ylab="Relative frequency")
matpoints(meshpts,t(nt1[imax,]), type="l",lty=2, xlim=c(1.8,U),col="red")
matpoints(meshpts,t(nt1[imean,]), type="l",lty=3,xlim=c(1.8,U),col="blue")
legend("topleft",bty="n",legend=c("min(17)","max(18)","mean(9)"),lty=c(1,2,3),col=c("black",
			"red","blue"),lwd=2)
add_panel_label("b")
dev.copy2eps(file="../../c5/figures/Soay1D.eps")
