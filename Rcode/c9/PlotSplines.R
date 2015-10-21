rm(list=ls(all=TRUE)); graphics.off(); 
library(fda) 

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c9",sep="")); 

source("../utilities/Standard Graphical Pars.R");

dev.new(width=8,height=4); 
set_graph_pars("panel2"); par(yaxs="i")

############################################################## 
# plot two B-spline bases: d=4, and d=6 with no intercept 
##############################################################
px <- seq(0,1,length=200);
B3 <- create.bspline.basis(rangeval=c(0,1), nbasis=4, norder=4)
py3 <- eval.basis(px,B3);
matplot(px,py3,type="l",lty=1,col="black",xlab="Individual state z",
ylab="Basis function value")

B8 <- create.bspline.basis(rangeval=c(0,1),nbasis=7,norder=3,dropind=1)
py8 <- eval.basis(px,B8);
matpoints(px,py8,type="l",lty=2,col="black");
add_panel_label("a"); 

############################################################## 
# Illustrate function approximation 
##############################################################
px <- seq(0,1,length=25);

## Functions to approximate 
py1 <- exp(-(3*px)^2); # Gaussian
py2 <- px^1.5 # Allometric and NOT a polynomial
py3 <- pmax(0, 3.6*px/(.2+px)-2);
par(xpd=TRUE)
matplot(px,cbind(py1,py2,py3),type="p",pch=1,col="black",
    xlab="Individual state z", ylab="Trait-valued function")

## Approximate with d=6
## Create the basis, and evaluate the functions at values in px  
B <- create.bspline.basis(rangeval=c(0,1), nbasis=6, norder=4)
X <- eval.basis(px,B);

# Spline approx = linear combination of basis functions
# so we can find the least-squares best fit by linear regression
# on the basis function values at the 'data' points
# note: no intercept ("X-1") because spline basis includes an intercept  
cj1 <- coef(lm(py1~X-1)); 
cj2 <- coef(lm(py2~X-1));
cj3 <- coef(lm(py3~X-1));

# Plot the fitted functions at a denser set of points
x <- seq(0,1,length=200)
X2 <- eval.basis(x,B);
y <- X2%*%cbind(cj1,cj2,cj3);
matpoints(x,y,type="l",lty=1,col="black")

# refit 3rd function with a cleverly placed break
B <- create.bspline.basis(rangeval=c(0,1), nbasis=13, dropind=1,norder=2)
X <- eval.basis(px,B)
fit3 <- lm(py3~X-1);
X2 <- eval.basis(x,B);
points(x,X2%*%fit3$coef,type="l",lty=2,lwd=2)
add_panel_label("b"); 
dev.copy2eps(file="../../c6/figures/PlotSplines.eps"); 





