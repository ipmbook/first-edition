mean = function(z) 2.751+0.407*z
sd=function (z) sqrt(9.119 *exp(-0.228*mean(z))) #not what's in the paper 
zvals=seq(-4,6,length=100); 
X=cbind(mean(zvals)-1.96*sd(zvals),mean(zvals),mean(zvals)+1.96*sd(zvals))
matplot(zvals,X,lty=c(2,1,2),type="l",col="black",xlim=c(-4,7),
        ylim=c(-4,8)); 
                     
                     
# looks like Fig 1b, but it shouldn't.

#z=runif(300,-2,4); z1=rnorm(300,mean(z),sd(z)); 
#points(z,z1,pch=16,cex=0.5);
# points are much too scattered 
z=runif(300,-3,7); z1=rnorm(300,mean(z),1.35/(1+0.1*z)); 
points(z,z1,pch=16,cex=0.5,col="red");


matplot(zvals,cbind(sd(zvals),sd(0)/(1+0.11*zvals+0.0065*zvals^2)));
