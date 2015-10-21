setwd("~/Repos/ipm_book/Rcode/c5")
source("../utilities/Standard Graphical Pars.R")



spp="PSSP"
dcritG = seq(5,15,length=30);
alphaG=1/dcritG^2

# alphaG=seq(0.005,0.03,length=20)


# survival --------------------------------------------------
survD=read.csv("PSSP_survD.csv")
distD=read.csv("PSSP_genet_xy.csv")
sD=subset(survD,allEdge==0)
sD$year=as.factor(sD$year)
sD$logarea=log(sD$area)

# calculate crowding 
sD$quad=as.character(sD$quad)
distD=distD[,c("quad","year","trackID","area","x","y")]
survGW=matrix(NA,dim(sD)[1],length(alphaG))
for(i in 1:dim(sD)[1]){
if (i%%100==1) cat(i,"\n"); 
 tmpD=subset(distD,year==sD$year[i] & quad==sD$quad[i])
 focal=which(tmpD$trackID==sD$trackID[i])
 xx=tmpD$x[focal] ; yy=tmpD$y[focal]
 tmpD$distance=sqrt((xx-tmpD$x)^2+(yy-tmpD$y)^2)
 tmpD=subset(tmpD,distance>0)
 for(j in 1:length(alphaG)){
    if(dim(tmpD)[1]>0){
      survGW[i,j]=sum(exp(-1*alphaG[j]*tmpD$distance^2)*tmpD$area)
     }else{
      survGW[i,j]=0
     }
   }   
} # next record

# fit with Gaussian form of alpha
survAICg=matrix(NA,length(alphaG),5)
survlogLikg=matrix(NA,length(alphaG),5)
for(j in 1:length(alphaG)){
   cat(j,"\n"); 
   sD$crowd=survGW[,j]
   out=glm(survives~Group+logarea+crowd,data=sD,family=binomial)
   survAICg[j,1]=AIC(out)
   survlogLikg[j,1]=logLik(out)[1]
   out=glm(survives~Group+logarea+crowd+year,data=sD,family=binomial)
   survAICg[j,2]=AIC(out)
   survlogLikg[j,2]=logLik(out)[1]
   out=glm(survives~Group+logarea+crowd+year+logarea:year,data=sD,family=binomial)
   survAICg[j,3]=AIC(out)
   survlogLikg[j,3]=logLik(out)[1]
   out=glm(survives~Group+logarea+crowd+year+logarea:year+crowd:year,data=sD,family=binomial)
   survAICg[j,4]=AIC(out)
   survlogLikg[j,4]=logLik(out)[1]
   out=glm(survives~Group+logarea*crowd+year+logarea:year,data=sD,family=binomial)
   survAICg[j,5]=AIC(out)
   survlogLikg[j,5]=logLik(out)[1]
}

source("../utilities/Standard Graphical Pars.R")

add_panel_label <- function(ltype="a",cex=1)
{
    text <- paste(LETTERS[letters==ltype], ")", sep="")
    mtext(text=text, side=3, adj=0, cex=cex)
}

set_graph_pars("panel1")
par(cex.axis=1.35,cex.lab=1.65,cex.main=1.35,mgp=c(2.25,1,0),bty="l")
matplot(dcritG,survAICg[,-1],type="l",xlab="Competition radius (cm)",ylab="AIC",
lty=1:4,col="black",main="Pseudoroegneria",lwd=c(1,2,2,2))
add_panel_label("b",cex=2); 
dev.copy2eps(file="../../c5/figures/PSSPSurvivalRadius.eps")

for(j in 2:5) {
	out=which(survAICg[,j]==min(survAICg[,j])); 
	cat(1/sqrt(alphaG[out]),survAICg[out,j],"\n");
}	