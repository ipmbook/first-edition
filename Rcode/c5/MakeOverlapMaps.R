graphics.off(); 

# data files from Peter B. Adler, not included here
PSSP=read.csv("PSSP/PSSP_genet_xy.csv"); 
ARTR=read.csv("ARTR/ARTR_genet_xy.csv");
POSE=read.csv("POSE/POSE_genet_xy.csv");
HECO=read.csv("HECO/HECO_genet_xy.csv");

speciesList=c("Artemisia tripartita","Hesperostipa comata","Poa secunda","Pseudoroegneria spicata")
myCols=gray(c(0.95,0.8,0.5,0.25)); 
dev.new(width=8,height=8); par(mfrow=c(2,2),mar=c(2,3,4,3),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xpd=FALSE,pty="s"); 
quadrats=c("Q1","Q8","Q13","Q15"); 
years=c(40,50,35,57); 

for (i in 1:4){
    PSSPi=PSSP[as.character(PSSP$quad)==quadrats[i],]
	POSEi=POSE[as.character(POSE$quad)==quadrats[i],]
    ARTRi=ARTR[as.character(ARTR$quad)==quadrats[i],]
	HECOi=HECO[as.character(HECO$quad)==quadrats[i],]
	
    PSSPij=subset(PSSPi,year==years[i]);
    POSEij=subset(POSEi,year==years[i]);
    ARTRij=subset(ARTRi,year==years[i]);
	HECOij=subset(HECOi,year==years[i]);    

# (95,70,50,30)=ARTR,PSSP,POSE,HECO 
    xg=PSSPij$x; yg=PSSPij$y; rg=sqrt(PSSPij$area/pi); 
    if(length(xg)>0) {
        symbols(x=xg,y=yg,circles=rg,xlab=" ",ylab=" ",xlim=c(0,100),ylim=c(0,100),
        inches=FALSE,fg="black",bg="grey25",xaxt="n",yaxt="n");
		title(main=paste(quadrats[i]," ",1900+years[i])); 
    }

    xg=ARTRij$x; yg=ARTRij$y; rg=sqrt(ARTRij$area/pi);
    if(length(xg)>0){
		symbols(x=xg,y=yg,circles=rg,xlab="",ylab="",xlim=c(0,100),ylim=c(0,100),
		fg="black",bg="grey95",inches=FALSE,add=TRUE);
	} 

    xg=POSEij$x; yg=POSEij$y; rg=sqrt(POSEij$area/pi);
    if(length(xg)>0) {
        symbols(x=xg,y=yg,circles=rg,xlab="",ylab="",xlim=c(0,100),ylim=c(0,100),
        inches=FALSE,fg="black",bg="grey50",add=TRUE);
    }
	
    xg=HECOij$x; yg=HECOij$y; rg=sqrt(HECOij$area/pi);
    if(length(xg)>0) {
        symbols(x=xg,y=yg,circles=rg,xlab="",ylab="",xlim=c(0,100),ylim=c(0,100),
        inches=FALSE,fg="black",bg="grey80",add=TRUE);
    }

######## re-add plot of the perimeters only
   xg=POSEij$x; yg=POSEij$y; rg=sqrt(POSEij$area/pi);
    if(length(xg)>0) {
        symbols(x=xg,y=yg,circles=rg,xlab="",ylab="",xlim=c(0,100),ylim=c(0,100),
        inches=FALSE,fg="black",bg="grey50",add=TRUE);
    }


   xg=PSSPij$x; yg=PSSPij$y; rg=sqrt(PSSPij$area/pi); 
   if(length(xg)>0) {
        symbols(x=xg,y=yg,circles=rg,xlab=" ",ylab=" ",xlim=c(0,100),ylim=c(0,100),inches=FALSE,fg="black",bg="grey25",add=TRUE);
   }
	lines(c(0,100,100,0,0),c(0,0,100,100,0))
	
	par(xpd=NA);
	if(i==1) legend(80,0.01,legend=speciesList,fill=c(myCols),cex=1.2,bty="n")
	par(xpd=FALSE); 

	par(new=FALSE); 
  
  
}
dev.copy2eps(file="IdahoMaps.eps");
