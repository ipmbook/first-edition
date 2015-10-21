library(maptools)
library(sp)


speciesList=c("Artemisia tripartita","Hesperostipa comata","Poa secunda","Pseudoroegneria spicata")
myCols=gray(c(0.95,0.8,0.5,0.25)); 

quadrats=c("Q1","Q8","Q13","Q15"); 
years=c(40,50,35,57); 

# Data files from Peter B. Adler, not distributed here
infiles=c("Q1\\Q1_40_C","Q8\\Q8_50_C","Q13\\Q13_35_C","Q15\\Q15_57_C"); 
quadrats=c("Q1","Q8","Q13","Q15"); 

graphics.off(); 
dev.new(width=8,height=8); par(mfrow=c(2,2),mar=c(1,1,3,1),xaxt="n",yaxt="n",xpd=NA,pty="s"); 

for(j in 1:4) {
infile=infiles[j]; 
polyD=readShapePoly(infile)
polyD=polyD[which(is.element(polyD$species,speciesList)),] # subset to species of interest

# set up colors vector
colV=as.character(polyD$species)
for(i in 1:length(speciesList)) colV[colV==speciesList[i]]=myCols[i]

# make figure
plot(polyD, xlim=c(0,1),ylim=c(0,1),col=colV) 
lines(c(0,1,1,0,0),c(0,0,1,1,0))

if(j==1) legend(0.8,0.01,legend=speciesList,fill=c(myCols),cex=1.2,bty="n")

title(main=paste(quadrats[j]," ",1900+years[j])); 

}

dev.copy2eps(file="IdahoShapefileMaps.eps");
savePlot(file="IdahoShapefileMaps.pdf",type="pdf"); 
