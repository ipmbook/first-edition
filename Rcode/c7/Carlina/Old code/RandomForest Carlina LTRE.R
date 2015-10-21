rm(list=ls(all=TRUE))
require(randomForest)

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c2",sep="")); 
source("../utilities/Standard Graphical Pars.R");
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c7/Carlina",sep="")); 

load("Big params plus lambda.Rdata")

# Format it as a data frame with only the parameters that vary over time. 
X <- data.frame(t(params.plus.lambda)); 
names(X) <- row.names(params.plus.lambda)
vars<- which(apply(X,2,sd)>0);
X <- X[,vars];  

# 'Tune' the random forest parameter mtry
out<- tuneRF(x=subset(X,select=-lambda.i), y=X$lambda.i,ntreeTry=500,stepfactor=1,improve=0.02); 
plot(out); ## tells us to use mtry=4 

# Fit the model, computing importance measures
x=subset(X,select=-lambda.i)
lambdaRF <- randomForest(x=x,y=X$lambda.i,ntree=500,importance=TRUE,mtry=4)

graphics.off(); 
dev.new(width=8,height=5); set_graph_pars("panel2");
# verify that the number of trees is enough
# the plot of R^2 vs. number of trees used should saturate 
plot(lambdaRF$rsq,type="l",xlab="Number of trees",ylab="Model r^2"); add_panel_label("a")

# Plot the importance results 
varImpPlot(lambdaRF,type=1,scale=FALSE,main="")
add_panel_label("b")

### Alternative importance measure. This takes some patience, to fit all the RF models
nvars <- ncol(x); 
mseDrop1 <- numeric(nvars);
for(j in 1:nvars) {
	RFi <- randomForest(x=x[,-j],y=X$lambda.i,keep.forest=FALSE,ntree=500,mtry=4)
	mseDrop1[j] <- min(RFi$mse)
    cat(j,"\n")	
} 
mseFull <- min(lambdaRF$mse);
deltaMSE <- (mseDrop1 - mseFull)/var(X$lambda.i); 

e <- order(deltaMSE); 
dsort <- deltaMSE[e];
names(dsort) <- names(x)[e];
dev.new(); dotchart(dsort,xlab="Relative increase in MSE",xlim=c(0,max(dsort)))
