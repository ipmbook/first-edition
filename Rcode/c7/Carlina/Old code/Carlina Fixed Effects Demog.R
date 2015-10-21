## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Evolutionary calculations for the Carlina stochastic IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

library(doBy)
library(lme4)
library(parallel)
set.seed(53241986)

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c2",sep="")); 

source("../utilities/Standard Graphical Pars.R");

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book/Rcode/c7/Carlina",sep="")); 


#Read in the Carlina survival data
dataSurv<-read.csv("Carlina survive or died.txt",header=T)

#Convert death to survival and make year a factor
dataSurv$Surv <- dataSurv$die==0
dataSurv$Yeart <- factor(dataSurv$yeart)

#fit some survival models

mod.Surv <- glm(Surv ~ Yeart  , family = binomial, data = dataSurv)
mod.Surv.1 <- glm(Surv ~ Yeart + z  , family = binomial, data = dataSurv)
mod.Surv.2 <- glm(Surv ~ Yeart * z  , family = binomial, data = dataSurv)

anova(mod.Surv,mod.Surv.1,mod.Surv.2,test="Chisq")
AIC(mod.Surv,mod.Surv.1,mod.Surv.2)

#hmmm, there is some evidence of an interaction - we in Rose et al 2002 say there isn't, there is some 
#overdispersion but not much and the interaction is still significant in a quasibinomial model
#Rees and Ellner 2009 assumed both intercept and slope varied between years

#for now let's go with the additive model, refit removing intercept

mod.Surv <- glm(Surv ~ Yeart + z -1 , family = binomial, data = dataSurv)

summary(mod.Surv)

#Read in flowering data

dataFlow<-read.csv("Carlina survive or flower.fre",header=T)

dataFlow$Flow <- dataFlow$flow
dataFlow$Yeart <- factor(dataFlow$yeart)

#fit some flowering models

mod.Flow <- glm(Flow ~ Yeart  , family = binomial, data = dataFlow)
mod.Flow.1 <- glm(Flow ~ Yeart + z  , family = binomial, data = dataFlow)
mod.Flow.2 <- glm(Flow ~ Yeart * z  , family = binomial, data = dataFlow)

anova(mod.Flow,mod.Flow.1,mod.Flow.2,test="Chisq")
AIC(mod.Flow,mod.Flow.1,mod.Flow.2)

#No interaction term - we have data on age at flowering but for now will ignore this and just 
#fit a size dependent model

mod.Flow <- glm(Flow ~ Yeart + z -1 , family = binomial, data = dataFlow)

#Read in flowering data

dataGrow<-read.csv("CRLNAGRW.txt",header=T)

#fit some growth models


with(dataGrow,plot(lst,lst1))

surv.ps <- summaryBy(z + Surv ~ z.classes, data = sim.data.noRepr, na.rm = TRUE)






dataf<-read.csv("Carlina survive or flower.fre",header=T)

#fit model ignoring year effects, used in ESS analysis

fit.flow<-glm(flow~lst,family=binomial,data=dataf)

fit.flow.r<-glmer(flow~lst+(1|yeart),data=dataf,family=binomial)

s.dataf<-dataf[dataf$flow==1,]

size.fl<-with(s.dataf,lst[flow==1])

s.dataf$flow<-0
names(s.dataf)[5]<-"die"

dataf<-read.csv("Carlina survive or died.txt",header=T)

dataf<-rbind(dataf,s.dataf)

minsize<-1.4481; 
maxsize<-5.170528;

