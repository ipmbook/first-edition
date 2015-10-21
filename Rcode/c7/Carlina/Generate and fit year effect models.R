#Fitting bivariate models

rm(list=ls(all=TRUE))

require(lattice)
require(lme4)
require(nlme)
require(MASS)
library(MCMCglmm)
library(rjags)
set.seed(53241986)

#recruit parameters
a.rec  <- 0.6
b.rec <- 0.02

sd.rec <- 0.1


#growthparameters	 
a.g         <- 1
b.g         <-0.6

sd.grow <- 0.2


#Var-Covar matrix for yearly variation in intercepts 
Sigma <- matrix(c(0.3,0.4,0.4,0.6),2,2)
cov2cor(Sigma)

#Simulation parameters
n.samp    <- 10000
n.years   <- 100
year.codes       <- factor(paste("year",1:n.years,sep="."))
     
years <- gl(n.years , n.samp/n.years, labels=year.codes)

#year effects on growth and recruitment, here they are independent but they can be correlated using mvrnorm in MASS

year.intercepts <- mvrnorm(n.years,c(a.g,a.rec),Sigma)

year.effects.g <- year.intercepts[,1]
year.effects.r <- year.intercepts[,2]

year.g <- year.effects.g[match(years,year.codes)]

year.r <- year.effects.r[match(years,year.codes)]

size <- rnorm(n.samp,1,0.3)

size1 <-  year.g  + b.g * size  + rnorm(n.samp,0,sd.grow)

plot(size1~size)

rec <- year.r + b.rec * size  +  rnorm(n.samp,0,sd.rec)

#do univariate fits to check the data are what we expect.
#In a real example we would use these to inform the construction of the multivariate model

tmp.g=groupedData(size1~size|year.g)

intervals(lme(size1 ~ size,random=~ 1|year.g, data=tmp.g))

tmp.r=groupedData(rec ~ 1|year.g)

intervals(lme(rec ~ size,random=~ 1|year.g, data=tmp.r))

#well that works OK

both <- data.frame(size1=c(size1,rec),size=c(size,size),demo.p=gl(2,n.samp),year=c(years,years))

both.g <-groupedData(size1 ~ size|year,data=both)

#fit.trait.1=lme(size1~demo.p + size -1 ,random = ~ (demo.p-1)|year , data=both.g,na.action=na.omit)

#allows the year effects to be correlated

fit.trait.1=lme(size1~demo.p / size -1 ,random = ~ (demo.p-1)|year , weights=varIdent(form = ~ 1|demo.p) , data=both.g)

summary(fit.trait.1)

#assumes the year effects are independent

fit.trait.2=lme(size1~demo.p / size -1 ,random=list(year=pdDiag(~ demo.p-1)) , weights=varIdent(form = ~ 1|demo.p) , data=both.g,na.action=na.omit)

summary(fit.trait.2)

anova(fit.trait.1,fit.trait.2)

getVarCov(fit.trait.1)

cat("grow.sd =",sd.grow," estimate ",summary(fit.trait.2)$sigma)

est.rec.sd <- coef(fit.trait.2$modelStruct$varStruct, uncons = FALSE) * (summary(fit.trait.2)$sigma)

cat("rec.sd =",sd.rec," estimate ",est.rec.sd )


#all the parameters look OK

random.E=ranef(fit.trait.2)

random.E <- random.E[order(as.numeric(row.names(random.E))),]

par(mfrow=c(1,3),pty="s",bty="l",pch=19)

plot(random.E[,1],year.effects.g)

plot(random.E[,2],year.effects.r)

plot(random.E[,1],random.E[,2])

#seems to get the year effects OK

############################################
#MCMCglmm version
mcmc.data <- data.frame(size1=size1, size=size, rec=rec, year=years)

#assumes the year effects are correlated

fit.trait.1.mcmc <- MCMCglmm(cbind(size1,rec)~trait / size -1 , family = c("gaussian","gaussian"), random = ~ us(trait):year, rcov=~idh(trait):units, data=mcmc.data)

summary(fit.trait.1.mcmc)
plot(fit.trait.1.mcmc)

#assumes the year effects are independent

fit.trait.2.mcmc <- MCMCglmm(cbind(size1,rec)~trait / size -1 , family = c("gaussian","gaussian"), random = ~ idh(trait):year, rcov=~idh(trait):units, data=mcmc.data)

summary(fit.trait.2.mcmc)

#allow individual level correlation

fit.trait.3.mcmc <- MCMCglmm(cbind(size1,rec)~trait / size -1 , family = c("gaussian","gaussian"), random = ~ us(trait):year, rcov=~us(trait):units, data=mcmc.data)

summary(fit.trait.3.mcmc)
plot(fit.trait.3.mcmc)







