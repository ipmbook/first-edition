load(file="Soay_data.rda")
library(lme4)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Demographic functions
## 
## Covariates / subsetting variables in demogData
## - obsY:      year of capture
## - NtF:       number of females
## - ageY:      an individual's age (-1 = spring lamb, 0 = summer lamb, 1 = yearling, etc)
## - capWgt:    body mass at capture
## - knowFate:  individuals with known fates == 1 (use only these)
## - isRecMum:  isMum recovered dead (only use isRecMum == 0 when working with spring lambs) 
## - ageMum:    an individual's mother's age in summer prior to year of birth
## - capWgtMum: an individual's mother's capture weight in summer prior to year of birth
## - isTwnMat:  indicator of natal twining envirnment (0 = singleton, 1 = twin)
## - isSurv:    survival indicator (0/1)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## centering of the years
yrcentre <- 1996.5

## /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ ##
##                    survival
## \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ ##

## age == -1
nrow(data <- subset(demogData,
                    knowFate==1 & isRecMum==0 & ageMum>=1 & ageY==-1 &
                    !is.na(capWgt) & capWgt>-1 & !is.na(isSurv) & !is.na(isTwnMat)))

summary(mod.surv.0 <- glmer(isSurv ~ 1 + NtFm1 + I(obsY-yrcentre) + (1|obsY), family=binomial, data=data))

## age >= 0
nrow(data <- subset(demogData, knowFate==1 & ageY>=1 &
                    !is.na(capWgt) & capWgt>1.75 & !is.na(isSurv)))

summary(mod.surv.1 <- glmer(isSurv ~ capWgt + NtF + I(obsY-yrcentre) + (1|obsY), family=binomial, data=data))

## /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ ##
##                    GROWTH
## \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ ##

## age=-1
nrow(data <- subset(demogData, isRecMum==0 & ageMum>=1 & ageY==-1 & !is.na(Ntm1) & capWgt>-1 &
                    !is.na(isTwnMat) & isTwnMat==0 & !is.na(capWgt) & !is.na(capWgtp1) & !is.na(capWgtMum)))
form <- capWgtp1 ~ capWgtMum + NtFm1 + I(obsY-yrcentre) + (1|obsY)
summary(mod.grow.0 <- lmer(form, data=data))

## age>=0
nrow(data <- subset(demogData, ageY>=0 & !is.na(capWgt) & !is.na(capWgtp1)))
form <- capWgtp1 ~ capWgt + NtF + I(obsY-yrcentre) + (1|obsY)
summary(mod.grow.1 <- lmer(form, data=data))

## /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ ##
##                REPRODUCTION
## \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ ##

##
## REPRODUCED SUCCESSFULLY
##

## age >= 0
nrow(data <- subset(demogData, ageY>=0 & !is.na(capWgt)))
form <- isReprop1 ~ capWgt + NtF + I(obsY-yrcentre)+ (1|obsY)
summary(mod.repro <- glmer(form, family=binomial, data=data))

##
## TWINNING RATE
##

## age >= 1
nrow(data <- subset(demogData, ageY>=1 & !is.na(offProd) & !is.na(isTwn)))
form <- isTwnp1 ~ capWgt + NtF + I(obsY-yrcentre) + (1|obsY)
summary(mod.twin <- glmer(form, family=binomial, data=data))