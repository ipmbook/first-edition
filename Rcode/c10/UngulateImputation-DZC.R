## working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/Repos"); 
setwd(paste(root,"/ipm_book",sep="")); 

## run the utility functions
source("./Rcode/utilities/Standard Graphical Pars.R")
## run the ungulate IBM
source("./Rcode/c2/Ungulate Demog Funs.R")
## add an extra parameter for the size-dependent growth s.d. 
m.par.true["grow.sd.c"] <- -0.07
m.par.true["grow.sd.d"] <- +0.05


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Section 1: Functions to generate the "data" and do the imputation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

############################################################# 
## Growth kernel: assumed to be already fitted 
#############################################################
g_z1z <- function(z1, z, m.par){
  E.z1 <- m.par["grow.int"] + m.par["grow.z"] * z  
  sd.z1 <- m.par["grow.sd.c"] + m.par["grow.sd.d"] * z                   
  p.den.grow <- dnorm(z1, mean = E.z1, sd = sd.z1) 
  return(p.den.grow)
}
###############################################################
## Function to simulate individual-based population model 
###############################################################
sim_growth <- function(m.par, cohortSize, maxA) {
  ## id vector
  id  <- factor(seq_len(cohortSize))
  ## initial size distribution
  z <- rnorm(cohortSize, 
             mean = m.par["rcsz.int"] +  m.par["rcsz.z"] * 3.2, sd = m.par["rcsz.sd"])
  ## set up the storage object
  cohortData <- list()
  cohortData[[1]] <- data.frame(id, a=0, z)
  ## iterate the cohort
  a <- 0
  repeat { 
    ## generate binomial random number for survival
    surv <- rbinom(n=cohortSize, prob=s_z(z, m.par), size=1)
    if (a==maxA) surv <- rep(0, length(surv)) # max age is 15
    ## index to select survivors
    i.subset <- which(surv == 1)
    ## update id and size vectors
    id <- id[i.subset]; z <- z[i.subset] 
    ## update cohort size
    cohortSize <- length(i.subset)
    ## generate new sizes
    if (cohortSize > 0) {
      ## generate the size of surviving individuals next age
      E.z1 <- m.par["grow.int"] + m.par["grow.z"] * z
      sd.z1 <- m.par["grow.sd.c"] + m.par["grow.sd.d"] * z
      z <- rnorm(n = cohortSize, mean = E.z1, sd = sd.z1) 
      ## store the cohort data for the current age
      a <- a + 1
      cohortData[[a+1]] <- data.frame(id, z, a = a)
    } else break
  }
  cohortData <- do.call("rbind", cohortData)
  cohortData <- transform(cohortData, zorig = z, zimpt = NA)
  return(cohortData)
}

#############################################################################
# Function to remove leading and trailing NAs in data on each individual
#############################################################################
prune_data <- function(cohortData) {
  idSet <- as.character(unique(cohortData$id))
  i <- 0; prunedData <- list()
  for (idNow in idSet) {
    indivData <- subset(cohortData, id == idNow)
    isMissing <- any(is.na(indivData$z))
    if (isMissing) {
      runLengthsNA <- rle(is.na(indivData$z))
      if (sum(!runLengthsNA$values)>=2) {
        # 1 - remove NAs at the 'head' and 'tail'
        headIsNA <- head(runLengthsNA$val, 1); headLen <- head(runLengthsNA$len, 1)
        tailIsNA <- tail(runLengthsNA$val, 1); tailLen <- tail(runLengthsNA$len, 1)
        if (headIsNA) {
          indivData <- tail(indivData, -headLen)
        } 
        if (tailIsNA) {
          indivData <- head(indivData, -tailLen)
        }
        # 2 - break up the seperate NA runs
        runLengthsNA <- rle(is.na(indivData$z))
        posRuns <- cumsum(runLengthsNA$len)
        i1 <- head(posRuns[!runLengthsNA$values], -1)
        i2 <- posRuns[ runLengthsNA$values]+1L
        for (iRun in seq_along(i1)) {
          i <- i + 1
          prunedData[[i]] <- indivData[seq.int(i1[iRun], i2[iRun]),]
        }
      }
    }  
  }
  return(prunedData)
}

#########################################################################
# Function to run the imputation scheme
#########################################################################
impute_missing <- function(imputeData, m.par, L, U, m) {
  # usual midpoint rule calcs
  h <- (U - L) / m
  meshpts <- L + ((1:m) - 1/2) * h
  
  # lengths of the (size, NA, NA, ..., NA, size) blocks to be filled in
  maxT <- sapply(imputeData, nrow)
  
  # make the set of G^T matrices
  G <- h * outer(meshpts, meshpts, g_z1z, m.par)
  for(i in 1:m) G[,i] <- G[,i]/sum(G[,i]); # scale to set eviction=0  
  G.T <- list()
  G.T[[1]] <- G; for (i in seq.int(2, max(maxT))) G.T[[i]] <- G %*% G.T[[i-1]]
  
  # pre- and post-NA values for each NA block 
  b <- sapply(imputeData, function(df) head(df$z, 1))  
  B <- sapply(imputeData, function(df) tail(df$z, 1))
  
  # Apply the methodology in Box 10.4
  for (iGroup in seq_along(imputeData)) {
    # find the meshpoint indices of observed sizes at start and end of the block  
    kappa.b <- suppressWarnings(which.min(abs(meshpts-b[iGroup])))
    kappa.B <- suppressWarnings(which.min(abs(meshpts-B[iGroup])))
    
    # Fill in the NAs, which are in rows 2 to maxT-1 of the block 
    for (iT in seq.int(2, maxT[iGroup]-1)) { 
      i <- iT-1;                      # row of most recent observed or imputed size;   
      dT <- maxT[iGroup]-i;           # time to next observed size; T in the book.         
      p1 <- G.T[[dT-1]] [kappa.B, ] * G.T[[1]] [, kappa.b] / G.T[[dT]] [kappa.B, kappa.b]
      kappa.b <- sample(length(meshpts), 1, prob=p1) # update the prior 'observed' size index
      imputeData[[iGroup]][iT, "zimpt"] <- meshpts[kappa.b]
    }  

  }
  return(do.call("rbind", imputeData)) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Section 2: Generate the "data" and do the imputation. 
#   First, with data missing at random, as the method assumes; then with
#   P(missing) as a function of size (decreasing or increasing)  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set.seed(270815)

### missing is random wrt to size
cohortData1 <- cohortData <- sim_growth(m.par.true, cohortSize = 3200, maxA = 15)
ii <- sample(nrow(cohortData), 3200, replace=FALSE)
cohortData$z[ii] <- NA
imputeData <- prune_data(cohortData)
allDataMCAR <- impute_missing(imputeData, m.par.true, L = 1.7, U = 4.2, m = 500)

### p(missing) decreases with size
cohortData <- sim_growth(m.par.true, cohortSize = 5000, maxA = 15)
pmiss <- 1/(1+exp(-(-12 * (cohortData$z-mean(cohortData$z)))))
ii <- sample(nrow(cohortData), 6000, replace=FALSE, prob=pmiss)
cohortData$z[ii] <- NA
imputeData <- prune_data(cohortData)
allDataNMAR1 <- impute_missing(imputeData, m.par.true, L = 1.7, U = 4.2, m = 500)

# p(missing) increases with size
cohortData <- sim_growth(m.par.true, cohortSize = 2600, maxA = 15)
pmiss <- 1/(1+exp(-(+12 * (cohortData$z-mean(cohortData$z)))))
ii <- sample(nrow(cohortData), 2600, replace=FALSE, prob=pmiss)
cohortData$z[ii] <- NA
imputeData <- prune_data(cohortData)
allDataNMAR2 <- impute_missing(imputeData, m.par.true, L = 1.7, U = 4.2, m = 500)

# check we have approximately the same number of cases
c(sum(!is.na(allDataMCAR$zimpt)), 
  sum(!is.na(allDataNMAR1$zimpt)),
  sum(!is.na(allDataNMAR2$zimpt)))

graphics.off(); set_graph_pars("panel4"); 
z <- sample(cohortData1$z,500,replace=FALSE); 
E.z1 <- m.par.true["grow.int"] + m.par.true["grow.z"] * z  
sd.z1 <- m.par.true["grow.sd.c"] + m.par.true["grow.sd.d"] * z 
z1 <- rnorm(length(z),mean=E.z1,sd=sd.z1)
plot(z,z1,xlab="Initial size",ylab="Final size",cex=0.75);
zfit <- lm(z1~z); abline(zfit); add_panel_label("a"); 
 
plot(zimpt~zorig, data=allDataMCAR, 
     xlab="True Size", ylab="Imputed Size",
     cex=0.7, pch=20, xlim=c(2.3, 3.8), ylim=c(2.3, 3.8), col=rgb(1,0,0, alpha=0.3))
abline(0, 1, lty=2); add_panel_label("b"); 
# fitimp <- lm(zimpt~zorig,data=allDataMCAR); abline(fitimp,col="red");  

# with(allDataNMAR2, hist(zimpt-zorig))

plot(zimpt,0~zorig, data=allDataNMAR1, 
     xlab="True Size", ylab="Imputed Size",
     cex=0.7, pch=20, xlim=c(2.3, 3.8), ylim=c(2.3, 3.8), col=rgb(0,1,0, alpha=0.3))
abline(0, 1, lty=2)
fitimp <- lm(zimpt~zorig,data=allDataNMAR1); abline(fitimp,col="red"); 

plot(zimpt~zorig, data=allDataNMAR2, 
     xlab="True Size", ylab="Imputed Size",
     cex=0.7, pch=20, xlim=c(2.3, 3.8), ylim=c(2.3, 3.8), col=rgb(0,0,1, alpha=0.3))
abline(0, 1, lty=2)
fitimp <- lm(zimpt~zorig,data=allDataNMAR2); abline(fitimp,col="red"); 

dev.new(); par(mfrow=c(2,2)); 
# bit clearer with density plots
plot(with(allDataMCAR,   density(na.exclude(zimpt-zorig), bw=0.05)),
     col="black", main="")
lines(with(allDataNMAR1,  density(na.exclude(zimpt-zorig), bw=0.05)),
      col="blue")
lines(with(allDataNMAR2,  density(na.exclude(zimpt-zorig), bw=0.05)),
      col="red")
abline(v=0, lty=2)






