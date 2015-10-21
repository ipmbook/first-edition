## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## initial size distribution (assuming everyone is a new recuit)

z <- rnorm(init.pop.size, mean = m.par["rcsz.int"] + m.par["rcsz.z"] * 
             3.2, sd = m.par["rcsz.sd"])

## vectors to store pop size and mean size

pop.size.t <- mean.z.t <- mean.z.repr.t <- numeric(n.yrs)

## Iterate the model using the 'true' parameters and store data in a
## data.frame
yr <- 1
while (yr < n.yrs & length(z) < 5000) {
  
  ## calculate current population size
  pop.size <- length(z)
  
  ## generate binomial random number for survival, where survival depends
  ## on your size z, this is a vector of 0's and 1's, you get a 1 if you
  ## survive
  surv <- rbinom(n = pop.size, prob = s_z(z, m.par), size = 1)
  
  ## generate the size of surviving individuals next year
  i.subset <- which(surv == 1)
  z1 <- rep(NA, pop.size)
  E.z1 <- m.par["grow.int"] + m.par["grow.z"] * z[i.subset]
  z1[i.subset] <- rnorm(n = length(i.subset), mean = E.z1, sd = m.par["grow.sd"])
  
  ## generate a binomial random number for reproduction from surviving
  ## individuals
  
  repr <- rep(NA, pop.size)
  repr[i.subset] <- rbinom(n = length(i.subset), prob = pb_z(z[i.subset], 
                                                             m.par), size = 1)
  
  ## generate a binomial random number for offspring sex (female==1)
  i.subset <- which(surv == 1 & repr == 1)
  osex <- rep(NA, pop.size)
  osex[i.subset] <- rbinom(n = length(i.subset), prob = 1/2, size = 1)
  
  ## generate a binomial random number for offspring recruitment from
  ## surviving / reproducing individuals
  i.subset <- which(surv == 1 & repr == 1 & osex == 1)
  recr <- rep(NA, pop.size)
  recr[i.subset] <- rbinom(n = length(i.subset), prob = pr_z(m.par), 
                           size = 1)
  
  ## generate the size of new recruits
  i.subset <- which(surv == 1 & repr == 1 & osex == 1 & recr == 1)
  z1.rec <- rep(NA, pop.size)
  E.rec.z1 <- m.par["rcsz.int"] + m.par["rcsz.z"] * z[i.subset]
  z1.rec[i.subset] <- rnorm(n = length(i.subset), mean = E.rec.z1, sd = m.par["rcsz.sd"])
  
  ## store the simulation data for the current year, we'll use this later
  sim.data <- data.frame(z, surv, z1, repr, osex, recr, z1.rec)
  
  ## store the population size and mean body mass
  pop.size.t[yr] <- length(z)
  mean.z.t[yr] <- mean(z)
  mean.z.repr.t[yr] <- mean(z[repr == 1], na.rm = TRUE)
  
  ## create new population body size vector
  z <- c(z1.rec[which(surv == 1 & repr == 1 & osex == 1 & recr == 1)], 
         z1[which(surv == 1)])
  
  
  ## iterate the year
  yr <- yr + 1
  
  ## print some info cat(length(z), ' ' ,mean(z), '\n')
  
}

## trim the population size vector to remove the zeros at the end
pop.size.t <- pop.size.t[pop.size.t > 0]

## reassign column names to match the text of chapter 1
names(sim.data) <- c("z", "Surv", "z1", "Repr", "Sex", "Recr", "Rcsz")

