## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Simulate the IPM 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

init.pop.size <- 25000
n.yrs <-100
z <- rnorm(init.pop.size, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])
age <- rep(0,init.pop.size)

## Calculate initial pop size and mean size

pop.size.t <- init.pop.size
mean.z.t <- mean(z)
mean.age.t <- mean(age)


## iterate the model using the "true" parameters and store data in a data.frame
yr <- 1
while(yr!= n.yrs & length(z) < 1500000) {

    ## calculate population size
    pop.size <- length(z)

    ## generate binomial random number for the probability of flowering, where the probability of flowering
    ## depends on your size z, this is a vector of 0's and 1's, you get a 1 if you flower
    Repr <- rbinom(n=pop.size, prob=p_bz(z, m.par), size=1)

    ## number of plants that flowered
    num.Repr <- sum(Repr)

    ## calculate seed production
    Seeds <- rep(NA, pop.size)

    ## we'll assume plant make a Poisson distributed number of seeds with a mean given by
    ## exp(params["seed.int"]+params["seed.size"] * z)
    ## rpois generated Poisson distributed random numbers
    Seeds[Repr==1] <- rpois(num.Repr, b_z(z[Repr==1], m.par))

    ## generate the number of recruits
    Recr <- ifelse(num.Repr==0,0,rbinom(1,sum(Seeds, na.rm=TRUE),m.par["p.r"])) 

    ## generate new recruit sizes
    ## rnorm generated normally distributed random numbers
    Rcsz <- rnorm(Recr, mean = m.par["rcsz.int"], sd = m.par["rcsz.sd"])

    ## for the non-reproductive plants generate random number for survival
    Surv <- rep(NA, pop.size)
    Surv[Repr==0] <- rbinom(n = pop.size - num.Repr, prob = s_z(z[Repr==0], m.par), size = 1)
    num.die <- sum(Surv==0, na.rm=TRUE)

    ## index for individuals that did not flower and survived
    i.subset <- which(Repr==0 & Surv==1)

    ## let them grow
    z1 <- rep(NA, pop.size)
    E.z1 <- m.par["grow.int"]+m.par["grow.z"]*z[i.subset]
    z1[i.subset] <- rnorm(n = pop.size - num.Repr - num.die, mean = E.z1, sd = m.par["grow.sd"])
	alive=(Repr==0 & Surv==1)
	
    ## store the simulation data, we'll use this later
   
    if(yr==1) {
        sim.data <- data.frame(z=z, Repr, Seeds, Surv, z1=z1, age=age, alive=alive,yr=1)
    }else{
        sim.data <- rbind(sim.data,data.frame(z=z, Repr, Seeds, Surv, z1=z1, age=age, alive=alive,yr=yr))
    }
    
    ## create new population size vector
    z1 <- c(Rcsz, z1[i.subset])
    age1 <- c(rep(0,Recr),age[i.subset]+1)

    pop.size.t          <- c(pop.size.t, length(z1))
    mean.z.death.t      <- if (yr==1) mean(z[!(Repr==0 & Surv==1)]) else 
                              c(mean.z.death.t,mean(z[!(Repr==0 & Surv==1)]))
    mean.age.death.t    <- if (yr==1) mean(age[!(Repr==0 & Surv==1)]) else 
                              c(mean.age.death.t,mean(age[!(Repr==0 & Surv==1)]))
    mean.fl.z.t         <- if (yr==1) mean(z[Repr==1]) else c(mean.fl.z.t,mean(z[Repr==1])) 
	mean.fl.age.t       <- if (yr==1) mean(age[Repr==1]) else c(mean.fl.age.t,mean(age[Repr==1]))
	
	
    z <- z1
    age <- age1
    
    cat(paste(yr, mean.age.death.t [yr], mean.fl.age.t[yr], "\n", sep=" "))

    yr <- yr+1

}














