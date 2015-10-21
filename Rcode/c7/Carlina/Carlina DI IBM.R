## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Simulate the IBM 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#store yearly parameters in

store.params.yr <- matrix(NA,nrow=n.yrs,ncol=12)
colnames(store.params.yr) <- names(m.par.true)

# initial population sizes and ages
z <- rnorm(init.pop.size, mean = init.mean.z, sd = init.sd.z)

## calculate initial pop size and mean size
pop.size <- init.pop.size
Recr <- rep(0,init.pop.size)

pop.size.t <- rep(NA,n.yrs)
m.par.sim <- m.par.true


## iterate the model using the "true" parameters and store data in a data.frame
yr <- 1

while(yr<= n.yrs & pop.size<2000000 & pop.size>0) {

	## Simulate yearly parameters
	m.par.year=m.par.sim + qnorm(runif(12,0.001,0.999))*m.par.sd.true
	m.par.year[c("grow.int","rcsz.int")] <- m.par.sim[c("grow.int","rcsz.int")]+matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 

	
    ## calculate population size
    pop.size <- length(z)
    
    if(pop.size==0) {cat("Extinct","\n")}
    
    ## generate random number for survival - this comes first
    Surv <- rep(NA, pop.size)
    Surv <- rbinom(n = pop.size , prob = s_z(z, m.par.year), size = 1)
    num.die <- sum(Surv==0, na.rm=TRUE)

    ## generate binomial random number for the probability of flowering, where the probability of flowering
    ## depends on your size z, this is a vector of 0's and 1's, you get a 1 if you flower
    i.subset <- which(Surv == 1)
    Repr <- rep(NA, pop.size)
    Repr[i.subset] <- rbinom(n=pop.size - num.die, prob=p_bz(z[i.subset], m.par.year ), size=1)

    ## number of plants that flowered
    num.Repr <- sum(Repr,na.rm=TRUE)

    ## calculate seed production
    #Seeds <- rep(NA, pop.size)

    ## we'll assume plant make a Poisson distributed number of seeds with a mean given by
    ## exp(params["seed.int"]+params["seed.size"] * z)
    ## rpois generated Poisson distributed random numbers
    i.subset <- which(Surv == 1 & Repr==1)
    Seeds <- rep(NA, pop.size)
    Seeds[i.subset]<- rpois(num.Repr, b_z(z[i.subset], m.par.year))
    Total.seeds <- sum(as.numeric(Seeds),na.rm=TRUE)
    
    #store parameters for the year
    store.params.yr[yr,] <- m.par.year
      
     
     if(Total.seeds>0) {
    
    ##Generate total number fo recruits	
    Recr.year <- rbinom(n=1,prob=m.par.year["p.r"],size=Total.seeds)
        
    ## generate new recruit sizes
    ## rnorm generated normally distributed random numbers
    Rcsz <- rnorm(Recr.year, mean = m.par.year["rcsz.int"], sd = m.par.year["rcsz.sd"])
    }
    
    ## index for individuals that did not flower and survived
    i.subset <- which(Repr==0 & Surv==1)

    ## let them grow
    z1 <- rep(NA, pop.size)
    E.z1 <- m.par.year["grow.int"]+m.par.year["grow.z"]*z[i.subset]
    z1[i.subset] <- rnorm(n = pop.size - num.Repr - num.die, mean = E.z1, sd = m.par.year["grow.sd"])
	
	if(reps==1) {
		
	year.data <- data.frame(z=z, Flow=Repr, Seeds, Surv, z1=z1, Yeart=yr, Recr=Recr)    
	
	if(pop.size>1000)	{
		samp <- sample(1:pop.size,1000)
		year.data <- year.data[samp,]
	}

    if(yr==1) {
        sim.data <- year.data
     } else {
        if(pop.size>0) sim.data <- rbind(sim.data,year.data)
        }
        
    }

	if(Total.seeds>0) {
       z1 <- c(Rcsz, z1[i.subset])
       Recr <- c(rep(1,length(Rcsz)),rep(0,length(z1[i.subset])))
     } else {
     	z1 <- z1[i.subset]
     	Recr <- rep(0,length(z1[i.subset]))
     }	
       
	
    
     if(length(z)>0){
     min.z.t               <- if (yr==1) min(z) else min(min.z.t,min(z)) 
     max.z.t               <- if (yr==1) max(z) else max(max.z.t,max(z))      
     mean.z.t              <- if (yr==1) mean(z) else c(mean.z.t,mean(z))                   
     mean.fl.z.t           <- if (yr==1) mean(exp(z[Repr==1]),na.rm=TRUE) else 
                                                    c(mean.fl.z.t,mean(exp(z[Repr==1]),na.rm=TRUE))
   	 pop.size.t[yr]    <- pop.size
	}
	
	z <- z1
	
	
	if(yr%%2==0)cat(paste(yr, pop.size,"\n", sep=" "))

    yr <- yr+1

}














