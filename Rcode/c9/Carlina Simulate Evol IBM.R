## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Simulate the IBM 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#store yearly parameters in

store.params.yr <- matrix(NA,nrow=n.yrs,ncol=12)
colnames(store.params.yr) <- names(m.par.true)

# initial population sizes and ages
z <- rnorm(init.pop.size, mean = init.mean.z, sd = init.sd.z)
flow.int.ind <- rnorm(init.pop.size, mean = init.mean.flow.int, sd = init.beta.sd)

## calculate initial pop size and mean size
pop.size.t <- init.pop.size

mean.z.t <- mean(z)
mean.flow.int <- mean (flow.int.ind)
var.flow.int <- var (flow.int.ind)

min.z.t <- min(z)
max.z.t <- max(z)

min.flow.int <- min(flow.int.ind)
max.flow.int <- max(flow.int.ind)


#In the model each plant has a size and it's beta_0, so when we calculate the flowering probability below
#we pass three things your size, your individual beta_0 and the yearly deviation an individual's beta_0
#in year t is then m.par["flow.int"] + flow.ints 

#probability of flowering function depends on individual beta and z
p_bz_ind<- function(z,flow.ints,m.par) {
	linear.p <- m.par["flow.int"] + flow.ints + m.par["flow.z"] * z      # linear predictor
    p <- 1/(1+exp(-linear.p))                                                                     # logistic transformation to probability
    return(p)
}

m.par.sim <- m.par.true

#set mean flow.int <- 0 as each individual has it's own beta_0 and so this is the yearly deviation
m.par.sim["flow.int"] <- 0

## iterate the model using the "true" parameters and store data in a data.frame
yr <- 1

while(yr< n.yrs & length(z) < 1500000) {

	## Simulate yearly parameters
	m.par.year=m.par.sim + qnorm(runif(12,0.001,0.999))*m.par.sd.true
	m.par.year[c("grow.int","rcsz.int")] <- m.par.sim[c("grow.int","rcsz.int")]+matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 

	Recr.year <- sample(nrec,1)*Rec.mult
	
    ## calculate population size
    pop.size <- length(z)
    
    if(pop.size==0) {cat("Extinct","\n")}
    
    ## generate random number for survival - this comes first
    Surv <- rep(NA, pop.size)
    Surv <- rbinom(n = pop.size , prob = s_z(z, m.par.year), size = 1)
    num.die <- sum(Surv==0, na.rm=TRUE)

    ## generate binomial random number for the probability of flowering, where the probability of flowering
    ## depends on your size z, this is a vector of 0's and 1's, you get a 1 if you flower
    Repr <- rep(0, pop.size)
    Repr[Surv==1] <- rbinom(n=pop.size - num.die, prob=p_bz_ind(z[Surv==1], flow.int.ind[Surv==1],m.par.year ), size=1)

    ## number of plants that flowered
    num.Repr <- sum(Repr,na.rm=TRUE)

    ## calculate seed production
    #Seeds <- rep(NA, pop.size)

    ## we'll assume plant make a Poisson distributed number of seeds with a mean given by
    ## exp(params["seed.int"]+params["seed.size"] * z)
    ## rpois generated Poisson distributed random numbers
    Seeds<- rpois(num.Repr, b_z(z[Repr==1], m.par.year))
    Total.seeds <- sum(as.numeric(Seeds),na.rm=TRUE)
    
    p.r <- num.Repr/Total.seeds
    
    #store parameters for the year
    store.params.yr[yr,] <- m.par.year
    #the flow intercept is made up of the yearly deviation and mean for the population
     store.params.yr[yr,"flow.int"] <- store.params.yr[yr,"flow.int"] + mean(flow.int.ind)
    #store.params.yr[yr,"flow.int"] <- mean(flow.int.ind)
    #the probability of establishment is num.Repr/Total.seeds
     store.params.yr[yr,"p.r"] <- p.r
    
    # if(Total.seeds>50000){
    	# Seeds <- round(50000*Seeds/Total.seeds)
    	# Total.seeds <- sum(Seeds,na.rm=TRUE)
    # }
    
    Recr.year <- ifelse(Total.seeds<Recr.year,Total.seeds,Recr.year)
    
    #SPE 
        p.r <- Recr.year/Total.seeds
        store.params.yr[yr,"p.r"] <- p.r
     
     
     if(Total.seeds>0) {
     	
    make.seeds <- Repr
    make.seeds[Repr==1] <- make.seeds[Repr==1]*Seeds>0
    Seeds.greater.0 <- Seeds[Seeds>0]
    Flow.ints.rec <- rep(flow.int.ind[make.seeds==1],Seeds.greater.0)[sample(1:Total.seeds,Recr.year)]
    Flow.ints.rec <- rnorm(Recr.year,Flow.ints.rec,beta.off.sd)
        
    ## generate new recruit sizes
    ## rnorm generated normally distributed random numbers
    Rcsz <- rnorm(Recr.year, mean = m.par.year["rcsz.int"], sd = m.par.year["rcsz.sd"])
    }
    
    ## index for individuals that did not flower and survived
    i.subset <- which(Repr==0 & Surv==1)

    ## let them grow
    
    E.z1 <- m.par.year["grow.int"]+m.par.year["grow.z"]*z[i.subset]
    z1 <- rnorm(n = pop.size - num.Repr - num.die, mean = E.z1, sd = m.par.year["grow.sd"])
	
	if(Total.seeds>0) {
    z1 <- c(Rcsz, z1)
    flow.int.ind <- c(Flow.ints.rec,flow.int.ind[i.subset])
    } else {
    z1 <- z1
    flow.int.ind <- flow.int.ind[i.subset]
    }
    
	z <- z1
    
     min.z.t            <- min(min.z.t,min(z)) 
     max.z.t            <- max(max.z.t,max(z))                         
        	
	mean.flow.int <- c(mean.flow.int,mean(flow.int.ind))
	min.flow.int  <- min(min.flow.int,min(flow.int.ind))
	max.flow.int  <- max(max.flow.int,max(flow.int.ind))
	var.flow.int  <- c(var.flow.int,var(flow.int.ind))
	
	if(yr%%10==0)cat(paste(yr, pop.size, mean.flow.int[yr],"  ",var.flow.int[yr],"\n", sep=" "))

    yr <- yr+1

}














