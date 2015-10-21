## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Simulate the IPM 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# initial population sizes and ages
z <- rnorm(init.pop.size, mean = 2.9, sd = m.par.true["rcsz.sd"])
flow.int.ind <- rnorm(init.pop.size, mean = init.mean.flow.int, sd = init.beta.sd)

## calculate initial pop size and mean size

pop.size.t <- init.pop.size
mean.z.t <- mean(z)
mean.flow.int <- mean (flow.int.ind)

#probability of flowering function depends on individual beta and z
p_bz_ind<- function(z,flow.ints,m.par) {
	linear.p <- flow.ints + m.par["flow.z"] * z      # linear predictor
    p <- 1/(1+exp(-linear.p))                        # logistic transformation to probability
    return(p)
}

## iterate the model using the "true" parameters and store data in a data.frame
yr <- 1
while(yr<= n.yrs & length(z) < 1500000) {

    ## calculate population size
    pop.size <- length(z)

    ## generate binomial random number for the probability of flowering, where the probability of flowering
    ## depends on your size z, this is a vector of 0's and 1's, you get a 1 if you flower
    Repr <- rbinom(n=pop.size, prob=p_bz_ind(z, flow.int.ind,m.par.true ), size=1)

    ## number of plants that flowered
    num.Repr <- sum(Repr)

    ## calculate seed production
    #Seeds <- rep(NA, pop.size)

    ## we'll assume plant make a Poisson distributed number of seeds with a mean given by
    ## exp(params["seed.int"]+params["seed.size"] * z)
    ## rpois generated Poisson distributed random numbers
    Seeds<- rpois(num.Repr, b_z(z[Repr==1], m.par.true))
    
    Total.seeds <- sum(Seeds,na.rm=TRUE)
    
    Flow.ints.rec <- rep(flow.int.ind[Repr==1],Seeds)[sample(1:Total.seeds,Recr)]
    Flow.ints.rec <- rnorm(Recr,Flow.ints.rec,beta.off.sd)
    ## generate the number of recruits
    
    
    ## generate new recruit sizes
    ## rnorm generated normally distributed random numbers
    Rcsz <- rnorm(Recr, mean = m.par.true["rcsz.int"], sd = m.par.true["rcsz.sd"])

    ## for the non-reproductive plants generate random number for survival
    Surv <- rep(NA, pop.size)
    Surv[Repr==0] <- rbinom(n = pop.size - num.Repr, prob = s_z(z[Repr==0], m.par.true), size = 1)
    num.die <- sum(Surv==0, na.rm=TRUE)

    ## index for individuals that did not flower and survived
    i.subset <- which(Repr==0 & Surv==1)

    ## let them grow
    
    E.z1 <- m.par.true["grow.int"]+m.par.true["grow.z"]*z[i.subset]
    z1 <- rnorm(n = pop.size - num.Repr - num.die, mean = E.z1, sd = m.par.true["grow.sd"])
	
    z1 <- c(Rcsz, z1)
    flow.int.ind <- c(Flow.ints.rec,flow.int.ind[i.subset])
    
	z <- z1
    
     min.z.t            <- if (yr==1) min(z) else min(min.z.t,min(z)) 
     max.z.t            <- if (yr==1) max(z) else max(max.z.t,max(z))                         
     mean.fl.z.t         <- if (yr==1) mean(exp(z[Repr==1])) else c(mean.fl.z.t,mean(exp(z[Repr==1])))
   	
	mean.flow.int <- if (yr==1) mean(flow.int.ind) else c(mean.flow.int,mean(flow.int.ind))
	min.flow.int <- if (yr==1) min(flow.int.ind) else min(min.flow.int,min(flow.int.ind))
	max.flow.int <- if (yr==1) max(flow.int.ind) else max(max.flow.int,max(flow.int.ind))
	var.flow.int <- if (yr==1) var(flow.int.ind) else c(var.flow.int,var(flow.int.ind))
	
	if(yr%%10==1)cat(paste(yr, mean.flow.int[yr],"  ",var.flow.int[yr],"\n", sep=" "))

    yr <- yr+1

}














