#############################################################################################
# Function to take one time step forward in the monocarp IBM, with initial population state z
# This uses functions in ipm_book/Rcode/c2/Monocarp Demog Funs.R, which must be source'd first. 
#############################################################################################
oneIBMstep<- function(z,m.par.true,constantSeeds=FALSE) {
    ## calculate population size
    pop.size <- length(z)

    ## generate binomial random number for the probability of flowering, where the probability of flowering
    Repr <- rbinom(n=pop.size, prob=p_bz(z, m.par.true), size=1)

    ## number of plants that flowered
    num.Repr <- sum(Repr)

    ## calculate seed production
    Seeds <- rep(NA, pop.size)

    ## Poisson distributed number of seeds
    Seeds[Repr==1] <- rpois(num.Repr, b_z(z[Repr==1], m.par.true))
    
    ## generate the number of recruits
    Recr <- ifelse(num.Repr==0,0,rbinom(1,sum(Seeds, na.rm=TRUE),m.par.true["p.r"])) 

    if(constantSeeds) {
            Seeds[Repr==1] <- b_z(z[Repr==1], m.par.true)    
            Recr <- ifelse(num.Repr==0,0,round(m.par.true["p.r"]*sum(Seeds,na.rm=TRUE))); 
            # cat(Recr," recruits", "\n"); 
    }

    ## generate new recruit sizes
    Rcsz <- rnorm(Recr, mean = m.par.true["rcsz.int"], sd = m.par.true["rcsz.sd"])

    ## for the non-reproductive plants generate random number for survival
    Surv <- rep(NA, pop.size)
    Surv[Repr==0] <- rbinom(n = pop.size - num.Repr, prob = s_z(z[Repr==0], m.par.true), size = 1)
    num.die <- sum(Surv==0, na.rm=TRUE)

    ## index for individuals that did not flower and survived
    i.subset <- which(Repr==0 & Surv==1)

    ## let them grow
    z1 <- rep(NA, pop.size)
    E.z1 <- m.par.true["grow.int"]+m.par.true["grow.z"]*z[i.subset]
    z1[i.subset] <- rnorm(n = pop.size - num.Repr - num.die, mean = E.z1, sd = m.par.true["grow.sd"])
    z1 <- sort(c(Rcsz, z1[i.subset]));
    return(z1)
}    

