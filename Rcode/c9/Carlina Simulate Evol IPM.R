###################################################################################################
#  Function to iterate Carlina stochastic IPM with size and flowering-intercept variation 
#  This simulates one stochastic realisation of the IPM from each starting point 
#  The environment process is not one of those used in the IBM simulations 
###################################################################################################
iterate_model <- function(params,n.iter,init.beta.mean,init.beta.sd,meshpts.beta,meshpts.z) {
    
    nBigMatrix.z <- length(meshpts.z); 
    nBigMatrix.beta <- length(meshpts.beta); 
    
    nt <- matrix(NA,nrow=nBigMatrix.z,ncol=nBigMatrix.beta)
	
	params.yr <- matrix(NA,ncol=length(params),nrow=n.iter)
    params.yr <- data.frame(params.yr); 
    names(params.yr) <- names(params); 
	
	z.s <- dnorm(meshpts.z, mean = init.mean.z, sd = init.sd.z)
    flow.int.s <- dnorm(meshpts.beta, mean = init.beta.mean, sd = init.beta.sd)
    nt <- init.pop.size*matrix(outer(z.s,flow.int.s),nrow=nBigMatrix.z)
    
    	
	## Matrix to store distribution of flowering intercepts 
	nBeta <- length(meshpts.beta);
	betaDist <- matrix(0, nBeta,n.iter)
	betaDist[,1] = apply(nt,2,sum); betaDist[,1]=betaDist[,1]/sum(betaDist[,1])
	
	#array of P matrices for each beta
    P.beta <- array(NA,c(nBigMatrix.beta,nBigMatrix.z,nBigMatrix.z))

	#array to store 
    Seeds.beta <- array(NA,c(nBigMatrix.beta,nBigMatrix.z))

	#Number of seeds by a size z plant
    Seeds_z <- function ( z, m.par) {
            return( p_bz(z, m.par) * b_z(z, m.par) )
     }

	#Inheritance kernel
	M <- matrix(NA,ncol=nBigMatrix.beta,nrow=nBigMatrix.beta)
	for(i in 1:nBigMatrix.beta){
 		 M[,i] <- dnorm(meshpts.beta,mean=meshpts.beta[i],sd=beta.off.sd)
  	}

	#pdf of offspring sizes
	beta2 <- meshpts.beta*meshpts.beta
	mean.beta <- rep(NA,n.iter)
	mean.beta[1] <- init.beta.mean
	var.beta <- rep(NA,n.iter)
	var.beta[1] <- init.beta.sd*init.beta.sd
	prop.Recr <- rep(NA,n.iter)
	mean.beta2 <- rep(NA,n.iter)
	mean.beta2[1] <- sum(flow.int.s*beta2)/sum(flow.int.s)
	var.beta2 <- rep(NA,n.iter);
	
	########### Compute var(beta2) using formula Var=E[(X-Xbar)^2], not Var=E[X^2]-E[X]^2. 
	beta2Dev <- beta2-mean.beta2[1]; 
	var.beta2[1] <- sum(flow.int.s*beta2Dev^2)/sum(flow.int.s)
	
	index <- 1:nBigMatrix.beta
	
	#Iterate the model
	for(gen in 2:n.iter){
		
	#generate yearly parameters
	m.par.year=m.par.sim + qnorm(runif(12,0.001,0.999))*m.par.sd.true
	m.par.year[c("grow.int","rcsz.int")] <- m.par.sim[c("grow.int","rcsz.int")]+matrix(qnorm(runif(2,0.001,0.999)),1,2) %*% chol.VarCovar.grow.rec 
	Recr.year <- sample(nrec,1)*200
	
    # build kernel functions
	params <- m.par.year; 
	for(i in 1:nBigMatrix.beta){
		params["flow.int"] <- m.par.year["flow.int"] + meshpts.beta[i]
		P.beta[i,1:nBigMatrix.z,1:nBigMatrix.z] <- h.z * (outer(meshpts.z, meshpts.z, P_z1z, m.par = 	params))
		Seeds.beta[i,1:nBigMatrix.z] <- Seeds_z(meshpts.z, params)
	}

	off.pdf <- c_0z1(meshpts.z,m.par.year)
		
	#calculate seeds produced by beta_i plants
	seeds.from.betai <- sapply(index,function(i) h.z*sum(Seeds.beta[i,] %*% nt[,i]))

	#redistribute seeds according to inheritance kernel
	seeds.with.betai <- h.beta * (M %*% seeds.from.betai)

	f.recruits.with.betai <- seeds.with.betai/sum(seeds.from.betai)

	nt1 <- sapply(index,function(i) P.beta[i,,] %*% nt[,i] + Recr.year * off.pdf * f.recruits.with.betai[i])

	n.surv <- sapply(index,function(i) P.beta[i,,] %*% nt[,i])
	n.surv <- sum(n.surv) * h.z * h.beta
	nt <- nt1
	
	betaDist[,gen] = apply(nt,2,sum); betaDist[,gen] <- betaDist[,gen]/sum(betaDist[,gen]); 
	betaFreq <- betaDist[,gen] 
	
	params.yr[gen-1,] <- m.par.year; 
	params.yr[gen-1,"p.r"] <- Recr.year/sum(h.beta*seeds.with.betai);  
	
	#### Compute variances of beta and beta^2 using formula Var=E[(X-Xbar)^2] 
	mean.beta[gen] <- sum(betaFreq*meshpts.beta); 
	betaDev <- meshpts.beta - mean.beta[gen];
	var.beta[gen] <- sum(betaFreq*betaDev^2); 
	
	mean.beta2[gen] <- sum(betaFreq*beta2); #beta2 = meshpts.beta^2 
	beta2Dev <- beta2 - mean.beta2[gen] 
	var.beta2[gen] <- sum(betaFreq*beta2Dev^2); 
	
	if(gen%%10==0) cat(gen,"  ",mean.beta[gen],"\n")
}

params.yr <- as.matrix(params.yr); 

return(list(mean.beta=mean.beta,var.beta=var.beta,mean.beta2=mean.beta2,var.beta2=var.beta2,
			nt=nt,betaDist=betaDist,params.yr=params.yr))
}
