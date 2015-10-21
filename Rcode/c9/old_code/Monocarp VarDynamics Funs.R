# Implement a 2-dimensional IPM 
# Each column of nt is the density of z for a given flowering intercept

iterate_model <- function(params,n.iter,init.beta.mean,init.beta.sd,meshpts.beta,meshpts.z) {

    nt <- matrix(NA,nrow=nBigMatrix.z,ncol=nBigMatrix.beta)
    z.s <- dnorm(meshpts.z, mean = 2.9, sd = m.par.true["rcsz.sd"])
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

	for(i in 1:nBigMatrix.beta){
		params["flow.int"] <- meshpts.beta[i]
		P.beta[i,1:nBigMatrix.z,1:nBigMatrix.z] <- h.z * (outer(meshpts.z, meshpts.z, P_z1z, m.par = 			params))
		Seeds.beta[i,1:nBigMatrix.z] <- Seeds_z(meshpts.z, params)
	}

	#Inheritance kernel
	M <- matrix(NA,ncol=nBigMatrix.beta,nrow=nBigMatrix.beta)
	for(i in 1:nBigMatrix.beta){
 		 M[,i] <- dnorm(meshpts.beta,mean=meshpts.beta[i],sd=beta.off.sd)
  	}

	#pdf of offspring sizes
	off.pdf <- c_0z1(meshpts.z,m.par.true)
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
	#calculate seeds produced by beta_i plants
	seeds.from.betai <- sapply(index,function(i) h.z*sum(Seeds.beta[i,] %*% nt[,i]))

	#redistribute seeds according to inheritance kernel
	seeds.with.betai <- h.beta * (M %*% seeds.from.betai)

	f.recruits.with.betai <- seeds.with.betai/sum(seeds.from.betai)

	nt1 <- sapply(index,function(i) P.beta[i,,] %*% nt[,i] + Recr * off.pdf * f.recruits.with.betai[i])

	n.surv <- sapply(index,function(i) P.beta[i,,] %*% nt[,i])
	n.surv <- sum(n.surv) * h.z * h.beta
	nt <- nt1
	
	prop.Recr[gen] <- Recr/(Recr+n.surv)
	betaDist[,gen] = apply(nt,2,sum); betaDist[,gen] <- betaDist[,gen]/sum(betaDist[,gen]); 
	betaFreq <- betaDist[,gen] 
	
	#### Compute variances of beta and beta^2 using formula Var=E[(X-Xbar)^2] 
	mean.beta[gen] <- sum(betaFreq*meshpts.beta); 
	betaDev <- meshpts.beta - mean.beta[gen];
	var.beta[gen] <- sum(betaFreq*betaDev^2); 
	
	mean.beta2[gen] <- sum(betaFreq*beta2); #beta2 = meshpts.beta^2 
	beta2Dev <- beta2 - mean.beta2[gen] 
	var.beta2[gen] <- sum(betaFreq*beta2Dev^2); 
	
	if(gen%%100==0) cat(gen,"  ",mean.beta[gen],"   ",Recr/(Recr+n.surv),"\n")
}


return(list(mean.beta=mean.beta,var.beta=var.beta,mean.beta2=mean.beta2,var.beta2=var.beta2,
			prop.Recr=prop.Recr,nt=nt,betaDist=betaDist))
}


###################################################################################################
#Approximate dynamics using Iwasa et al. 1991 Evolution
####################################################################################################
R0_calc <- function (params) {
	IPM.true <- mk_K(nBigMatrix, params, L.z, U.z)
	# to keep close to the formulae in the text next we define the F and P iteration matricies
    P <- IPM.true$P;  F <- IPM.true$F;

    # Fundamental operator 
    N <- solve(diag(nBigMatrix)-P); 

    # Compute R0 as dominant eigenvalue of FN
    R <- F %*% N
    R0 <- abs(eigen(R,only.values=TRUE)$values[1])
    return(R0)
}

Approx_dynamics <- function(params,init.mean.flow.int,n.iter,delta,var.beta.t,prop.Recr) {

# params=m.par.true; init.mean.flow.int=-20; n.iter=10; delta=0.01; 
# var.beta.t=sol.2D.minus20$var.beta;prop.Recr=sol.2D.minus20$prop.Recr;

    params["flow.int"]<- init.mean.flow.int
	beta.t.mean <- rep(NA,n.iter)
	beta.t.var <- rep(NA,n.iter)
	S.t.mean <- rep(NA,n.iter) #this is mean of S=beta^2
	
	beta.t.mean[1] <- init.mean.flow.int
	beta.t.var[1] <- init.beta.sd^2
	S.t.mean[1] <- init.beta.sd^2 + init.mean.flow.int^2
	
for(gen in 2:n.iter){
    #Calculate p.r so the current mean strategy is at equilibrium lambda=1
	params["p.r"] <- 1
	equ.p.r <- 1/R0_calc(params)
	params["p.r"] <- equ.p.r

# Calculate the derivative of lambda with respect to beta
	int.0 <- params["flow.int"]; 
	params["flow.int"] <- int.0 + delta
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	lambda.up <- Re(eigen(IPM.kernel$K,only.values = TRUE)$values[1])
	params["flow.int"] <- int.0 - delta
	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	lambda.down <- Re(eigen(IPM.kernel$K,only.values = TRUE)$values[1])
	params["flow.int"] <- int.0

	d.lambda.d.beta <- (lambda.up - lambda.down)/(2*delta)

# Calculate the mean intercept
	#We're cheating here as we're plugging in the var beta from the 2D IPM, var.beta.t[gen-1] 
	beta.t.mean[gen] <- beta.t.mean[gen-1] + d.lambda.d.beta * var.beta.t[gen-1]

# Genetic variance of beta t-1
	G.beta <- S.t.mean[gen-1] - beta.t.mean[gen-1]^2
	# G.beta <- var.beta.t[gen-1]
	
#Genetic variance of beta^2
#	G.S <- beta.t.mean[gen-1]^4 + 6*(beta.t.mean[gen-1] ^2)*G.beta + 3*G.beta*G.beta - S.t.mean[gen-1]^2
    G.S <- 4 * (beta.t.mean[gen-1]^2)*G.beta; 
	
#Mean S at time t
	den <- (-2*sqrt(abs(S.t.mean[gen-1])))  # CHANGE: sqrt(sbar) is not the same as xbar. 
	S.t.mean[gen] <- S.t.mean[gen-1] + G.S*d.lambda.d.beta/(2*beta.t.mean[gen-1]) 
	params["flow.int"] <- beta.t.mean[gen] 
	if(gen%%10==0) cat(gen,"  ",beta.t.mean[gen-1],"  ",S.t.mean[gen-1],"\n")
		
}

return(list(beta.t.mean=beta.t.mean,S.t.mean=S.t.mean))
}

Approx_dynamics2 <- function(params,init.mean.flow.int,n.iter,delta,var.beta.t,prop.Recr) {

	params["flow.int"]<- init.mean.flow.int
	beta.t.mean <- rep(NA,n.iter)
	beta.t.var <- rep(NA,n.iter)
	S.t.mean <- rep(NA,n.iter) #this is mean S
	
	beta.t.mean[1] <- init.mean.flow.int
	beta.t.var[1] <- init.beta.sd*init.beta.sd
	S.t.mean[1] <- init.beta.sd*init.beta.sd + init.mean.flow.int*init.mean.flow.int
	

for(gen in 2:n.iter){

#Calculate p.r so the current mean strategy is at equilibrium lambda=1
	params["p.r"] <- 1
	equ.p.r <- 1/R0_calc(params)
	params["p.r"] <- equ.p.r

	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	b0 <- params["flow.int"]; 
	params["flow.int"] <- b0 + delta
	Kplus <- mk_K(nBigMatrix, params, L.z, U.z)
	params["flow.int"] <- b0 - delta
	Kminus <- mk_K(nBigMatrix, params, L.z, U.z)
	params["flow.int"] <- b0 	

#Calculate the derivative of lambda with respect to beta
	lambda.up <- Re(eigen(Kplus$K,only.values = TRUE)$values[1])
	lambda.down <- Re(eigen(Kminus$K,only.values = TRUE)$values[1])
	d.lambda.d.beta <- (lambda.up - lambda.down)/(2*delta)
	
#Calculate the derivative of W - survival only - with respect to beta
	eigen.sys <- eigen(IPM.kernel$K)
	w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))
	ave.surv.Up <- sum(Kplus$P %*% w)
	ave.surv.Down <- sum(Kminus$P %*% w)
	d.W.surv.d.beta <- (ave.surv.Up - ave.surv.Down)/(2*delta)
	W.surv <- sum(IPM.kernel$P %*% w)
	
#Calculate the derivative of W - fecundity only - with respect to beta
	ave.fec.Up <- sum(Kplus$F %*% w)
	ave.fec.Down <- sum(Kminus$F %*% w)
	d.W.fec.d.beta <- (ave.fec.Up - ave.fec.Down)/(2*delta)
	W.fec <- sum(IPM.kernel$F%*%w)

#Calculate the mean intercept

	#We're cheating here as we're plugging in the var beta from the 2 D IPM var.beta.t[gen-1] 
	beta.t.mean[gen] <- beta.t.mean[gen-1] + d.lambda.d.beta*var.beta.t[gen-1]

#mean beta after survival
	beta.th.mean <- beta.t.mean[gen-1] + d.W.surv.d.beta * var.beta.t[gen-1] / W.surv

#Genetic variance of beta t-2	
	G.beta <- S.t.mean[gen-1] - beta.t.mean[gen-1] * beta.t.mean[gen-1]

#Genetic variance of beta^2
	# G.S <- beta.t.mean[gen-1]^4 + 6 * beta.t.mean[gen-1] ^2 * G.beta + 3 * G.beta * G.beta - S.t.mean[gen-1] * S.t.mean[gen-1]
    # G.S <- 4 * (beta.t.mean[gen-1]^2)*G.beta; 
	 G.S <- 4*(var.beta.t[gen-1])*beta.t.mean[gen]; #cheating 
	
#mean S in recruits	
	S.t.mean.Recr <- S.t.mean[gen-1] + G.S*d.W.fec.d.beta/(2*beta.t.mean[gen-1]*W.fec)

#mean S in survivors	
	S.t.mean.Surv <- S.t.mean[gen-1] + G.S*d.W.surv.d.beta/(2*beta.t.mean[gen-1]*W.surv)

#mean S in population	
	S.t.mean[gen] <- prop.Recr[gen] * S.t.mean.Recr + (1-prop.Recr[gen]) * S.t.mean.Surv

#Variance of beta in recruits, note addition of mutational variance at end
	Var.Recr <- S.t.mean.Recr - beta.t.mean[gen]*beta.t.mean[gen] + beta.off.sd * beta.off.sd

#Variance of beta in survivors
	Var.Surv <- S.t.mean.Surv - beta.th.mean*beta.th.mean

#overall variance of beta at time t
	beta.t.var[gen] <- prop.Recr[gen] * Var.Recr + (1-prop.Recr[gen]) * Var.Surv +
	 				prop.Recr[gen]*(1-prop.Recr[gen])*(beta.t.mean[gen]-beta.th.mean)^2
		
	params["flow.int"] <- beta.t.mean[gen] 
	
	cat(gen,"  ",beta.t.mean[gen],"  ",beta.t.var[gen],"  ",var.beta.t[gen],"\n")
		
}

return(list(beta.t.mean=beta.t.mean,beta.t.var=beta.t.var,S.t.mean=S.t.mean))
}


Approx_dynamics2 <- function(params,init.mean.flow.int,n.iter,delta,var.beta.t,prop.Recr) {

	params["flow.int"]<- init.mean.flow.int
	beta.t.mean <- rep(NA,n.iter)
	beta.t.var <- rep(NA,n.iter)
	S.t.mean <- rep(NA,n.iter) #this is mean S
	
	beta.t.mean[1] <- init.mean.flow.int
	beta.t.var[1] <- init.beta.sd*init.beta.sd
	S.t.mean[1] <- init.beta.sd*init.beta.sd + init.mean.flow.int*init.mean.flow.int
	

for(gen in 2:n.iter){

#Calculate p.r so the current mean strategy is at equilibrium lambda=1
	params["p.r"] <- 1
	equ.p.r <- 1/R0_calc(params)
	params["p.r"] <- equ.p.r

	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	b0 <- params["flow.int"]; 
	params["flow.int"] <- b0 + delta
	Kplus <- mk_K(nBigMatrix, params, L.z, U.z)
	params["flow.int"] <- b0 - delta
	Kminus <- mk_K(nBigMatrix, params, L.z, U.z)
	params["flow.int"] <- b0 	

#Calculate the derivative of lambda with respect to beta
	lambda.up <- Re(eigen(Kplus$K,only.values = TRUE)$values[1])
	lambda.down <- Re(eigen(Kminus$K,only.values = TRUE)$values[1])
	d.lambda.d.beta <- (lambda.up - lambda.down)/(2*delta)
	
#Calculate the derivative of W - survival only - with respect to beta
	eigen.sys <- eigen(IPM.kernel$K)
	w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))
	ave.surv.Up <- sum(Kplus$P %*% w)
	ave.surv.Down <- sum(Kminus$P %*% w)
	d.W.surv.d.beta <- (ave.surv.Up - ave.surv.Down)/(2*delta)
	W.surv <- sum(IPM.kernel$P %*% w)
	
#Calculate the derivative of W - fecundity only - with respect to beta
	ave.fec.Up <- sum(Kplus$F %*% w)
	ave.fec.Down <- sum(Kminus$F %*% w)
	d.W.fec.d.beta <- (ave.fec.Up - ave.fec.Down)/(2*delta)
	W.fec <- sum(IPM.kernel$F%*%w)

#Calculate the mean intercept

	#We're cheating here as we're plugging in the var beta from the 2 D IPM var.beta.t[gen-1] 
	beta.t.mean[gen] <- beta.t.mean[gen-1] + d.lambda.d.beta*var.beta.t[gen-1]

#mean beta after survival
	beta.th.mean <- beta.t.mean[gen-1] + d.W.surv.d.beta*var.beta.t[gen-1]/W.surv

#Genetic variance of beta t	
	G.beta <- S.t.mean[gen-1] - beta.t.mean[gen-1]*beta.t.mean[gen-1]

#Genetic variance of beta^2
	# G.S <- beta.t.mean[gen-1]^4 + 6 * beta.t.mean[gen-1] ^2 * G.beta + 3 * G.beta * G.beta - S.t.mean[gen-1] * S.t.mean[gen-1]
    # G.S <- 4 * (beta.t.mean[gen-1]^2)*G.beta; 
	 G.S <- 4*(var.beta.t[gen-1])*beta.t.mean[gen]; #cheating 
	
#mean S in recruits	
	S.t.mean.Recr <- S.t.mean[gen-1] + G.S*d.W.fec.d.beta/(2*beta.t.mean[gen-1]*W.fec)

#mean S in survivors	
	S.t.mean.Surv <- S.t.mean[gen-1] + G.S*d.W.surv.d.beta/(2*beta.t.mean[gen-1]*W.surv)

#mean S in population	
	S.t.mean[gen] <- prop.Recr[gen] * S.t.mean.Recr + (1-prop.Recr[gen]) * S.t.mean.Surv

#Variance of beta in recruits, note addition of mutational variance at end
	Var.Recr <- S.t.mean.Recr - beta.t.mean[gen]*beta.t.mean[gen] + beta.off.sd * beta.off.sd

#Variance of beta in survivors
	Var.Surv <- S.t.mean.Surv - beta.th.mean*beta.th.mean

#overall variance of beta at time t
	beta.t.var[gen] <- prop.Recr[gen] * Var.Recr + (1-prop.Recr[gen]) * Var.Surv +
	 				prop.Recr[gen]*(1-prop.Recr[gen])*(beta.t.mean[gen]-beta.th.mean)^2
		
	params["flow.int"] <- beta.t.mean[gen] 
	
	cat(gen,"  ",beta.t.mean[gen],"  ",beta.t.var[gen],"  ",var.beta.t[gen],"\n")
		
}

return(list(beta.t.mean=beta.t.mean,beta.t.var=beta.t.var,S.t.mean=S.t.mean))
}

Approx_dynamics_Var <- function(params,init.mean.flow.int,n.iter,delta,var.beta.t,prop.Recr) {

	params["flow.int"]<- init.mean.flow.int
	beta.t.mean <- rep(NA,n.iter)
	beta.t.var <- rep(NA,n.iter)
	S.t.mean <- rep(NA,n.iter) #this is mean S
	
	beta.t.mean[1] <- init.mean.flow.int
	beta.t.var[1] <- init.beta.sd*init.beta.sd
	S.t.mean[1] <- init.beta.sd*init.beta.sd + init.mean.flow.int*init.mean.flow.int
	
	dlam <- dsurv <- dfec <- rep(NA,n.iter)
	
	change.P.w <- change.F.w <- P.change.w <- F.change.w <- rep(NA,n.iter)

for(gen in 2:n.iter){

#Calculate p.r so the current mean strategy is at equilibrium lambda=1
	params["p.r"] <- 1
	equ.p.r <- 1/R0_calc(params)
	params["p.r"] <- equ.p.r

	IPM.kernel <- mk_K(nBigMatrix, params, L.z, U.z)
	b0 <- params["flow.int"]; 
	params["flow.int"] <- b0 + delta
	Kplus <- mk_K(nBigMatrix, params, L.z, U.z)
	params["flow.int"] <- b0 - delta
	Kminus <- mk_K(nBigMatrix, params, L.z, U.z)
	params["flow.int"] <- b0 	

#Calculate the derivative of lambda with respect to beta
	lambda.up <- Re(eigen(Kplus$K,only.values = TRUE)$values[1])
	lambda.down <- Re(eigen(Kminus$K,only.values = TRUE)$values[1])
	d.lambda.d.beta <- (lambda.up - lambda.down)/(2*delta)
	
	eigen.sys <- eigen(IPM.kernel$K)
	w <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))
	
	eigen.sys <- eigen(Kplus$K)
	wplus <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))

	eigen.sys <- eigen(Kminus$K)
	wminus <- Re(eigen.sys$vectors[,1]/sum(eigen.sys$vectors[,1]))

	
#Calculate the derivative of W - survival only - with respect to beta
	ave.surv <- sum(IPM.kernel$P %*% w)
	ave.surv.Up <- sum(Kplus$P %*% wplus)
	ave.surv.Down <- sum(Kminus$P %*% wminus)
	d.W.surv.d.beta <- (ave.surv.Up - ave.surv.Down)/(2*delta)
	d2.W.surv.d.beta2 <- (ave.surv.Up - 2 * ave.surv + ave.surv.Down)/(delta*delta)
	
	partial.surv <- (Kplus$P - Kminus$P)/(2*delta)

	
#Calculate the derivative of W - fecundity only - with respect to beta
	ave.fec <- sum(IPM.kernel$F %*% w)
	ave.fec.Up <- sum(Kplus$F %*% wplus)
	ave.fec.Down <- sum(Kminus$F %*% wminus)
	d.W.fec.d.beta <- (ave.fec.Up - ave.fec.Down)/(2*delta)
	d2.W.fec.d.beta2 <- (ave.fec.Up - 2 * ave.fec + ave.fec.Down)/(delta*delta)
	
	partial.fec <- (Kplus$F - Kminus$F)/(2*delta)
	
##############################################
change.P.w[gen] <- sum(partial.surv %*% w)
change.F.w[gen] <- sum(partial.fec %*% w)

partial.w <- (wplus - wminus)/(2*delta)

P.change.w[gen] <- sum(IPM.kernel$P %*% partial.w)
F.change.w[gen] <- sum(IPM.kernel$F %*% partial.w)

#store selection coef

dlam[gen]  <- d.lambda.d.beta
dsurv[gen] <- d.W.surv.d.beta
dfec[gen]   <- d.W.fec.d.beta

#Calculate the mean intercept

	#We're cheating here as we're plugging in the var beta from the 2-D IPM var.beta.t[gen-1] 
	beta.t.mean[gen] <- beta.t.mean[gen-1] + d.lambda.d.beta * var.beta.t[gen-1]

#mean beta of survivors
	beta.t.mean.surv <- beta.t.mean[gen-1] + d.W.surv.d.beta*var.beta.t[gen-1]/ave.surv
	
#mean beta in recruits
	beta.t.mean.rec  <- beta.t.mean[gen-1] + d.W.fec.d.beta*var.beta.t[gen-1]/ave.fec
	
	
#sigma 2 in survivors
	beta.t.sigma2.surv <- beta.t.var[gen-1] + (d2.W.surv.d.beta2/ave.surv - (d.W.surv.d.beta/ave.surv)^2) * 	
			beta.t.var[gen-1]*beta.t.var[gen-1]
	
#sigma 2 in recruits
	beta.t.sigma2.rec <- beta.t.var[gen-1] + (d2.W.fec.d.beta2/ave.fec   - (d.W.fec.d.beta/ave.fec)^2)   * 
			beta.t.var[gen-1]*beta.t.var[gen-1]/ave.fec +
            beta.off.sd*beta.off.sd

#overall variance of beta at time t
	beta.t.var[gen] <- prop.Recr[gen] * beta.t.sigma2.rec + (1-prop.Recr[gen]) * beta.t.sigma2.surv +
	 				prop.Recr[gen]*(1-prop.Recr[gen])*(beta.t.mean.surv - beta.t.mean.rec)^2
		
	params["flow.int"] <- beta.t.mean[gen] 
	
	#cat(gen,"  ",beta.t.mean[gen],"  ",beta.t.mean.surv,"  ",beta.t.mean.rec,"\n")
	
	if(gen %% 100 == 0) cat(gen,"  ",beta.t.mean[gen],"  ",beta.t.var[gen],"  ",var.beta.t[gen],"  ",d2.W.surv.d.beta2,"  ",d2.W.fec.d.beta2,"\n")
		
}

return(list(beta.t.mean=beta.t.mean,beta.t.var=beta.t.var,dlam=dlam,dsurv=dsurv,dfec=dfec,
change.P.w=change.P.w,change.F.w=change.F.w,P.change.w=P.change.w,F.change.w=F.change.w))
}





