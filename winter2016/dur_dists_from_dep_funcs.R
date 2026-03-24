
###########################################################
# Functions that take the dependence parameters estimated #
# using the different approaches and simulate tail chains #
###########################################################

SimMCThetaCombo <- function(nrep,n,thresh,gamma,sim.exc.type=1,n.sig.u=5,n.xi=0.2){
  # Compute values of theta and thetac given in equation (6) of heatwaves paper
  # given logistic dependence parameter. Can implement importance sampling.
  #
  # Args:
  #   nrep: Number of chains simulated to calculate theta and thetac
  #   n: The length of each chain to be simulated
  #   thresh: Threshold on the Gumbel scale
  #   gamma: Logistic dependence parameter
  #   sim.exc.type: Set as 1 for regular method of generating exceedances from
  #                 an exponential distribution; set as 2 for importance sampling
  #                 exceedances drawn from distribution with heavier tail
  #   n.sig.u: If sim.exc.type==2, scale paramter of GPD to simulate exceedance from
  #   n.xi: If sim.exc.type==2, shape paramter of GPD to simulate exceedance from
  #
  # Returns:
  #   Estimates of theta and thetac for the given logistic dependence parameter.
  #
  # Error handling:
  if (gamma>1 && gamma<0){
    stop("Logistic dependence parameter must lie in the region (0,1)")
  }
  if (nrep<= 0){
    stop("Must have positive number of replications")
  }
  if (n<= 0){
    stop("Chain replications have positive length")
  }
  # Function START
  X <- array(NA,dim=c(n,nrep))
  if (sim.exc.type == 1){   # No importance sampling
    X[1,] <- thresh + rexp(nrep)
  } else if (sim.exc.type == 2){   # With importance sampling
    sim.exc <- rgpd(nrep,scale=n.sig.u,shape=n.xi)
    X[1,] <- thresh + sim.exc
    gpd.correct <- dexp(sim.exc,rate=1)/dgpd(sim.exc,scale=n.sig.u,shape=n.xi)
  }
  
  for (i in 1:(n-1)){
    X[i+1,] <- FwdStepSM(X[i,],dep=gamma,nrep=nrep)
  }
  
  phi <- apply(X=X>thresh,MARGIN=2,FUN=sum)
  consec.vals <- apply(X=X>thresh,MARGIN=2,FUN=decluster.runs,r=1)
  chi <- unlist(lapply(X=lapply(X=consec.vals,FUN=SizeFromList),FUN=MaxValIfFirst))
  
  if (sim.exc.type == 1){
    theta <- c(tabulate(phi),rep(0,length.out=n-length(tabulate(phi))))/nrep
    thetac <- c(tabulate(chi),rep(0,length.out=n-length(tabulate(chi))))/nrep
  } else if (sim.exc.type == 2){
    theta <- thetac <- numeric(n)
    for (k in 1:n){
      theta[k] <- sum(gpd.correct[phi==k])/nrep
      thetac[k] <- sum(gpd.correct[chi==k])/nrep
    }
  }
  
  return( cbind(theta,thetac) )
  
}


Dep2PiSM <- function(dep,nrep,chain.length,thresh,exc.type){
  # Compute estimates for pi and pic with logistic dependence structure 
  # using method from Smith (1992). Run the function 'SimMCThetaCombo' 
  # to compute values of theta and thetac and apply PAV method
  # to obtain decreasing function which can be used to obtain pi and pic.
  #
  # Args:
  #   dep: Logistic dependence parameter
  #   nrep: Number of chains simulated to calculate theta and thetac
  #   chain.length: The length of each chain to be simulated
  #   thresh: Threshold on the Gumbel scale
  #   exc.thresh: Variable to be passed through to 'SimMCThetaCombo', 
  #               set as 2 for importance sampling, 1 otherwise.
  #
  # Returns:
  #   Estimates of the distributions pi and pic for given dependence parameter
  #
  # Error handling:
  if (dep>1 && dep<0){
    stop("Logistic dependence parameter must lie in the region (0,1)")
  }
  if (nrep<= 0){
    stop("Must have positive number of replications")
  }
  if (chain.length<= 0){
    stop("Chain replications have positive length")
  }
  # Function START
  # Generate the theta values
  theta.mat <- SimMCThetaCombo(nrep,chain.length,thresh,dep,exc.type)  
  theta.mat <- apply(X=theta.mat,MARGIN=2,FUN=pava,decreasing=TRUE)
  # Need to have check in place as PAV technique can lead to negative values
  theta.mat[theta.mat<0] <- 0
  pi.dists <- array(0,dim=c(chain.length,2))  
  for (k in 1:(chain.length-1)){
    # Apply differencing formulas in equation (8) of heatwave paper
    pi.dists[k,1] <- (theta.mat[k,1]-theta.mat[k+1,1])/(theta.mat[1,1])
    pi.dists[k,2] <- (theta.mat[k,2]-theta.mat[k+1,2])/(theta.mat[1,2])
  }  
  
  return( pi.dists )
  
}


FwdStepSM <- function(x,dep,nrep){
  # Function that steps a MC chain forwards with transition probability
  # given in eqn (20) of my paper
  #
  # Args:
  #   x: The value of the chain at the current time
  #   dep: The logistic dependence parameter
  #   nrep: The number of replications which controls how many new uniform values to simulate
  #
  # Returns:
  #   out.put: The value of the chain at the next timestep
  #
  # Function START
  out.put <- x - (dep)*log(runif(nrep)^(1/(dep-1))-1)
  
  return( out.put )
  
}


MaxValIfFirst <- function(x){
  # Function that takes the size of the consecutive clusters defined by the runs
  # estimator and returns a value according to the event C^(i) on p4 of my paper
  #
  # Args:
  #   x: Set of sizes of the consecutive clusters
  #
  # Returns:
  #   The first value unless a longer cluster is observed after the first event
  #   which results in an output of zero
  #
  # Function START
  
  if (max(x)>=x[1]){
    out.put <- x[1]
  } else {
    out.put <- 0
  }
  return( out.put )
}


SizeFromList <- function(x){
  # Function that takes a decluster.runs object and gives just the 
  # size of each of the clusters
  #
  # Args:
  #   x: decluster.runs object
  #
  # Returns:
  #   Array containing the size of each cluster defined using the runs estimator
  #
  # Function START
  
  return( x$size )
}


SimHTThetaCombo <- function(htpar,htGz,nrep,n,thresh,npsm=FALSE,sim.exc.type=1,n.sig.u=2,n.xi=0.4){
  # Compute values of theta and thetac given in equation (6) of heatwaves paper
  # using the Heffernan and Tawn or non-parameteric approach 
  #
  # Args:
  #   htpar: Heffernan and Tawn dependence parameters
  #   htGz: Non-parameteric estimate to the distribution G
  #   nrep: Number of chains simulated to calculate theta and thetac
  #   n: The length of each chain to be simulated
  #   thresh: Threshold on the Gumbel scale
  #   npsm: FALSE to use conditional extremes, TRUE for non-parametric approach.
  #   sim.exc.type: Set as 2 for importance sampling, 1 (default) otherwise
  #   n.sig.u: If sim.exc.type==2, scale paramter of GPD to simulate exceedance from
  #   n.xi: If sim.exc.type==2, shape paramter of GPD to simulate exceedance from
  #
  # Returns:
  #   Estimates of theta and thetac for the given H+T dependence parameters
  #
  # Function START
  
  thresh.l <- qlaplace(thresh)
  
  X <- array(NA,dim=c(n,nrep))
  
  if (sim.exc.type == 1){
    X[1,] <- thresh.l + rexp(nrep)
  } else if (sim.exc.type == 2){
    sim.exc <- rgpd(nrep,scale=n.sig.u,shape=n.xi)
    X[1,] <- thresh.l + sim.exc
    gpd.correct <- dexp(sim.exc,rate=1)/dgpd(sim.exc,scale=n.sig.u,shape=n.xi)
  }
  
  if (npsm==FALSE){
    for (i in 1:(n-1)){
      X[i+1,] <- FwdStepHT(X[i,],dep=htpar,z.samp=htGz,nrep=nrep)[[1]]
    }
  } else if (npsm==TRUE){
    for (i in 1:(n-1)){
      X[i+1,] <- FwdStepNP(X[i,],dep=htpar,z.samp=htGz,nrep=nrep)
    }
  }
  
  phi <- apply(X=X>thresh.l,MARGIN=2,FUN=sum)
  consec.vals <- apply(X=X>thresh.l,MARGIN=2,FUN=decluster.runs,r=1)
  chi <- unlist(lapply(X=lapply(X=consec.vals,FUN=SizeFromList),FUN=MaxValIfFirst))
  
  if (sim.exc.type == 1){
    theta <- c(tabulate(phi),rep(0,length.out=n-length(tabulate(phi))))/nrep
    thetac <- c(tabulate(chi),rep(0,length.out=n-length(tabulate(chi))))/nrep
  } else if (sim.exc.type == 2){
    theta <- thetac <- numeric(n)
    for (k in 1:n){
      theta[k] <- sum(gpd.correct[phi==k])/nrep
      thetac[k] <- sum(gpd.correct[chi==k])/nrep
    }
  }
  
  return( cbind(theta,thetac) )
  
}


Dep2PiHT <- function(dep,zsamp,nrep,chain.length,thresh,exc.type,npsm=FALSE,lag.2=FALSE){
  # Compute estimates for pi and pic with logistic dependence structure using method 
  # from Heffernan and Tawn (2004). Run the function SimHTThetaCombo to compute 
  # values of theta and thetac and apply PAV method
  # to obtain decreasing function which can be used to obtain pi and pic.
  #
  # Args:
  #   dep: Dependence parameters (alpha,beta,mu,sigma)
  #   zsamp: Non-parametric estimate to distribution G
  #   nrep: Number of chains simulated to calculate theta and thetac
  #   chain.length: The length of each chain to be simulated
  #   thresh: Threshold on Uniform scale
  #   exc.thresh: Variable to be passed through to 'SimHTThetaCombo', 
  #               set as 2 for importance sampling, 1 otherwise.
  #   npsm: Variable to be passed through to 'SimHTThetaCombo', set
  #         as FALSE to use conditional extremes, TRUE for 
  #         non-parametric approach.
  #   lag.2: Should lag-2 information be incoporated?
  #
  # Returns:
  #   Estimates of the distributions pi and pic for given dependence parameter
  #
  # Error handling:
  if (nrep<= 0){
    stop("Must have positive number of replications")
  }
  if (chain.length<= 0){
    stop("Chain replications have positive length")
  }
  # Function START
  
  if (lag.2==TRUE){
    theta.mat <- SimHTThetaComboWithLag2(dep,zsamp,nrep,chain.length,thresh,
                                         sim.exc.type=exc.type,npsm=npsm)
  } else {
    theta.mat <- SimHTThetaCombo(dep,zsamp,nrep,chain.length,thresh,
                                 sim.exc.type=exc.type,npsm=npsm)
  }
  theta.mat <- apply(X=theta.mat,MARGIN=2,FUN=pava,decreasing=TRUE)
  # Need to have check in place as PAV technique can lead to negative values
  theta.mat[theta.mat<0] <- 0 
  pi.dists <- array(0,dim=c(chain.length,2))  
  for (k in 1:(chain.length-1)){
    # Apply differencing formulas in equation (8) of heatwave paper
    pi.dists[k,1] <- (theta.mat[k,1]-theta.mat[k+1,1])/(theta.mat[1,1])
    pi.dists[k,2] <- (theta.mat[k,2]-theta.mat[k+1,2])/(theta.mat[1,2])
  }  
  
  return( pi.dists )
  
}


FwdStepHT <- function(x,dep,z.samp,nrep){
  # Function that steps a MC forwards as in algorithm 3 of my paper
  # (uses Heffernan and Tawn methodology)
  #
  # Args:
  #   x: The value of the chain at the current time
  #   dep: The Heffernan and Tawn dependence parameters (alpha,beta)
  #   z.samp: Non-parametric estimate of the distribution G
  #   nrep: The number of replications which controls how many new uniform values to simulate
  #
  # Returns:
  #   out.put: The value of the chain at the next timestep
  #
  # Function START
  
  z.vals <- qlaplace(sample(x=z.samp,size=sum(x>0),replace=TRUE))
  out.put <- list()
  out.put[[1]] <- numeric(nrep)
  out.put[[1]][x>0] <- dep[1] * pmax(x[x>0],0) + (pmax(x[x>0],0)^(dep[2])) * z.vals
  out.put[[1]][x<=0] <- dep[1] * pmax(x[x<=0],0)
  out.put[[2]] <- numeric(nrep)
  out.put[[2]][x>0] <- z.vals
  out.put[[2]][x<=0] <- NA
  
  return( out.put )
}


FwdStepNP <- function(x,dep,z.samp,nrep){
  # Function that steps a MC forwards as in algorithm 2 of my paper
  # (uses non-parametric methodology)
  #
  # Args:
  #   x: The value of the chain at the current time
  #   dep: The Heffernan and Tawn dependence parameters (alpha=1,beta=0)
  #   z.samp: Non-parametric estimate of the distribution G obtained as differences of original chain at consecutive time steps
  #   nrep: The number of replications which controls how many new uniform values to simulate
  #
  # Returns:
  #   out.put: The value of the chain at the next timestep
  #
  # Function START
  
  #if (dep[1]!=1 || dep[2]!=0){
  #  stop("Need dependence parameters alpha=1 and beta=0 for asymptotically dependent method")
  #}
  
  out.put <- numeric(nrep)
  out.put[x>0] <- dep[1] * pmax(x[x>0],0) + (pmax(x[x>0],0)^(dep[2])) * sample(x=z.samp,size=sum(x>0),replace=TRUE)
  out.put[x<=0] <- dep[1] * pmax(x[x<=0],0)
  
  return( out.put )
}


Pi2ExceedProb <- function(distr,exceed.val,p,n,u,v){
  # Function that takes a duration distribution and gives the probability
  # of observing at least one event with length exceed.val in the time
  # period n
  #
  # Args:
  #   distr: Duration distribution of interest
  #   exceed.val: Interested in probability of observing at least one event with duration exceed.val
  #   p: Set of marginal parameters in the form (sig.u,xi,lambda.u)
  #   n: Time period of interest in days (92 days will correspond to one year since only interested in JJA)
  #   u: Modelling threshold
  #   v: Critical level
  #
  # Returns:
  #   The probability of observing at least one event with length exceed.val in the time period n.
  #   Also outputs the extremal index and mean number of clusters at level v.
  #
  # Function START  
  ext.ind.v <- 1/(sum((1:length(distr)*distr)))
  tau.v <- ext.ind.v*n*p[3]*(1+p[2]*((v-u)/p[1]))^(-1/p[2])
  exceed.prob <- 1 - exp(-tau.v*sum(distr[exceed.val:length(distr)]))  # What is the probability of observing at least one cluster with exceedVal exceedances or larger?
  
  if (n==92){
    ret.period <- 1/exceed.prob
  } else {
    print("Looking at probability over more than one year will cause NA return period")
    ret.period <- NA
  }
  
  return( c(ext.ind.v,tau.v,exceed.prob,ret.period) )
  
}


FwdBwdSimWithIntensity <- function(fwdpar,fwdGz,bwdpar,bwdGz,n,thresh,peakVal,nrep){
  # Compute probabilities P(N=i|M=Mcl) and P(NC=i|M=Mcl) using 
  # peak value estimation technique
  #
  # Args:
  #   fwdpar: H+T dependence parameters from a fit of X_{t+1}|X_{t}>thresh
  #   fwdGz: Empirical estimate of the distribution G for the above forward fit
  #   bwdpar: H+T dependence parameters from a fit of X_{t}|X_{t+1}>thresh
  #   bwdGz: Empirical estimate of the distribution G for the above backward fit
  #   n: The length of each chain to be simulated, final output clusters will have length 2n-1
  #   thresh: Critical level on uniform scale
  #   peakVal: Maximum value of the replicate clusters to be simulated
  #   nrep: The number of replicate clusters simulated to generate P(N=i|M=Mcl) and P(NC=i|M=Mcl)
  #
  # Returns:
  #   Probabilities P(N=i|M=Mcl) and P(NC=i|M=Mcl) computed from the replicate clusters
  #
  # Function START
  
  thresh.l <- qlaplace(thresh)  # transform the threshold from the uniform scale onto the laplace scale
  
  Xtfwd <- Xtbwd <- z.out.fwd <- z.out.bwd <- array(NA,dim=c(n,nrep))
  
  peakVal.l <- qlaplace(peakVal)   # Convert the cluster maximum to the Laplace margins that we wish to work on
  
  Xtfwd[1,] <- peakVal.l
  
  for (i in 1:(n-1)){
    hold.list <- list()
    hold.list <- FwdStepHT(Xtfwd[i,],dep=fwdpar,z.samp=fwdGz,nrep=nrep)
    Xtfwd[i+1,] <- hold.list[[1]]
    z.out.fwd[i,] <- hold.list[[2]]
  }
  
  while(any(apply(X=Xtfwd,MARGIN=2,FUN=max)>peakVal.l)){
    fail.nos <- which(apply(X=Xtfwd,MARGIN=2,FUN=max)>peakVal.l)
    Xtfwd[,fail.nos] <- 0
    z.out.fwd[,fail.nos] <- NA
    Xtfwd[1,fail.nos] <- peakVal.l
    for (i in 1:(n-1)){
      hold.list <- list()
      hold.list <- FwdStepHT(Xtfwd[i,fail.nos],dep=fwdpar,z.samp=fwdGz,nrep=length(Xtfwd[i,fail.nos]))
      Xtfwd[i+1,fail.nos] <- hold.list[[1]]
      z.out.fwd[i,fail.nos] <- hold.list[[2]]
    }
  }
  
  Xtbwd[1,] <- peakVal.l
  for (i in 1:(n-1)){
    hold.list <- list()
    hold.list <- FwdStepHT(Xtbwd[i,],dep=bwdpar,z.samp=bwdGz,nrep=nrep)
    Xtbwd[i+1,] <- hold.list[[1]]
    z.out.bwd[i,] <- hold.list[[2]]
  }
  
  while(any(apply(X=Xtbwd,MARGIN=2,FUN=max)>peakVal.l)){
    fail.nos <- which(apply(X=Xtbwd,MARGIN=2,FUN=max)>peakVal.l)
    Xtbwd[,fail.nos] <- 0
    z.out.bwd[,fail.nos] <- NA
    Xtbwd[1,fail.nos] <- peakVal.l
    for (i in 1:(n-1)){
      hold.list <- list()
      hold.list <- FwdStepHT(Xtbwd[i,fail.nos],dep=bwdpar,z.samp=bwdGz,nrep=length(Xtbwd[i,fail.nos]))
      Xtbwd[i+1,fail.nos] <- hold.list[[1]]
      z.out.bwd[i,fail.nos] <- hold.list[[2]]
    }
  }
  
  Xtsim <- rbind(apply(X=Xtbwd,MARGIN=2,FUN=rev),Xtfwd[-1,])
  
  phi <- apply(X=Xtsim>thresh.l,MARGIN=2,FUN=sum)
  consec.vals <- apply(X=Xtsim>thresh.l,MARGIN=2,FUN=decluster.runs,r=1)   # Only want the cluster sizes coming out here!
  chi <- unlist(lapply(X=lapply(X=consec.vals,FUN=SizeFromList),FUN=max))
  
  pNcM <- c(tabulate(phi),rep(0,length.out=2*n-length(tabulate(phi))))/nrep
  pNCcM <- c(tabulate(chi),rep(0,length.out=2*n-length(tabulate(chi))))/nrep
  
  out.put <- list()
  out.put[[1]] <- pNcM
  out.put[[2]] <- pNCcM
  out.put[[3]] <- Xtsim
  out.put[[4]] <- z.out.fwd
  out.put[[5]] <- z.out.bwd
  
  
  return( out.put )   # Remember to return the output of the code!
  
}


FwdBwdSimWithIntensityNP <- function(fwdpar,fwdGz,bwdpar,bwdGz,n,thresh,peakVal,nrep){
  # Compute probabilities P(N=i|M=Mcl) and P(NC=i|M=Mcl) using 
  # peak value estimation technique
  #
  # Args:
  #   fwdpar: H+T dependence parameters from a fit of X_{t+1}|X_{t}>thresh
  #   fwdGz: Empirical estimate of the distribution G for the above forward fit
  #   bwdpar: H+T dependence parameters from a fit of X_{t}|X_{t+1}>thresh
  #   bwdGz: Empirical estimate of the distribution G for the above backward fit
  #   n: The length of each chain to be simulated, final output clusters will have length 2n-1
  #   thresh: Critical level on uniform scale
  #   peakVal: Maximum value of the replicate clusters to be simulated
  #   nrep: The number of replicate clusters simulated to generate P(N=i|M=Mcl) and P(NC=i|M=Mcl)
  #
  # Returns:
  #   Probabilities P(N=i|M=Mcl) and P(NC=i|M=Mcl) computed from the replicate clusters
  #
  # Function START
  
  thresh.l <- qlaplace(thresh)  # transform the threshold from the uniform scale onto the laplace scale
  
  Xtfwd <- Xtbwd <- array(NA,dim=c(n,nrep))
  
  peakVal.l <- qlaplace(peakVal)   # Convert the cluster maximum to the Laplace margins that we wish to work on
  
  Xtfwd[1,] <- peakVal.l
  for (i in 1:(n-1)){
    Xtfwd[i+1,] <- FwdStepNP(Xtfwd[i,],dep=fwdpar,z.samp=fwdGz,nrep=nrep)
  }
  
  while(any(apply(X=Xtfwd,MARGIN=2,FUN=max)>peakVal.l)){
    fail.nos <- which(apply(X=Xtfwd,MARGIN=2,FUN=max)>peakVal.l)
    Xtfwd[,fail.nos] <- 0
    Xtfwd[1,fail.nos] <- peakVal.l
    for (i in 1:(n-1)){
      Xtfwd[i+1,fail.nos] <- FwdStepNP(Xtfwd[i,fail.nos],dep=fwdpar,z.samp=fwdGz,nrep=length(Xtfwd[i,fail.nos]))
    }
  }
  
  Xtbwd[1,] <- peakVal.l
  for (i in 1:(n-1)){
    Xtbwd[i+1,] <- FwdStepNP(Xtbwd[i,],dep=bwdpar,z.samp=bwdGz,nrep=nrep)
  }
  
  while(any(apply(X=Xtbwd,MARGIN=2,FUN=max)>peakVal.l)){
    fail.nos <- which(apply(X=Xtbwd,MARGIN=2,FUN=max)>peakVal.l)
    Xtbwd[,fail.nos] <- 0
    Xtbwd[1,fail.nos] <- peakVal.l
    for (i in 1:(n-1)){
      Xtbwd[i+1,fail.nos] <- FwdStepNP(Xtbwd[i,fail.nos],dep=bwdpar,z.samp=bwdGz,nrep=length(Xtbwd[i,fail.nos]))
    }
  }
  
  Xtsim <- rbind(apply(X=Xtbwd,MARGIN=2,FUN=rev),Xtfwd[-1,])
  
  phi <- apply(X=Xtsim>thresh.l,MARGIN=2,FUN=sum)
  consec.vals <- apply(X=Xtsim>thresh.l,MARGIN=2,FUN=decluster.runs,r=1)   # Only want the cluster sizes coming out here!
  chi <- unlist(lapply(X=lapply(X=consec.vals,FUN=SizeFromList),FUN=max))
  
  pNcM <- c(tabulate(phi),rep(0,length.out=2*n-length(tabulate(phi))))/nrep
  pNCcM <- c(tabulate(chi),rep(0,length.out=2*n-length(tabulate(chi))))/nrep
  
  out.put <- list()
  out.put[[1]] <- cbind(pNcM,pNCcM)
  out.put[[2]] <- Xtsim
  
  return( out.put )   # Remember to return the output of the code!
  
}


FwdBwdSimWithIntensitySmith <- function(fwdpar,fwdGz,bwdpar,bwdGz,n,thresh,peakVal,nrep){
  # Compute probabilities P(N=i|M=Mcl) and P(NC=i|M=Mcl) using 
  # peak value estimation technique
  #
  # Here use the parametric approach of Smith (1992)
  #
  # Args:
  #   fwdpar: H+T dependence parameters from a fit of X_{t+1}|X_{t}>thresh
  #   fwdGz: Empirical estimate of the distribution G for the above forward fit
  #   bwdpar: H+T dependence parameters from a fit of X_{t}|X_{t+1}>thresh
  #   bwdGz: Empirical estimate of the distribution G for the above backward fit
  #   n: The length of each chain to be simulated, final output clusters will have length 2n-1
  #   thresh: Critical level on uniform scale
  #   peakVal: Maximum value of the replicate clusters to be simulated
  #   nrep: The number of replicate clusters simulated to generate P(N=i|M=Mcl) and P(NC=i|M=Mcl)
  #
  # Returns:
  #   Probabilities P(N=i|M=Mcl) and P(NC=i|M=Mcl) computed from the replicate clusters
  #
  # Function START
  
  thresh.g <- -log(-log(thresh))  # transform the threshold from the uniform scale onto the Gumbel scale
  
  Xtfwd <- Xtbwd <- array(NA,dim=c(n,nrep))
  
  peakVal.g <- -log(-log(peakVal))  # Convert the cluster maximum to the Gumbel margins that we wish to work on
  
  Xtfwd[1,] <- peakVal.g
  for (i in 1:(n-1)){
    Xtfwd[i+1,] <- FwdStepSM(x = Xtfwd[i,],dep=fwdpar,nrep=nrep)
  }
  
  while(any(apply(X=Xtfwd,MARGIN=2,FUN=max)>peakVal.g)){
    fail.nos <- which(apply(X=Xtfwd,MARGIN=2,FUN=max)>peakVal.g)
    Xtfwd[,fail.nos] <- 0
    Xtfwd[1,fail.nos] <- peakVal.g
    for (i in 1:(n-1)){
      Xtfwd[i+1,fail.nos] <- FwdStepSM(x = Xtfwd[i,fail.nos],dep = fwdpar,nrep = length(Xtfwd[i,fail.nos]))
    }
  }
  
  Xtbwd[1,] <- peakVal.g
  for (i in 1:(n-1)){
    Xtbwd[i+1,] <- FwdStepSM(Xtbwd[i,],dep=bwdpar,nrep=nrep)
  }
  
  while(any(apply(X=Xtbwd,MARGIN=2,FUN=max)>peakVal.g)){
    fail.nos <- which(apply(X=Xtbwd,MARGIN=2,FUN=max)>peakVal.g)
    Xtbwd[,fail.nos] <- 0
    Xtbwd[1,fail.nos] <- peakVal.g
    for (i in 1:(n-1)){
      Xtbwd[i+1,fail.nos] <- FwdStepSM(x = Xtbwd[i,fail.nos],dep = bwdpar,nrep = length(Xtbwd[i,fail.nos]))
    }
  }
  
  Xtsim <- rbind(apply(X=Xtbwd,MARGIN=2,FUN=rev),Xtfwd[-1,])
  
  phi <- apply(X=Xtsim>thresh.g,MARGIN=2,FUN=sum)
  consec.vals <- apply(X=Xtsim>thresh.g,MARGIN=2,FUN=decluster.runs,r=1)   # Only want the cluster sizes coming out here!
  chi <- unlist(lapply(X=lapply(X=consec.vals,FUN=SizeFromList),FUN=max))
  
  pNcM <- c(tabulate(phi),rep(0,length.out=2*n-length(tabulate(phi))))/nrep
  pNCcM <- c(tabulate(chi),rep(0,length.out=2*n-length(tabulate(chi))))/nrep
  
  out.put <- list()
  out.put[[1]] <- cbind(pNcM,pNCcM)
  out.put[[2]] <- Xtsim
  
  return( out.put )   # Remember to return the output of the code!
  
}
