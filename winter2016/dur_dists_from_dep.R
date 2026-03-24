
##############################################################
# Take the dependence parameters for the original data       #
# and generate duration distributions for each approach.     #                                                
##############################################################

# (1) Using the forward tail chain algorithm from Rootzen (1988)

# Need the threshold on Gumbel margin to apply method
# of Smith (1992)
v.threshg <- -log(-log(v.threshu))
# Obtain the critical level on original margins
v.thresh <- u.thresh + (sig.u/xi)*((((1-v.threshu)/lambda.u)^(-xi))-1)
nrepp <- 100

# Method 1: Conditional extremes approach
piHT.orig <- array(NA,dim=c(40,2,nrepp))

for ( i in 1:nrepp ){
  print(i)
  piHT.orig[,,i] <- Dep2PiHT(dep=ht.dep.fit$par[1:2],zsamp=ht.dep.fit$origz,nrep=20000,
                             chain.length=40,thresh=v.threshu,exc.type=1)
}
piHT.orig <- apply(X = piHT.orig, MARGIN = c(1,2), FUN = mean)

# Method 2: Non-parametric approach
piNP.orig <- array(NA,dim=c(40,2,nrepp))

for ( i in 1:nrepp ){
  print(i)
  piNP.orig[,,i] <- Dep2PiHT(dep=c(1,0),zsamp=npsmGz,nrep=20000,
                             chain.length=40,thresh=v.threshu,exc.type=1,npsm=TRUE)
}
piNP.orig <- apply(X = piNP.orig, MARGIN = c(1,2), FUN = mean)

# Method 3: Parametric approach
piSM.orig <- array(NA,dim=c(40,2,nrepp))
for ( i in 1:nrepp ){
  print(i)
  piSM.orig[,,i] <- Dep2PiSM(dep=gamma,nrep=20000,chain.length=40,thresh=v.threshg,exc.type=1)
}
piSM.orig <- apply(X = piSM.orig, MARGIN = c(1,2), FUN = mean)


# (2) Using the peak value tail chain estimation from Smith et al. (1997)
# Only use the conditional extremes approach in this example

# Obtain the value of P(N=i|M=eta) for the original sample
retPeriod <- 1
v.threshu <- 1 - 1/(92*retPeriod)   # Critical level v on uniform scale

# With a peak value equivalent to 1-in-5 year event
peak.val.rp <- 5
peak.val.u <- 1 - 1/(92*peak.val.rp)   # Peak value on uniform scale

fwdbwdHT5.orig <- FwdBwdSimWithIntensity(fwdpar = ht.dep.fit$par[1:2], fwdGz = ht.dep.fit$origz, 
                                         bwdpar = ht.dep.fit.bwd$par[1:2], bwdGz = ht.dep.fit.bwd$origz, 
                                         n = 40, thresh = v.threshu, peakVal = peak.val.u, 
                                         nrep = 50000)[[1]]

# With a peak value equivalent to 1-in-50 year event
peak.val.rp <- 50
peak.val.u <- 1 - 1/(92*peak.val.rp)   # Peak value on uniform scale

fwdbwdHT50.orig <- FwdBwdSimWithIntensity(fwdpar = ht.dep.fit$par[1:2], fwdGz = ht.dep.fit$origz, 
                                          bwdpar = ht.dep.fit.bwd$par[1:2], bwdGz = ht.dep.fit.bwd$origz, 
                                          n = 40, thresh = v.threshu, peakVal = peak.val.u, 
                                          nrep = 50000)[[1]]

# With a peak value equivalent to 1-in-1000 year event
peak.val.rp <- 1000
peak.val.u <- 1 - 1/(92*peak.val.rp)   # Peak value on uniform scale

fwdbwdHT1000.orig <- FwdBwdSimWithIntensity(fwdpar = ht.dep.fit$par[1:2], fwdGz = ht.dep.fit$origz, 
                                            bwdpar = ht.dep.fit.bwd$par[1:2], bwdGz = ht.dep.fit.bwd$origz, 
                                            n = 40, thresh = v.threshu, peakVal = peak.val.u, 
                                            nrep = 50000)[[1]]


# Convert from within cluster probability to over cluster probability
exceedVal <- 11; nlength <- 92
v.thresh <- u.thresh + (sig.u/xi)*((((1-v.threshu)/lambda.u)^(-xi))-1)
peak.val <- 39.9
peak.val.u <- EcdfWithGpd(data=peak.val,p=c(sig.u,xi,lambda.u),u=u.thresh)
sig.v <- sig.u+xi*(v.thresh-u.thresh)

# What is the value of P(X>eta|X>v_{1})
cond.prob.mv <- 1-pgpd(q = 39.9, loc = v.thresh, scale = sig.v, shape = xi)

# Get the over cluster probability for the original data sample
# Changing 1 to 2 leads to a change between N and N_{C}
Pi2ExceedProb(distr = piSM.orig[,1], exceed.val = exceedVal, p=c(sig.u,xi,lambda.u),
              n=nlength,u=u.thresh,v=v.thresh)[3]
Pi2ExceedProb(distr = piHT.orig[,1], exceed.val = exceedVal, p=c(sig.u,xi,lambda.u),
              n=nlength,u=u.thresh,v=v.thresh)[3]
Pi2ExceedProb(distr = piNP.orig[,1], exceed.val = exceedVal, p=c(sig.u,xi,lambda.u),
              n=nlength,u=u.thresh,v=v.thresh)[3]


# Now use peak value tail chain estimation to estimate P(N>11|M=eta)
# For the conditional extremes approach
fwd.bwd.HT <- FwdBwdSimWithIntensity(fwdpar = ht.dep.fit$par[1:2], fwdGz = ht.dep.fit$origz, 
                                     bwdpar = ht.dep.fit.bwd$par[1:2], bwdGz = ht.dep.fit.bwd$origz, 
                                     n = 40, thresh = v.threshu, peakVal = peak.val.u, nrep = 50000)

# For the non-parametric approach
fwd.bwd.NP <- FwdBwdSimWithIntensityNP(fwdpar = c(1,0), fwdGz = npsmGz, 
                                       bwdpar = c(1,0), bwdGz = npsmGz.bwd, 
                                       n = 40, thresh = v.threshu, peakVal = peak.val.u, nrep = 50000)

# For the parametric approach
fwd.bwd.SM <- FwdBwdSimWithIntensitySmith(fwdpar = gamma, fwdGz = NULL, 
                                          bwdpar = gamma, bwdGz = NULL, 
                                          n = 40, thresh = v.threshu, peakVal = peak.val.u, nrep = 50000)

# Estimate the probability of having more than 11 (consecutive) exceedances
# given the maximum value fixed at 39.9 degrees C
sum(fwd.bwd.HT[[1]][11:length(fwd.bwd.HT[[1]])])
sum(fwd.bwd.NP[[1]][11:length(fwd.bwd.NP[[1]][,1]),1])
sum(fwd.bwd.SM[[1]][11:length(fwd.bwd.SM[[1]][,1]),1])


# Use Monte-Carlo approaches to estimate P(N>11|M>eta)
num.ticks <- 10   # Have used 5 before but could try to use more
# Need to choose how to cover the space M>eta, either put interpolation
# points on fixed grid on uniform or Laplace scale
peak.vals.u <- plaplace(seq(qlaplace(peak.val.u),qlaplace(1-(1e-8)),length.out =num.ticks))
#peak.vals.u <- seq(peak.val.u,(1-(1e-8)),length.out = 5)
Pgt11cM <- numeric(length(peak.vals.u))
for ( i in 1:length(peak.vals.u) ){ 
  print(i)
  # Conditional extremes approach
  hold.list <- FwdBwdSimWithIntensity(fwdpar = ht.dep.fit$par[1:2], fwdGz = ht.dep.fit$origz, 
                                      bwdpar = ht.dep.fit.bwd$par[1:2], bwdGz = ht.dep.fit.bwd$origz, 
                                      n = 40, thresh = v.threshu, peakVal = peak.vals.u[i], nrep = 50000)
  # Non-parametric approach
  #hold.list <- FwdBwdSimWithIntensityNP(fwdpar = c(1,0), fwdGz = npsmGz, 
  #                                      bwdpar = c(1,0), bwdGz = npsmGz.bwd, 
  #                                      n = 40, thresh = v.threshu, peakVal = peak.val.u, nrep = 50000)
  # Parametric approach
  #hold.list <- FwdBwdSimWithIntensitySmith(fwdpar = gamma, fwdGz = NULL, 
  #                                         bwdpar = gamma, bwdGz = NULL, 
  #                                         n = 40, thresh = v.threshu, peakVal = peak.val.u, nrep = 50000)
  
  # Changing 1 to 2 leads to a change between N and N_{C}
  Pgt11cM[i] <- sum(hold.list[[2]][exceedVal:80])
}

# Calculate the probability of observing a set of 3 days with
# average value above the critical level
prob.exc.mort.event <- numeric(num.ticks)
for ( i in 2:num.ticks ){
  print(i)
  # Use the 'runmean' function from the 'caTools' package to calculate the
  # 3 day running means
  # (1) or (2) Use the lines below
  # NB: Set list number as 3 for conditional extremes and 2 for other approaches
  prob.exc.mort.event[i] <- sum(apply(X = eval(as.name(paste("fwd.bwd.HT.lss", i, sep = "")))[[3]] - 
                                        qlaplace(v.threshu), MARGIN = 2, 
                                      FUN = function(x){any(runmean(x, k = 3, alg = "fast")>0)}))/50000
  # (3) Use the lines below
  #prob.exc.mort.event[i] <- sum(apply(X = eval(as.name(paste("fwd.bwd.HT.lss", i, sep = "")))[[2]] - 
  #                                      (-log(-log(v.threshu))), MARGIN = 2, 
  #                                    FUN = function(x){any(runmean(x, k = 3, alg = "fast")>0)}))/50000
}

sig.pv <- sig.u+xi*(peak.val-u.thresh)
nsim <- 100000
# Simulate cluster maxima, if threshold not high enough then use theory
# from Eastoe and Tawn (2012) to obtain more realistic sample of 
# cluster maxima (not shown)
m.new <- rgpd(nsim,loc=peak.val,scale=sig.pv,shape=xi)
m.new.ugpd <- pgpd(q = m.new,loc=u.thresh,scale=sig.u,shape=xi)
m.new.u <- u.threshu + (1-u.threshu)*m.new.ugpd

# Linear interpolation for the duration longer than d days event
lo <- approx(x = peak.vals.u, y = Pgt11cM, xout = m.new.u[(m.new.u <= max(peak.vals.u))&(m.new.u > min(peak.vals.u))], 
             method = "linear")

# Linear interpolation for the 3-day average mortality style event
# lo <- approx(x = peak.vals.u, y = prob.exc.mort.event, xout = m.new.u[(m.new.u <= max(peak.vals.u))&(m.new.u > min(peak.vals.u))], 
#             method = "linear")

# Any values simulated above the upper end-point that simulations can be
# made to (since qlaplace(1)=Inf) are set to 'Pgt11cM[5]' and any values 
# that fall below are set to 'Pgt11cM[1]'
out.prob <- numeric(length(m.new.u))
out.prob[(m.new.u <= max(peak.vals.u))&(m.new.u > min(peak.vals.u))] <- lo$y
# For the d days event
out.prob[m.new.u > max(peak.vals.u)] <- Pgt11cM[num.ticks]
out.prob[m.new.u < min(peak.vals.u)] <- Pgt11cM[1]
# For the 3-day average mortality style event
# out.prob[m.new.u > max(peak.vals.u)] <- prob.exc.mort.event[num.ticks]
# out.prob[m.new.u < min(peak.vals.u)] <- prob.exc.mort.event[1]

# Take the median of the simulated sample
mean(x = out.prob)

# Use the output from the estimation of P(N>11|M=eta) to
# estimate the probability of at least one event in a year with
1 - exp(-PNgt11cM.probs[2]*mean(x = out.prob))

# Work out P(N>11,M>eta) by multiplying conditional probability above by P(M>eta)
cond.prob.mv*mean(x = out.prob)
1 - exp(-PNgt11cM.probs[2]*(cond.prob.mv*mean(x = out.prob)))
# Invert to find the return period associated with such a probability
1/((1 - exp(-PNgt11cM.probs[2]*(cond.prob.mv*mean(x = out.prob)))))
