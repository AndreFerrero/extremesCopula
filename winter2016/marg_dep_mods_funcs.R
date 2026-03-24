
######################################################
# Set of functions required to estimate marginal and #
# dependence structure features                      #
######################################################

# Load in packages required
library(VGAM); library(extRemes); library(abind); library(Iso); library(evd);
library(Hmisc); library(msm); library(MASS); library(zoo); library(caTools); library(ismev)

EcdfWithGpd <- function(data,p,u){
  # Transform data from original margins onto uniform margins using GPD above 
  # threshold whilst using empirical cdf to points that lie below the threshold
  #
  # Args:
  #   data: Single array of data given on original margins
  #   p: Vector of parameters (sigma,xi,lambda)
  #   u: Threshold on original margins above which to fit GPD
  #
  # Returns:
  #   Array of data values on uniform margins
  #
  # Error handling:
  if (length(p)!=3){
    stop("Length of parameter vector incorrect")
  }
  # Function START
  fit.cdf <- ecdf(data)
  cdfs <- as.vector(sapply(data, fit.cdf))
  Fex <- 1 - p[3]*pmax(0,(1+(p[2]*((data[data>u]-u)/p[1]))))^(-1/p[2])   # pmax function gives the parallel maximum (element by element maximum)
  cdfs[data>u] <- Fex
  return(cdfs)
}


FitGpd <- function(dat,u,start = c(0.1,0.1,0.1,0.1), optim.type=1){
  # Function that wraps the fitting of the logistic 
  # dependence structure and marginal GPD distribution
  # into a function that can be run more easily.
  # If a regular GPD fit to univariate data is required just use
  # 'gpd.fit' from 'ismev' package
  #
  # Have bounded the lowest value that the shape parameter can take
  # at -0.5 to ensure no problems with the numerical minimisation
  # of the negative log-likelihood function.
  #
  # Args:
  #   dat: 2-column data array containing bivariate data
  #   u: Threshold on original margins above which to fit GPD
  #   optim.type: Choose 1 for "L-BFGS-B" optimization and 2 for "Nelder-Mead"
  #
  # Returns:
  #   The parameter values that maximise the log-likelihood with other
  #   diagnostics of the fit
  #
  # Function START
  
  # Bounds and starting value for L-BFGS-B style optimization #
  bl <- c(0.001,0.001,-0.5,0.001)
  bu <- c(1,Inf,Inf,0.9999)     
  
  if ( optim.type==1 ){
    marg.gpd.fit <- optim(start,fn=PotNllhGpd,method="L-BFGS-B",data=dat,u=u,
                          lower=bl,upper=bu,hessian=TRUE)
  } else {
    marg.gpd.fit <- optim(start,fn=PotNllhGpd,method="Nelder-Mead",data=dat,u=u,
                          hessian=TRUE)
  }
  
  return( marg.gpd.fit )
  
}


PotNllhGpd <- function(p,data,u){
  # Compute the logistic dependence parameter gamma and marginal parameters (sigma,xi,lambda)
  # Use equation (3.2) in Ledford and Tawn (1996) to fit the marginal parameters
  # Function to be called inside 'optim' function in R.
  # Code below works in the stationary scenario at the moment.
  #
  # Args:
  #   p: Vector of parameters (gamma,sigma,xi,lambda) to be estimated
  #   data: Two column data frame with values given on original margins
  #   u: Threshold on original margin at which to fit model
  #
  # Returns:
  #   The negative log-likelihood
  #
  # Error handling:
  if (length(p)!=4){
    stop("Length of parameter vector incorrect")
  }
  if (dim(data)[2]!=2){
    stop("Need bivariate data")
  }
  
  # Function START
  y1 <- data[,1][data[,1]>u]
  y2 <- data[,2][data[,2]>u]
  
  n00 <- sum(data[,1]<u & data[,2]<u)   # How many data points fall in bottom left quadrant?
  
  # Look at data in terms of exceedances => GPD 
  # The steps of fitting a GPD are separated below 
  
  sy1 <- 1 + p[3]*(y1-u)/p[2]
  sy2 <- 1 + p[3]*(y2-u)/p[2]
  
  ty1 <- pmax(sy1,0)^(-1/p[3])
  ty2 <- pmax(sy2,0)^(-1/p[3])
  
  z1 <- -1/log(1-p[4]*ty1)
  z2 <- -1/log(1-p[4]*ty2)
  
  r <- -1/log(1-p[4])
  
  # Split into different regions 
  
  y10 <- data[,1][data[,1]>u & data[,2]<u]
  y01 <- data[,2][data[,2]>u & data[,1]<u]
  y111 <- data[,1][data[,1]>u & data[,2]>u]
  y112 <- data[,2][data[,1]>u & data[,2]>u]
  
  sy10 <- 1 + p[3]*(y10-u)/p[2]
  sy01 <- 1 + p[3]*(y01-u)/p[2]
  sy111 <- 1 + p[3]*(y111-u)/p[2]
  sy112 <- 1 + p[3]*(y112-u)/p[2]
  
  ty10 <- pmax(sy10,1e-4)^(-1/p[3])
  ty01 <- pmax(sy01,1e-4)^(-1/p[3])
  ty111 <- pmax(sy111,1e-4)^(-1/p[3])
  ty112 <- pmax(sy112,1e-4)^(-1/p[3])
  
  z10 <- -1/log(1-p[4]*ty10)
  z01 <- -1/log(1-p[4]*ty01)
  z111 <- -1/log(1-p[4]*ty111)
  z112 <- -1/log(1-p[4]*ty112)
  
  K10 <- (p[4])*(1/p[2])*(ty10^(1+p[3]))*(z10^2)*exp(1/z10)
  K01 <- (p[4])*(1/p[2])*(ty01^(1+p[3]))*(z01^2)*exp(1/z01)
  K111 <- (p[4])*(1/p[2])*(ty111^(1+p[3]))*(z111^2)*exp(1/z111)
  K112 <- (p[4])*(1/p[2])*(ty112^(1+p[3]))*(z112^2)*exp(1/z112)
  
  # Calculate the likelihood from points falling in each quadrant where:
  # 00 => bottom left, 01 => top left, 10 => bottom right, 11 => top right
  L00 <- ( exp(-(r^(-1/p[1])+r^(-1/p[1]))^(p[1])) )
  L01 <- ( ((r^(-1/p[1])+z01^(-1/p[1]))^(p[1]-1))*(z01^(-1/p[1]-1))*exp(-(r^(-1/p[1])+z01^(-1/p[1]))^(p[1]))*(K01))
  L10 <- ( ((z10^(-1/p[1])+r^(-1/p[1]))^(p[1]-1))*(z10^(-1/p[1]-1))*exp(-(z10^(-1/p[1])+r^(-1/p[1]))^(p[1]))*(K10) )
  L11 <- ( (((z111^(-1/p[1])+z112^(-1/p[1]))^(p[1]-1))*(z111^(-1/p[1]-1))*((z111^(-1/p[1])+z112^(-1/p[1]))^(p[1]-1))*(z112^(-1/p[1]-1))- ((p[1]-1)/p[1])*((z111^(-1/p[1])+z112^(-1/p[1]))^(p[1]-2))*((z111*z112)^(-1/p[1]-1)))*exp(-(z111^(-1/p[1])+z112^(-1/p[1]))^(p[1])) * (K111*K112) )
  
  vec.test <- 1+p[3]*(y2-u)/(p[2])
  
  ML <- ( (1/p[2])*(p[4])*((pmax(vec.test,1e-4))^(-1/p[3]-1)) )
  
  nexy <- length(data[,2])-length(y2)
  
  if (any(is.na(L00))==TRUE || any(is.na(L01))==TRUE || any(is.na(L10))==TRUE ||any(is.na(L11))==TRUE || any(is.na(ML))==TRUE){
    return(10^6)
  } else {
    return( - ( n00*log(L00) + sum(log(L01)) + sum(log(L10)) + sum(log(L11)) - nexy*log(1-p[4]) - sum(log(ML)) )  )
  }
}


GetTimeLagkData <- function(dat,k,seas.length=90,overlap=TRUE){
  # Function to obtain consecutive values for time lag k.
  # Is a general version of the function 'getTimeLag2Data'
  #
  # Args:
  #   dat: Data from which we wish to obtain values for t,t+1,...,t+k
  #   k: Time lag
  #   seas.length: If there is some overlap from year to year where there shouldn't 
  #                be (i.e. taking only summer months)
  #   overlap: Is there an overlap from year to year?
  #
  # Returns:
  #   A nx(k+1) array that contains the original data shifted to different 
  #   time lags 0,1,2,...,k
  #
  # Function START
  
  consec.dat <- array(NA,dim=c(length(dat)+k,k+1))
  for ( i in 1:(k+1) ){
    # Add NAs to the top and bottom of the data to allow time-lag shift
    consec.dat[,i] <- c(rep(NA,length.out=(k-i+1)),dat,rep(NA,length.out=(i-1)))
  }
  if (overlap==TRUE){   # Is there an overlap between years which is not natural?
    if ( length(seas.length) == 1){
      overlap.vals.list <- list()
      for ( j in 1:k ){
        # Take the points where the overlap occurs
        overlap.vals.list[[j]] <- ((2:(dim(consec.dat)[1]/seas.length)-1)*seas.length+j) 
      }
      overlap.vals <- unlist(x=overlap.vals.list)
      # Give the final data with the overlap and artifically inserted NAs removed
      out.dat <- consec.dat[-c(1:k,overlap.vals,dim(consec.dat)[1]-((k-1):0)),]
    } else {   # If different years have a different season length
      if ( sum(seas.length)!=length(dat) ){   # Need the sum of days in 'seas.length' to be equal to the length of the data set
        stop( "Length of seasons must be equal to length of original data set" )
      }
      overlap.vals.list <- list()
      for ( j in 1:k ){
        # Take the points where the overlap occurs
        overlap.vals.list[[j]] <- (cumsum(seas.length)+j)[-length(cumsum(seas.length)+j)] 
      }
      overlap.vals <- unlist(x=overlap.vals.list)
      # Give the final data with the overlap and artifically inserted NAs removed
      out.dat <- consec.dat[-c(1:k,overlap.vals,dim(consec.dat)[1]-((k-1):0)),]
    }
  } else {
    # Even if no overlap need to remove artificially inserted NAs
    out.dat <- consec.dat[-c(1:k,dim(consec.dat)[1]-((k-1):0)),]  
  }
  return( out.dat )
}


BveHTDepPen <- function(data, mod.thresh.u, crit.lev.u, nsim, sim.exc = rexp(1000), bwd = 0, lam.pen = 0, 
                         shrink = TRUE, zSampKnown = FALSE, noiseKnown = FALSE){
  # Compute the dependence parameters (alpha,beta) for the conditional extremes approach 
  # and output new simulated sample. Has additional penalty function to ensure that 
  # dependence isn't lost to the parameters of the Normal distribution taken as a 
  # false working assumption
  #
  # Args:
  #   data: Two column data frame with values given on uniform margins
  #   mod.thresh.u: Modelling threshold on uniform margins
  #   crit.lev.u: Critical level on uniform margins (now set as equal in x and y directions)
  #   nsim: Number of new data points to be simulated above the critical level
  #   bwd: Choice of bandwidth if kernal smoothing is added
  #   lam.pen: Penalty value to ensure that dependence is not lost to discarded Normal parameters
  #   shrink: Is kernal shrinkage to be included
  #   zSampKnown: If the z-sample to be picked is known then set this as a specific value instead of FALSE
  #   noiseKnown: If the kernal noise added needs to be specified then set this as a specific value instead of FALSE
  #   cond.list: The list that contains all the output
  #
  # Returns:
  #   Numerous outputs added to a list, including the dependence parameter values, full
  #   sample z and ones picked to simulate new exceedances.
  #
  # Function START
  
  x.dat.l <- qlaplace(data[, 1])
  y.dat.l <- qlaplace(data[, 2])
  dat.l <- cbind(x.dat.l, y.dat.l)
  thresh.l <- qlaplace(mod.thresh.u)
  dat.lu <- dat.l[dat.l[, 1] > thresh.l, ]
  x.dat.lu <- dat.lu[, 1]
  y.dat.lu <- dat.lu[, 2]
  cond.llik <- function(p, x, y) {
    return((length(x)/2) * (log(2 * pi)) + (length(x)/2) * log(p[4]) + 
             p[2] * sum(log(x)) + (1/((2 * p[4]))) * sum(((y - 
                                                             p[1] * x - (p[3] * x^(p[2])))^2)/(x^(2 * p[2]))) + lam.pen*(p[3]^2))
  }
  bl <- c(-1, -10, -10^6, 1e-04)
  bu <- c(1, 1, 10^6, 10^6)
  start <- c(0.1, 0.1, 0, 1)
  par.vec <- optim(start, cond.llik, method = "L-BFGS-B", x = x.dat.lu, 
                   y = y.dat.lu, lower = bl, upper = bu, hessian = TRUE)
  nw.par <- par.vec$par[1:4]
  nw.nllh <- par.vec$value
  z.samp <- (y.dat.lu - nw.par[1] * x.dat.lu)/(x.dat.lu^(nw.par[2]))
  nw.musig <- c(mean(z.samp), var(z.samp))
  h.thresh.xl <- qlaplace(crit.lev.u)
  h.thresh.yl <- qlaplace(crit.lev.u)
  xn.l <- h.thresh.xl + sim.exc
  if (zSampKnown == FALSE) {
    zn.ls <- sample(z.samp, nsim, replace = T)
  } else {
    zn.ls <- zSampKnown
  }
  if (bwd == "scott") {
    h <- 1.06 * (min(sd(zn.ls), (IQR(zn.ls)/1.34))) * (nsim^(-1/5))
  } else {
    h <- bwd
  }
  if (shrink == TRUE) {
    a <- sqrt(1 - (h^2))
    zn.lm <- mean(zn.ls)
    zn.lsd <- sd(zn.ls)
    zn.bsd <- h * zn.lsd
    if (noiseKnown == FALSE) {
      gNoise <- rnorm(nsim, mean = 0, sd = zn.bsd)
    } else {
      gNoise <- noiseKnown
    }
    zn.l <- a * zn.ls + (1 - a) * zn.lm + gNoise
  } else {
    zn.lm <- mean(zn.ls)
    zn.lsd <- sd(zn.ls)
    zn.bsd <- h * zn.lsd
    if (noiseKnown == FALSE) {
      gNoise <- rnorm(nsim, mean = 0, sd = zn.bsd)
    } else {
      gNoise <- noiseKnown
    }
    zn.l <- zn.ls + gNoise
  }
  yn.l <- nw.par[1] * xn.l + (xn.l^(nw.par[2])) * zn.l
  pAcX <- sum(yn.l > h.thresh.yl)/length(yn.l)
  pX <- 1 - crit.lev.u
  pA <- pAcX * pX
  
  cond.list <- list()
  cond.list$diag <- c(pX, pAcX, pA)
  cond.list$x <- plaplace(xn.l)
  cond.list$y <- plaplace(yn.l)
  cond.list$z <- plaplace(zn.l)
  cond.list$zsamp <- zn.ls
  cond.list$gNoise <- gNoise
  cond.list$origz <- plaplace(z.samp)
  cond.list$par <- nw.par
  cond.list$nllh <- nw.nllh
  cond.list$hessian <- par.vec$hessian
  
  return( cond.list )
  
}
