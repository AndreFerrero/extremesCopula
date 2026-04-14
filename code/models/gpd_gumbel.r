
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
