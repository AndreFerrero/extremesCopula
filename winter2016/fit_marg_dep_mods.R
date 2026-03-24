
##############################################################
# Script that fits marginal GPD model to data set made up of #
# consecutive pairs and fits logistic dependence structure   #
# and conditional extremes approach to the resulting         #
# data series on uniform margins                             #
##############################################################

# Need to have data loaded in from 'load_data.R'
source(file="winter2016/marg_dep_mods_funcs.R"); source(file="winter2016/dur_dists_from_dep_funcs.R")

u.threshu <- 0.9   # Modelling threshold u on uniform scale 
# The return period (in years) associated with the desired critical level (in the paper 
# v_j corresponds to critical level associated with j year return period)
retPeriod <- 1; seas.length <- 92
v.threshu <- 1 - 1/(seas.length*retPeriod)   # Critical level v on uniform scale

# The modelling threshold on the original margins
u.thresh <- quantile(JJA.data$TX*0.1, probs = u.threshu)

# Fit the GPD distribution to the margins and obtain logistic dependence parameter
dat <- j.data; marg.gpd.fit <- FitGpd(dat = dat, u = u.thresh)

# Store the parameter outputs from the optim call. Gamma is the logistic dependence parameter
# and (sig.u,xi,lambda.u) are the scale, shape and rate parameters respectively             
gamma <- marg.gpd.fit$par[1]; sig.u <- marg.gpd.fit$par[2]
xi <- marg.gpd.fit$par[3]; lambda.u <- marg.gpd.fit$par[4]

# To obtain dependence parameters for semi-parametric conditional extremes approach 
# need to transform data onto uniform margins by fitting eqn (10) of my paper above 
# u and the empirical cdf below u     
# The work around added below means that there won't be a problem with slightly 
# different ECDF functions being fitted below the threshold
# Remove the missing values before using GPD with ECDF to transform values 
udat1 <- EcdfWithGpd(data=(JJA.data$TX*0.1)[JJA.data$TX>-100],p=c(sig.u,xi,lambda.u),
                     u=u.thresh)
# Need to re-paste in the missing values before the lagging is done
udat2 <- JJA.data$TX*0.1
# Paste the unit values onto the original array
udat2[udat2>-100] <- udat1
# Create the lagged data and remove any overlapping values between years
udat <- GetTimeLagkData(dat=udat2,k=1,seas.length=92)
# Remember to remove the missing values again!
udat <- udat[apply(udat>-100,MARGIN=1,FUN=all),]

# The 1-year return level is 34.97 on the original scale
u.thresh+(sig.u/xi)*(((1-v.threshu)/(lambda.u))^(-xi)-1) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# QQ-plot for Figure 2 in Section 5.1
# Going to create the plot manually
# The value of the modelling threshold on the original
# scale for the plot
u.thresh.plot <- u.thresh+(sig.u/xi)*(((1-u.threshu)/(lambda.u))^(-xi)-1)
# Obtain 95% confidence bands around each of the values 
kval <- 1.36/sqrt(sum(dat[,1]>u.thresh.plot))  # Obtain the critical value of K-S test statistic
ecdf.fun <- ecdf(dat[dat[,1]>u.thresh.plot,1])   # Fit the ECDF to the set of exceedances
emp.dat.u <- ecdf.fun(dat[dat[,1]>u.thresh.plot,1])   # Get the exceedance values on uniform
# Calculate the intervals around the value of F_{n}(x) that are used
# to obtain the intervals
qqplot.int.u <- sapply(X = sort(emp.dat.u), FUN = function(x,k){c(x-k,x+k)}, k = kval)
# Any values that go outside the range of the data need to pulled back inside
# as we are using an empirical measure so can't go outside the bounds of the data
qqplot.int.u[qqplot.int.u < 0] <- 0; qqplot.int.u[qqplot.int.u > 1] <- 1
qqplot.int.o <- apply(X = qqplot.int.u, MARGIN = 1, FUN = quantile, x = dat[dat[,1]>u.thresh.plot,1])
# Do not want the line stuck at 40.3 for the plot as this
# is not done in other example of this type of plot
set.to.change <- qqplot.int.o == max(qqplot.int.o)
set.to.change[593,2] <- FALSE
qqplot.int.o[set.to.change] <- NA

# Not doing this so need to obtain quantiles from observed data
# and model values
par(mfrow=c(1,1), mar = c(5,5,2,2))
# Plot the points from empirical and model against one another 
plot(x = gpdq(a = c(sig.u,xi), u = u.thresh.plot, p = (1 - (1:length(dat[dat[,1]>u.thresh.plot,1])/(length(dat[dat[,1]>u.thresh.plot,1])+1)))),
     y = sort(dat[dat[,1]>u.thresh.plot,1]), xlim = c(u.thresh.plot, 40.5), ylim = c(u.thresh.plot, 40.5),
     xlab = "Model", ylab = "Empirical", cex.lab = 1.2, cex.axis = 1.2)
# Add a diagonal line which would signify perfect matching 
abline(0, 1, col = 4, lwd = 2)
# Add the 95% confidence intervals to the plot 
lines(x = gpdq(a = c(sig.u,xi), u = u.thresh.plot, p = (1 - (1:length(dat[dat[,1]>u.thresh.plot,1])/(length(dat[dat[,1]>u.thresh.plot,1])+1)))),
      y = qqplot.int.o[,1], lty = 3, lwd = 2)
lines(x = gpdq(a = c(sig.u,xi), u = u.thresh.plot, p = (1 - (1:length(dat[dat[,1]>u.thresh.plot,1])/(length(dat[dat[,1]>u.thresh.plot,1])+1)))),
      y = qqplot.int.o[,2], lty = 3, lwd = 2)
box(lwd = 3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Use hessian to obtain 95% CIs for each of the parameters #
# Fit the penalised version of HT
ht.dep.fit <- BveHTDepPen(data = udat, mod.thresh.u = u.threshu, crit.lev.u = v.threshu, nsim = 10000,
                          sim.exc = rexp(10000), bwd = 0 ,lam.pen = 30)
ht.dep.fit$par

# For the peak value simulation method will need to fit the H+T approach to the reversed data
ht.dep.fit.bwd <- BveHTDepPen(data = udat[,c(2,1)], mod.thresh.u = u.threshu, crit.lev.u = v.threshu,
                              nsim = 10000, sim.exc = rexp(10000), bwd = 0, lam.pen = 30)
ht.dep.fit.bwd$par

# Need to work out the empirical differences for the non-parametric method
udat.abovexthresh <- udat[udat[,1]>u.threshu,]
npsmGz <- qlaplace(udat.abovexthresh[,2])-qlaplace(udat.abovexthresh[,1])

udat.abovexthresh.bwd <- udat[udat[,2]>u.threshu,]
npsmGz.bwd <- qlaplace(udat.abovexthresh.bwd[,1])-qlaplace(udat.abovexthresh.bwd[,2])
