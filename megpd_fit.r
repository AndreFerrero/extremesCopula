# Load required libraries
library(mgcv)
library(egpd) # Ensure your specific EGPD package is loaded

fit_megpd_bivariate <- function(x, y, m_degree = NULL) {
  
  # 0. Data Preparation
  # The model is defined for positive intensities (e.g., river discharge)
  n <- length(x)
  R <- x + y               # Radial component (sum-norm)
  V <- log(x / y)          # Angular component (log-ratio)
  V <- V - mean(V) # Centering V to have mean 0 (as per the Z ~ zero-mean assumption)

  cat("Step 1: Fitting Radial EGPD via Bernstein Polynomials...\n")
  
  # 1. First Step: Radial component EGPD
  # The paper suggests m = floor(0.5 * n / log(n)) for Bernstein degree
  if (is.null(m_degree)) {
    m_degree <- floor(0.5 * n / log(n))
  }
  
  # Note: The 'egpd' package usually has a fit function. 
  # We specify the 'bernstein' model as described in Section 3.
  # If using the standard 'egpd' package:
  fit_r <- egpd::fitegpd(
     R,
     type = 1,
     family = "egpd",
     method = "bernstein",
     bernstein.m = m_degree
)
  
  cat("Step 2: Fitting Heteroscedastic Angular Component via mgcv...\n")
  
  # 2. Second Step: Logistic-heteroscedastic modeling
  # We model V | R=r ~ N(0, delta(r)^2)
  # According to Eq (19), log(delta(r)) is a linear combination of basis functions.
  # We use mgcv's 'gauslss' family to model the scale (sigma) specifically.
  
  # In 'gauslss', the first formula is for the mean (mu), the second for log(sigma).
  # The paper assumes a zero-mean exchangeable vector Z (Eq 11). 
  # Therefore, we fix intercept to 0 and no covariates for the mean.
  
  # K = 10 is the default basis dimension suggested in the paper (Section 4.1).
  fit_delta <- gam(list(
    V ~  1,              # Mean fixed at 0 (as per the Z ~ zero-mean assumption)
    ~ s(R, k = 10, bs = "cr") # log(sigma) modeled as a cubic spline of R
  ),
  family = gaulss(),
  data = data.frame(V = V, R = R))
  
  # 3. Return Results
  results <- list(
    radial_fit = fit_r,
    angular_fit = fit_delta
  )
  
  return(results)
}
