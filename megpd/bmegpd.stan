// megpd_biv.stan
functions {
  // Log-PDF of the EGPD distribution (Radial component)
  real egpd_lpdf(real r, real kappa, real sigma, real xi) {
    real H;
    real h;
    if (abs(xi) < 1e-10) {
      H = 1 - exp(-r / sigma);
      h = (1/sigma) * exp(-r / sigma);
    } else {
      real z = r / sigma;
      H = fmax(0, 1 - pow(1 + xi * z, -1/xi));
      h = (1/sigma) * pow(1 + xi * z, -1/xi - 1);
    }
    return log(kappa) + log(h) + (kappa - 1) * log(H);
  }
}

data {
  int<lower=0> N;          // Number of observations
  vector[N] x;             // Variable 1
  vector[N] y;             // Variable 2
  int<lower=0> K;          // Number of spline basis functions
  matrix[N, K] B;          // Pre-computed B-spline matrix
}

transformed data {
  vector[N] R = x + y;     // Radial component
  vector[N] V = log(y ./ x); // Angular component (log-ratio)
}

parameters {
  // EGPD Parameters
  real<lower=0> kappa;
  real<lower=0> sigma;
  real<lower=0, upper=1.0> xi; // Constraint for heavy-tailed stability
  
  // Angular Spline Parameters
  vector[K] beta;          // Spline coefficients
  real<lower=0> tau;       // Smoothing parameter (SD of Random Walk)
}

transformed parameters {
  // delta is the heteroscedastic scale of the log-ratio
  // We model log(delta) to ensure positivity
  vector[N] delta = exp(B * beta) + 0.01; 
}

model {
  // 1. Priors for EGPD
  kappa ~ lognormal(0.5, 0.5);
  sigma ~ lognormal(0.5, 0.5);
  xi ~ normal(0.2, 0.2);
  
  // 2. Bayesian P-Spline Prior (Random Walk 1)
  tau ~ exponential(1);
  beta[1] ~ normal(0, 1);
  for (k in 2:K) {
    beta[k] ~ normal(beta[k-1], tau);
  }
  
  // 3. Likelihood: f(x, y) = f_R(R) * f_V(V|R) * |Jacobian|
  // Jacobian of (x,y) -> (R, V) is R/(x*y)
  for (n in 1:N) {
    // Radial Likelihood
    target += egpd_lpdf(R[n] | kappa, sigma, xi);
    
    // Angular Likelihood (Gaussian log-ratio)
    target += normal_lpdf(V[n] | 0, delta[n]);
    
    // Jacobian Correction
    target += log(R[n]) - log(x[n]) - log(y[n]);
  }
}
