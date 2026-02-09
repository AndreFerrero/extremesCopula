functions {
  // ... [Keep your gpd_lpdf, gpd_lcdf, egpd_lpdf, and gumbel_copula_lpdf functions here] ...
  
  real gpd_lpdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-10) return -y - log(sigma);
    else {
      if (1 + xi * y <= 0) return -1e10;
      return -(1 + 1/xi) * log1p(xi * y) - log(sigma);
    }
  }

  real gpd_lcdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-10) return log1m_exp(-y);
    else return log1m(pow(1 + xi * y, -1.0/xi));
  }

  real egpd_lpdf(real x, real kappa, real sigma, real xi) {
    real log_G = gpd_lcdf(x | sigma, xi);
    real log_g = gpd_lpdf(x | sigma, xi);
    return log(kappa) + (kappa - 1) * log_G + log_g;
  }

  real gumbel_copula_lpdf(real u, real v, real theta) {
    if (u <= 1e-10 || u >= 1.0-1e-10 || v <= 1e-10 || v >= 1.0-1e-10) return -1e10;
    real w_u = -log(u);
    real w_v = -log(v);
    real r = pow(pow(w_u, theta) + pow(w_v, theta), 1.0/theta);
    real term_top = log(r + theta - 1.0);
    return -r + (theta - 1.0) * (log(w_u) + log(w_v)) + (1.0 - 2.0 * theta) * log(r) + term_top - (log(u) + log(v));
  }
}

data {
  int<lower=2> T;
  vector<lower=0>[T] x;
  vector[T] z; // The regressor (e.g., time or temperature)
}

parameters {
  real<lower=0.01> kappa;
  real<lower=-0.5, upper=0.5> xi;
  
  // Regressors for Sigma (log-link)
  real beta0_sigma; 
  real beta1_sigma;
  
  // Regressors for Theta (shifted log-link to ensure theta > 1)
  real beta0_theta;
  real beta1_theta;
}

transformed parameters {
  vector[T] sigma;
  vector[T] theta;
  
  for (t in 1:T) {
    sigma[t] = exp(beta0_sigma + beta1_sigma * z[t]);
    theta[t] = 1.0 + exp(beta0_theta + beta1_theta * z[t]);
  }
}

model {
  // Priors on coefficients (Weakly informative)
  kappa ~ lognormal(0, 1);
  xi ~ normal(0, 0.2);
  
  beta0_sigma ~ normal(0, 1);
  beta1_sigma ~ normal(0, 0.5); // Slopes should be tighter
  
  beta0_theta ~ normal(0, 1);
  beta1_theta ~ normal(0, 0.5);

  // Likelihood
  target += egpd_lpdf(x[1] | kappa, sigma[1], xi);
  for (t in 2:T) {
    real u = exp(gpd_lcdf(x[t]   | sigma[t], xi) * kappa);
    real v = exp(gpd_lcdf(x[t-1] | sigma[t-1], xi) * kappa);
    
    target += egpd_lpdf(x[t] | kappa, sigma[t], xi);
    target += gumbel_copula_lpdf(u| v , theta[t]);
  }
}
