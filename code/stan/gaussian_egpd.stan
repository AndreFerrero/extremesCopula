functions {
  // --- GPD & EGPD functions
  real gpd_lpdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-5) return -y - log(sigma);
    if (1 + xi * y <= 0) return -1e10;
    return -(1 + 1/xi) * log1p(xi * y) - log(sigma);
  }

  real gpd_lcdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-5) return log1m_exp(-y);
    if (1 + xi * y <= 0) return -1e10;
    return log1m_exp(-(1/xi) * log1p(xi * y));
  }

  real egpd_lpdf(real x, real kappa, real sigma, real xi) {
    real log_G = gpd_lcdf(x | sigma, xi);
    real log_g = gpd_lpdf(x | sigma, xi);
    return log(kappa) + (kappa - 1) * log_G + log_g;
  }

  real egpd_rng(real kappa, real sigma, real xi) {
    real p_gpd = pow(uniform_rng(0, 1), 1.0/kappa);
    if (abs(xi) < 1e-5) return -sigma * log1m(p_gpd);
    return (sigma / xi) * (pow(1 - p_gpd, -xi) - 1);
  }

  // --- Gaussian Copula Density ---
  real gaussian_copula_lpdf(real u, real v, real rho) {
    real x = inv_Phi(u);
    real y = inv_Phi(v);
    real rho2 = square(rho);
    // Standard bivariate normal copula log-density
    return -0.5 * log1m(rho2) - (rho2 * (square(x) + square(y)) - 2 * rho * x * y) / (2 * (1 - rho2));
  }
}

data {
  int<lower=2> T;
  vector<lower=0>[T] x;
  int<lower=0, upper=1> prior_check;
  int<lower=0, upper=1> run_ppc; 
}

parameters {
  real<lower=0, upper=min(x)> mu;
  real<lower=0.01> kappa;
  real<lower=0.01> sigma;
  real<lower=0, upper=0.5> xi; 
  real<lower=-0.99, upper=0.99> rho; // Correlation parameter
}

model {
  // Priors (Aligned with your weakly-informative setup)
  mu ~ normal(5, 2.5); 
  kappa ~ lognormal(2, 1);
  sigma ~ exponential(0.1);
  xi ~ gamma(2, 10);
  rho ~ normal(0, 0.5); // Prior on rho

  if (prior_check == 0) {
    // Likelihood: First observation
    target += egpd_lpdf(x[1] - mu | kappa, sigma, xi);

    // Markov loop
    for (t in 2:T) {
      real u = exp(gpd_lcdf(x[t] - mu | sigma, xi) * kappa);
      real v = exp(gpd_lcdf(x[t-1] - mu | sigma, xi) * kappa);

      target += egpd_lpdf(x[t] - mu | kappa, sigma, xi);
      target += gaussian_copula_lpdf(u| v, rho);
    }
  }
}

generated quantities {
  vector[T] x_rep;
  vector[T] log_lik;

  // 1. Pointwise log-likelihood for LOO-CV
  log_lik[1] = egpd_lpdf(x[1] - mu | kappa, sigma, xi);
  for (t in 2:T) {
    real u = exp(gpd_lcdf(x[t] - mu | sigma, xi) * kappa);
    real v = exp(gpd_lcdf(x[t-1] - mu | sigma, xi) * kappa);
    log_lik[t] = egpd_lpdf(x[t] - mu | kappa, sigma, xi) + gaussian_copula_lpdf(u| v, rho);
  }

  // 2. Fast Analytical PPC
  if (run_ppc == 1) {
    x_rep[1] = mu + egpd_rng(kappa, sigma, xi);
    for (t in 2:T) {
      real v_prev = exp(gpd_lcdf(x_rep[t-1] - mu | sigma, xi) * kappa);
      real w = uniform_rng(0, 1);
      
      // Inverse Rosenblatt for Gaussian is analytical:
      // u_next = Phi( rho * Phi^-1(v) + sqrt(1-rho^2) * Phi^-1(w) )
      real u_next = Phi(rho * inv_Phi(v_prev) + sqrt(1 - square(rho)) * inv_Phi(w));
      
      real g_inv_p = pow(u_next, 1.0/kappa);
      real excess = (abs(xi) < 1e-5) ? -sigma * log1m(g_inv_p) : (sigma / xi) * (pow(1 - g_inv_p, -xi) - 1);
      x_rep[t] = mu + excess;
    }
  }
}
