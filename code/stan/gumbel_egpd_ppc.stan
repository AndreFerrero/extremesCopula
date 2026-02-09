functions {
  // Standard GPD Log-PDF
  real gpd_lpdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-10) {
      return -y - log(sigma);
    } else {
      if (1 + xi * y <= 0) return -1e10;
      return -(1 + 1/xi) * log1p(xi * y) - log(sigma);
    }
  }

  // Standard GPD Log-CDF
  real gpd_lcdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-10) {
      return log1m_exp(-y);
    } else {
      return log1m(pow(1 + xi * y, -1.0/xi));
    }
  }
  
  // EGPD Quantile Function for Simulation
  real egpd_rng(real kappa, real sigma, real xi) {
    real p = pow(uniform_rng(0, 1), 1.0/kappa);
    if (abs(xi) < 1e-10) return -sigma * log1m(p);
    else return (sigma / xi) * (pow(1.0 - p, -xi) - 1.0);
  }

  // EGPD Log-PDF
  real egpd_lpdf(real x, real kappa, real sigma, real xi) {
    real log_G = gpd_lcdf(x | sigma, xi);
    real log_g = gpd_lpdf(x | sigma, xi);
    return log(kappa) + (kappa - 1) * log_G + log_g;
  }

  // Stable Gumbel Copula Log-Density (d=2)
  real gumbel_copula_lpdf(real u, real v, real theta) {
    if (u <= 1e-10 || u >= 1.0-1e-10 || v <= 1e-10 || v >= 1.0-1e-10) return -1e10;
    real w_u = -log(u);
    real w_v = -log(v);
    real r = pow(pow(w_u, theta) + pow(w_v, theta), 1.0/theta);
    real term_top = log(r + theta - 1.0);
    return -r + (theta - 1.0) * (log(w_u) + log(w_v)) + (1.0 - 2.0 * theta) * log(r) + term_top - (log(u) + log(v));
  }
  
  // Conditional Gumbel CDF (h-function) for simulation
  real gumbel_hfunc(real u, real v, real theta) {
    real ln_u = -log(u);
    real ln_v = -log(v);
    real term = pow(ln_u, theta) + pow(ln_v, theta);
    real C = exp(-pow(term, 1.0/theta));
    return C * pow(ln_v, theta - 1.0) * pow(term, 1.0/theta - 1.0) / v;
  }
}

data {
  int<lower=2> T;
  vector<lower=0>[T] x;
  int<lower=0, upper=1> prior_check; // 1 = Skip Likelihood, 0 = Use Likelihood
}

parameters {
  real<lower=0.01> kappa;
  real<lower=0.01> sigma;
  real<lower=-0.5, upper=0.5> xi; 
  real<lower=0> thetam1; 
}

transformed parameters {
  real theta = thetam1 + 1.0;
}

model {
  // --- Priors ---
  kappa ~ uniform(0.1, 10);
  sigma ~ uniform(0.1, 15);
  xi ~ uniform(-0.2, 0.5);
  thetam1 ~ uniform(0, 10);

  // --- Likelihood (Conditional on prior_check) ---
  if (prior_check == 0) {
    target += egpd_lpdf(x[1] | kappa, sigma, xi);
    for (t in 2:T) {
      real u = exp(gpd_lcdf(x[t]   | sigma, xi) * kappa);
      real v = exp(gpd_lcdf(x[t-1] | sigma, xi) * kappa);
      target += egpd_lpdf(x[t] | kappa, sigma, xi);
      target += gumbel_copula_lpdf(u| v, theta);
    }
  }
}

generated quantities {
  vector[T] x_rep;
  x_rep[1] = egpd_rng(kappa, sigma, xi);
  
  // Rosenblat Transform
  for (t in 2:T) {
    real v_prev = exp(gpd_lcdf(x_rep[t-1] | sigma, xi) * kappa);
    real w = uniform_rng(0, 1);
    
    // Bisection solver
    real low = 0.00001;
    real high = 0.99999;
    for (i in 1:16) {
      real mid = (low + high) / 2.0;
      if (gumbel_hfunc(mid, v_prev, theta) < w) low = mid;
      else high = mid;
    }
    real u_next = (low + high) / 2.0;
    
    // Inverse transformation
    real g_inv_p = pow(u_next, 1.0/kappa);
    if (abs(xi) < 1e-10) x_rep[t] = -sigma * log1m(g_inv_p);
    else x_rep[t] = (sigma / xi) * (pow(1.0 - g_inv_p, -xi) - 1.0);
  }
}
