functions {
  // --- GPD Functions ---
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

  // --- EGPD Functions ---
  real egpd_lcdf(real x, real kappa, real sigma, real xi) {
    return gpd_lcdf(x | sigma, xi) * kappa;
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

  // --- Gumbel Copula Functions ---
  real gumbel_copula_lpdf(real u, real v, real theta) {
    if (u <= 1e-10 || u >= 1.0-1e-10 || v <= 1e-10 || v >= 1.0-1e-10) return -1e10;
    real alpha = 1.0 / theta;
    real w_u = -log(u); real w_v = -log(v);
    real t = pow(w_u, theta) + pow(w_v, theta);
    real r = pow(t, alpha);
    return -r + (theta - 1) * (log(w_u) + log(w_v)) + (alpha - 2) * log(t) + log(r + theta - 1) - (log(u) + log(v));
  }
  
  real gumbel_hfunc(real u, real v, real theta) {
    real ln_u = -log(u); real ln_v = -log(v);
    real term = pow(ln_u, theta) + pow(ln_v, theta);
    return exp(-pow(term, 1.0/theta)) * pow(ln_v, theta - 1.0) * pow(term, 1.0/theta - 1.0) / v;
  }
}

data {
  int<lower=2> T;
  vector<lower=0>[T] x;
  int<lower=0, upper=1> prior_check; 
  int<lower=0, upper=1> run_ppc; 
  int<lower=16> I; 
}

parameters {
  real<lower=0, upper=min(x)> mu;
  real<lower=0.01> kappa;
  real<lower=0.01> sigma;
  real<lower=0, upper=0.5> xi; 
  real<lower=0> thetam1; 
}

transformed parameters {
  real theta = thetam1 + 1.0;
}

model {
  // Priors
  mu ~ normal(5, 2.5); 
  kappa ~ lognormal(2, 1);
  sigma ~ exponential(0.1);
  xi ~ gamma(2, 10);
  thetam1 ~ gamma(2, 1);

  if (prior_check == 0) {
    target += egpd_lpdf(x[1] - mu | kappa, sigma, xi);

    for (t in 2:T) {
      real u = exp(egpd_lcdf(x[t] - mu | kappa, sigma, xi));
      real v = exp(egpd_lcdf(x[t-1] - mu | kappa, sigma, xi));

      target += egpd_lpdf(x[t] - mu | kappa, sigma, xi);
      target += gumbel_copula_lpdf(u| v, theta);
    }
  }
}

generated quantities {
  vector[T] x_rep;
  vector[T] log_lik; // Added for LOO-CV model comparison

  if (run_ppc == 1) {
    x_rep[1] = mu + egpd_rng(kappa, sigma, xi);
    for (t in 2:T) {
      real v_prev = exp(egpd_lcdf(x_rep[t-1] - mu | kappa, sigma, xi));
      real w = uniform_rng(0, 1);
      real low = 1e-5; real high = 1 - 1e-5;
      for (i in 1:I) {
        real mid = (low + high) / 2.0;
        if (gumbel_hfunc(mid, v_prev, theta) < w) low = mid; else high = mid;
      }
      real u_next = (low + high) / 2.0;
      real g_inv_p = pow(u_next, 1.0/kappa);
      real excess = (abs(xi) < 1e-5) ? -sigma * log1m(g_inv_p) : (sigma / xi) * (pow(1 - g_inv_p, -xi) - 1);
      x_rep[t] = mu + excess;
    }
  }

  // Calculate pointwise log-likelihood for model selection (LOO-CV)
  log_lik[1] = egpd_lpdf(x[1] - mu | kappa, sigma, xi);
  for (t in 2:T) {
    real u = exp(egpd_lcdf(x[t] - mu | kappa, sigma, xi));
    real v = exp(egpd_lcdf(x[t-1] - mu | kappa, sigma, xi));
    log_lik[t] = egpd_lpdf(x[t] - mu | kappa, sigma, xi) + gumbel_copula_lpdf(u| v, theta);
  }
}
