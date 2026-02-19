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

  // --- Joe Copula Functions (Hofert et al. 2012 Notation) ---
  real joe_copula_lpdf(real u, real v, real theta) {
    if (u <= 1e-10 || u >= 1.0-1e-10 || v <= 1e-10 || v >= 1.0-1e-10) return -1e10;
    
    real alpha = 1.0 / theta;
    // Survival terms used in Hofert's hJ function
    real u_term = pow(1.0 - u, theta);
    real v_term = pow(1.0 - v, theta);
    
    // Joint survival component (hJ) and its complement (1 - hJ)
    // hJ = (1 - (1-u)^theta) * (1 - (1-v)^theta)
    // 1 - hJ = (1-u)^theta + (1-v)^theta - (1-u)^theta * (1-v)^theta
    real om_hJ = u_term + v_term - (u_term * v_term);
    real hJ = (1.0 - u_term) * (1.0 - v_term);
    
    // Bivariate Joe Polynomial: P_2,alpha(x) = 1 + (1 - alpha) * x
    // where x = hJ / (1 - hJ)
    real log_poly = log1p((1.0 - alpha) * (hJ / om_hJ));
    
    // Final Log Density (Corollary 1, Part 5)
    return log(theta) + (theta - 1.0) * (log1m(u) + log1m(v)) 
           - (1.0 - alpha) * log(om_hJ) + log_poly;
  }
  
  real joe_hfunc(real u, real v, real theta) {
    // Conditional CDF P(U <= u | V = v) derived from partial derivative
    real u_term = pow(1.0 - u, theta);
    real v_term = pow(1.0 - v, theta);
    real term_joint = u_term + v_term - (u_term * v_term);
    
    return pow(term_joint, (1.0/theta) - 1.0) * pow(1.0 - v, theta - 1.0);
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
  real<lower=0> xi; 
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
      target += joe_copula_lpdf(u| v, theta);
    }
  }
}

generated quantities {
  vector[T] x_rep;
  vector[T] log_lik;

  if (run_ppc == 1) {
    x_rep[1] = mu + egpd_rng(kappa, sigma, xi);
    for (t in 2:T) {
      real v_prev = exp(egpd_lcdf(x_rep[t-1] - mu | kappa, sigma, xi));
      real w = uniform_rng(0, 1);
      real low = 1e-5; real high = 1 - 1e-5;
      
      // Bisection solver for Joe Copula
      for (i in 1:I) {
        real mid = (low + high) / 2.0;
        if (joe_hfunc(mid, v_prev, theta) < w) low = mid; else high = mid;
      }
      real u_next = (low + high) / 2.0;
      
      real g_inv_p = pow(u_next, 1.0/kappa);
      real excess = (abs(xi) < 1e-5) ? -sigma * log1m(g_inv_p) : (sigma / xi) * (pow(1 - g_inv_p, -xi) - 1);
      x_rep[t] = mu + excess;
    }
  }

  // Pointwise log-likelihood for model selection
  log_lik[1] = egpd_lpdf(x[1] - mu | kappa, sigma, xi);
  for (t in 2:T) {
    real u = exp(egpd_lcdf(x[t] - mu | kappa, sigma, xi));
    real v = exp(egpd_lcdf(x[t-1] - mu | kappa, sigma, xi));
    log_lik[t] = egpd_lpdf(x[t] - mu | kappa, sigma, xi) + joe_copula_lpdf(u| v, theta);
  }
}