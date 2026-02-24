functions {
  // --- GPD Functions (Internal Engine) ---
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

  // --- EGPD Functions (Now with mu integrated) ---
  
  // Log-CDF: F(x) = [G(x - mu)]^kappa
  real egpd_lcdf(real x, real mu, real kappa, real sigma, real xi) {
    return gpd_lcdf(x - mu | sigma, xi) * kappa;
  }

  // Log-PDF: f(x) = kappa * [G(x - mu)]^(kappa-1) * g(x - mu)
  real egpd_lpdf(real x, real mu, real kappa, real sigma, real xi) {
    real log_G = gpd_lcdf(x - mu | sigma, xi);
    real log_g = gpd_lpdf(x - mu | sigma, xi);
    return log(kappa) + (kappa - 1) * log_G + log_g;
  }

  // Quantile function (Inverse CDF)
  real egpd_quantile(real u, real mu, real kappa, real sigma, real xi) {
    real p_gpd = pow(u, 1.0/kappa);
    real excess;
    if (abs(xi) < 1e-5) {
      excess = -sigma * log1m(p_gpd);
    } else {
      excess = (sigma / xi) * (pow(1.0 - p_gpd, -xi) - 1.0);
    }
    return mu + excess;
  }

  // Random Number Generator
  real egpd_rng(real mu, real kappa, real sigma, real xi) {
    return egpd_quantile(uniform_rng(0, 1), mu, kappa, sigma, xi);
  }

  // --- Joe Copula Functions (Hofert et al. 2012 Notation) ---
  real joe_copula_lpdf(real u, real v, real theta) {
    if (u <= 1e-10 || u >= 1.0-1e-10 || v <= 1e-10 || v >= 1.0-1e-10) return -1e10;
    
    real alpha = 1.0 / theta;
    real u_term = pow(1.0 - u, theta);
    real v_term = pow(1.0 - v, theta);
    
    real om_hJ = u_term + v_term - (u_term * v_term);
    real hJ = (1.0 - u_term) * (1.0 - v_term);
    
    real log_poly = log1p((1.0 - alpha) * (hJ / om_hJ));
    
    return log(theta) + (theta - 1.0) * (log1m(u) + log1m(v)) 
           - (1.0 - alpha) * log(om_hJ) + log_poly;
  }
  
  real joe_hfunc(real u, real v, real theta) {
    real u_term = pow(1.0 - u, theta);
    real v_term = pow(1.0 - v, theta);
    real K = u_term + v_term - (u_term * v_term);
    
    return pow(K, (1.0/theta) - 1.0) * pow(1.0 - v, theta - 1.0) * (1.0 - pow(1.0 - u, theta));
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
    // First observation likelihood
    target += egpd_lpdf(x[1] | mu, kappa, sigma, xi);

    for (t in 2:T) {
      // Probability Integral Transform (PIT)
      real u = exp(egpd_lcdf(x[t] | mu, kappa, sigma, xi));
      real v = exp(egpd_lcdf(x[t-1] | mu, kappa, sigma, xi));

      // Marginal Density + Copula Transition Density
      target += egpd_lpdf(x[t] | mu, kappa, sigma, xi);
      target += joe_copula_lpdf(u | v, theta);
    }
  }
}

generated quantities {
  vector[T] x_rep;
  vector[T] log_lik;

  if (run_ppc == 1) {
    // Generate initial state
    x_rep[1] = egpd_rng(mu, kappa, sigma, xi);
    
    for (t in 2:T) {
      real v_prev = exp(egpd_lcdf(x_rep[t-1] | mu, kappa, sigma, xi));
      real w = uniform_rng(0, 1);
      real low = 1e-5; 
      real high = 1 - 1e-5;
      
      // Bisection solver for Joe Copula inversion (Rosenblatt)
      for (i in 1:I) {
        real mid = (low + high) / 2.0;
        if (joe_hfunc(mid, v_prev, theta) < w) low = mid; else high = mid;
      }
      real u_next = (low + high) / 2.0;
      
      // Map uniform back to physical scale
      x_rep[t] = egpd_quantile(u_next, mu, kappa, sigma, xi);
    }
  }

  // Pointwise log-likelihood for LOO-CV
  log_lik[1] = egpd_lpdf(x[1] | mu, kappa, sigma, xi);
  for (t in 2:T) {
    real u = exp(egpd_lcdf(x[t] | mu, kappa, sigma, xi));
    real v = exp(egpd_lcdf(x[t-1] | mu, kappa, sigma, xi));
    log_lik[t] = egpd_lpdf(x[t] | mu, kappa, sigma, xi) + joe_copula_lpdf(u | v, theta);
  }
}
