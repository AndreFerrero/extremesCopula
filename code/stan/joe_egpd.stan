functions {
  // --- GPD Functions (Internal Engine) ---
  real gpd_lpdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-5) 
      return -y - log(sigma);
    if (1 + xi * y <= 0) 
      return -1e10;
    return -(1 + 1 / xi) * log1p(xi * y) - log(sigma);
  }
  
  real gpd_lcdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-5) 
      return log1m_exp(-y);
    if (1 + xi * y <= 0) 
      return -1e10;
    return log1m_exp(-(1 / xi) * log1p(xi * y));
  }
  
  // --- EGPD Functions
  // Log-CDF: F(x) = [G(x - mu)]^kappa
  real egpd_lcdf(real x, real mu, real kappa, real sigma, real xi) {
    return gpd_lcdf(x - mu | sigma, xi) * kappa;
  }
  
  vector egpd_lcdf_vec(vector x, real mu, real kappa, real sigma, real xi) {
    int N = num_elements(x);
    vector[N] out;
    
    for (n in 1 : N) 
      out[n] = kappa * gpd_lcdf(x[n] - mu | sigma, xi);
    
    return out;
  }
  
  // Log-PDF: f(x) = kappa * [G(x - mu)]^(kappa-1) * g(x - mu)
  real egpd_lpdf(real x, real mu, real kappa, real sigma, real xi) {
    real log_G = gpd_lcdf(x - mu | sigma, xi);
    real log_g = gpd_lpdf(x - mu | sigma, xi);
    return log(kappa) + (kappa - 1) * log_G + log_g;
  }
  
  real egpd_lpdf_vec(vector x, real mu, real kappa, real sigma, real xi) {
    int N = num_elements(x);
    vector[N] log_G;
    vector[N] log_g;
    
    for (n in 1 : N) {
      log_G[n] = gpd_lcdf(x[n] - mu | sigma, xi);
      log_g[n] = gpd_lpdf(x[n] - mu | sigma, xi);
    }
    
    return sum(log(kappa) + (kappa - 1) .* log_G + log_g);
  }
  
  // Quantile function (Inverse CDF)
  real egpd_quantile(real u, real mu, real kappa, real sigma, real xi) {
    real p_gpd = pow(u, 1.0 / kappa);
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
    if (u <= 1e-10 || u >= 1.0 - 1e-10 || v <= 1e-10 || v >= 1.0 - 1e-10) 
      return -1e10;
    
    real alpha = 1.0 / theta;
    real u_term = pow(1.0 - u, theta);
    real v_term = pow(1.0 - v, theta);
    
    real om_hJ = u_term + v_term - (u_term * v_term);
    real hJ = (1.0 - u_term) * (1.0 - v_term);
    
    real log_poly = log1p((1.0 - alpha) * (hJ / om_hJ));
    
    return log(theta) + (theta - 1.0) * (log1m(u) + log1m(v))
           - (1.0 - alpha) * log(om_hJ) + log_poly;
  }
  
  real joe_copula_lpdf_vec(vector u, vector v, real theta) {
    int N = num_elements(u);
    
    real alpha = 1.0 / theta;
    
    // Vectorised components
    vector[N] u_term = pow(1.0 - u, theta);
    vector[N] v_term = pow(1.0 - v, theta);
    
    vector[N] om_hJ = u_term + v_term - (u_term .* v_term);
    vector[N] hJ = (1.0 - u_term) .* (1.0 - v_term);
    
    vector[N] log_poly = log1p((1.0 - alpha) * (hJ ./ om_hJ));
    
    // Sum everything once
    return N * log(theta) + (theta - 1.0) * sum(log1m(u) + log1m(v))
           - (1.0 - alpha) * sum(log(om_hJ)) + sum(log_poly);
  }
  
  real joe_hfunc(real u, real v, real theta) {
    real u_term = pow(1.0 - u, theta);
    real v_term = pow(1.0 - v, theta);
    real K = u_term + v_term - (u_term * v_term);
    
    return pow(K, (1.0 / theta) - 1.0) * pow(1.0 - v, theta - 1.0)
           * (1.0 - pow(1.0 - u, theta));
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
  real<lower=0> kappa;
  real<lower=0> sigma;
  real xi;
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
  xi ~ normal(0, 0.3);
  thetam1 ~ gamma(2, 1);
  
  if (prior_check == 0) {
    // =========================
    // Vectorized marginal
    // =========================
    
    target += egpd_lpdf_vec(x, mu, kappa, sigma, xi);
    
    // =========================
    // Vectorized PIT
    // =========================
    
    vector[T] log_u = egpd_lcdf_vec(x, mu, kappa, sigma, xi);
    vector[T] u = exp(log_u);
    
    // =========================
    // Vectorized copula over pairs
    // =========================
    
    target += joe_copula_lpdf_vec(u[2 : T], u[1 : (T - 1)], theta);
  }
}
generated quantities {
  vector[T] x_rep;

  if (run_ppc == 1) {
    // Generate initial state
    x_rep[1] = egpd_rng(mu, kappa, sigma, xi);
    
    for (t in 2 : T) {
      real v_prev = exp(egpd_lcdf(x_rep[t - 1] | mu, kappa, sigma, xi));
      real w = uniform_rng(0, 1);
      real low = 1e-5;
      real high = 1 - 1e-5;
      
      // Bisection solver for Joe Copula inversion (Rosenblatt)
      for (i in 1 : I) {
        real mid = (low + high) / 2.0;
        if (joe_hfunc(mid, v_prev, theta) < w) 
          low = mid;
        else 
          high = mid;
      }
      real u_next = (low + high) / 2.0;
      
      // Map uniform back to physical scale
      x_rep[t] = egpd_quantile(u_next, mu, kappa, sigma, xi);
    }
  }
}
