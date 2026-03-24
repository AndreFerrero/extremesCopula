functions {
  // =========================
  // GPD Core
  // =========================
  real gpd_lpdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-5) return -y - log(sigma);
    if (1 + xi * y <= 0) return negative_infinity();
    return -(1 + 1/xi) * log1p(xi * y) - log(sigma);
  }

  // Element-wise GPD LCDF for vectors
  vector gpd_lcdf_vec(vector x, real mu, real sigma, real xi) {
    int N = num_elements(x);
    vector[N] y = (x - mu) / sigma;
    vector[N] out;
    if (abs(xi) < 1e-5) {
      for (n in 1:N) out[n] = log1m_exp(-y[n]);
    } else {
      for (n in 1:N) {
        if (1 + xi * y[n] <= 0) out[n] = negative_infinity();
        else out[n] = log1m_exp(-(1.0/xi) * log1p(xi * y[n]));
      }
    }
    return out;
  }

  // =========================
  // EGPD Model (ii) Vectorized
  // =========================

  // Vectorized Mixture LPDF: log( p*kappa1*G^(kappa1-1) + (1-p)*kappa2*G^(kappa2-1) ) + log_g
  real egpd_mod2_lpdf_vec(vector x, real mu, real p, real kappa1, real kappa2, real sigma, real xi) {
    int N = num_elements(x);
    vector[N] log_G = gpd_lcdf_vec(x, mu, sigma, xi);
    vector[N] log_g;
    vector[N] log_dens_G;

    for (n in 1:N) {
      log_g[n] = gpd_lpdf(x[n] - mu | sigma, xi);
      // Mixture of power law densities
      log_dens_G[n] = log(p * kappa1 * exp((kappa1 - 1) * log_G[n]) + (1 - p) * kappa2 * exp((kappa2 - 1) * log_G[n]));
    }
    return sum(log_dens_G + log_g);
  }

  // Vectorized Mixture LCDF: log( p*G^kappa1 + (1-p)*G^kappa2 )
  vector egpd_mod2_lcdf_vec(vector x, real mu, real p, real kappa1, real kappa2, real sigma, real xi) {
    int N = num_elements(x);
    vector[N] log_G = gpd_lcdf_vec(x, mu, sigma, xi);
    vector[N] out;
    for (n in 1:N) {
      out[n] = log(p * exp(kappa1 * log_G[n]) + (1 - p) * exp(kappa2 * log_G[n]));
    }
    return out;
  }

  // Scalar version for PPC initial step
  real egpd_mod2_lcdf_scalar(real x, real mu, real p, real kappa1, real kappa2, real sigma, real xi) {
    real log_G = log1m_exp(-(1.0/xi) * log1p(xi * (x - mu) / sigma));
    return log(p * exp(kappa1 * log_G) + (1 - p) * exp(kappa2 * log_G));
  }

  // Numerical Inversion: Solve F(x) = u for x
  real egpd_mod2_quantile(real u, real mu, real p, real kappa1, real kappa2, real sigma, real xi, int I) {
    real low_v = 0.0;
    real high_v = 1.0;
    // Step 1: Solve p*v^kappa1 + (1-p)*v^kappa2 = u for v (GPD CDF value)
    for (i in 1:I) {
      real mid = (low_v + high_v) / 2.0;
      if (p * pow(mid, kappa1) + (1.0 - p) * pow(mid, kappa2) < u) low_v = mid;
      else high_v = mid;
    }
    real v = (low_v + high_v) / 2.0;
    // Step 2: GPD Inverse CDF
    real excess = (sigma / xi) * (pow(1.0 - v, -xi) - 1.0);
    return mu + excess;
  }

  // =========================
  // Gumbel Copula
  // =========================
  real gumbel_copula_lpdf_vec(vector u, vector v, real theta) {
    int N = num_elements(u);
    real alpha = 1.0 / theta;
    vector[N] ln_u = -log(u);
    vector[N] ln_v = -log(v);
    vector[N] t = pow(ln_u, theta) + pow(ln_v, theta);
    vector[N] r = pow(t, alpha);
    return sum(-r + (theta - 1) * (log(ln_u) + log(ln_v)) + (alpha - 2) * log(t) + log(r + theta - 1) - (log(u) + log(v)));
  }

  real gumbel_hfunc(real u, real v, real theta) {
    real ln_u = -log(u);
    real ln_v = -log(v);
    real t = pow(ln_u, theta) + pow(ln_v, theta);
    return exp(-pow(t, 1.0/theta)) * pow(ln_v, theta - 1.0) * pow(t, 1.0/theta - 1.0) / v;
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
  real<lower=0, upper=1> p;       
  real<lower=0> kappa1;           
  real<lower=0> kappa_diff;       
  real<lower=0> sigma;
  real xi;
  real<lower=0> thetam1;
}

transformed parameters {
  real kappa2 = kappa1 + kappa_diff; 
  real theta = thetam1 + 1.0;
}

model {
  // Priors
  mu ~ normal(5, 2.5);
  p ~ beta(2, 2);
  kappa1 ~ lognormal(2, 1);
  kappa_diff ~ lognormal(2, 1);
  sigma ~ exponential(0.1);
  xi ~ normal(0, 0.3);
  thetam1 ~ gamma(2, 1);

  if (prior_check == 0) {
    // 1. Marginal Likelihood (Vectorized)
    target += egpd_mod2_lpdf_vec(x, mu, p, kappa1, kappa2, sigma, xi);

    // 2. Copula Likelihood via Vectorized PIT
    vector[T] u = exp(egpd_mod2_lcdf_vec(x, mu, p, kappa1, kappa2, sigma, xi));
    
    target += gumbel_copula_lpdf_vec(u[2:T], u[1:(T-1)], theta);
  }
}

generated quantities {
  vector[T] x_rep;

  if (run_ppc == 1) {
    // Initial draw from marginal
    x_rep[1] = egpd_mod2_quantile(uniform_rng(0, 1), mu, p, kappa1, kappa2, sigma, xi, I);
    
    for (t in 2:T) {
      // 1. Current PIT value
      real u_prev = exp(egpd_mod2_lcdf_scalar(x_rep[t-1], mu, p, kappa1, kappa2, sigma, xi));
      real w = uniform_rng(0, 1);
      
      // 2. Invert Copula H-Function for u_next
      real low_c = 1e-10; 
      real high_c = 1.0 - 1e-10;
      for (i in 1:I) {
        real mid = (low_c + high_c) / 2.0;
        if (gumbel_hfunc(mid, u_prev, theta) < w) low_c = mid; else high_c = mid;
      }
      real u_next = (low_c + high_c) / 2.0;
      
      // 3. Invert Marginal Mixture for x_next
      x_rep[t] = egpd_mod2_quantile(u_next, mu, p, kappa1, kappa2, sigma, xi, I);
    }
  }
}
