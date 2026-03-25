functions {
  // =========================
  // GPD Helpers
  // =========================
  real gpd_lpdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-5) return -y - log(sigma);
    if (1 + xi * y <= 0) return negative_infinity();
    return -(1 + 1/xi) * log1p(xi * y) - log(sigma);
  }

  vector gpd_lcdf_vec(vector x, real sigma, real xi) {
    int N = num_elements(x);
    vector[N] out;
    for (n in 1:N) {
      real y = x[n] / sigma;
      if (abs(xi) < 1e-5) out[n] = log1m_exp(-y);
      else out[n] = log1m_exp(-(1.0/xi) * log1p(xi * y));
    }
    return out;
  }

  // =========================
  // EGPD Model (iv) 
  // F(x) = [1 - Q_delta( (1-GPD_cdf)^delta )]^(kappa/2)
  // =========================
  
  // Vectorized LCDF
  vector egpd4_lcdf_vec(vector x, real kappa, real delta, real sigma, real xi) {
    int N = num_elements(x);
    vector[N] log_G = gpd_lcdf_vec(x, sigma, xi);
    vector[N] out;
    real a = 1.0 / delta;
    
    for (n in 1:N) {
      // W = (1 - GPD_cdf)^delta
      real W = exp(delta * log1m_exp(log_G[n])); 
      // Q_delta is Beta CDF with params (1/delta, 2)
      real Q = beta_cdf(W| a, 2.0);
      out[n] = (kappa / 2.0) * log1m(Q);
    }
    return out;
  }

  // Vectorized LPDF
  real egpd4_lpdf_vec(vector x, real kappa, real delta, real sigma, real xi) {
    int N = num_elements(x);
    vector[N] log_G = gpd_lcdf_vec(x, sigma, xi);
    vector[N] out;
    real a = 1.0 / delta;

    for (n in 1:N) {
      real G = exp(log_G[n]);
      real W = pow(1.0 - G, delta);
      real Q = beta_cdf(W| a, 2.0);
      real log_q = beta_lpdf(W | a, 2.0);
      real log_gpd = gpd_lpdf(x[n] | sigma, xi);
      
      // Chain rule: d/dx F(x)
      // log_pdf = log(kappa/2) + (kappa/2 - 1)*log(1-Q) + log(q_beta) + (delta-1)*log(1-G) + log(delta) + log_gpd
      out[n] = log(kappa / 2.0) + (kappa / 2.0 - 1.0) * log1m(Q) + log_q + 
               (delta - 1.0) * log1m(G) + log(delta) + log_gpd;
    }
    return sum(out);
  }

  // Quantile function for PPC (Numerical)
  real egpd4_quantile(real u, real kappa, real delta, real sigma, real xi, int I) {
    // 1. Invert outer power: V = 1 - Q_delta(W) = u^(2/kappa) => Q_delta(W) = 1 - u^(2/kappa)
    real target_Q = 1.0 - pow(u, 2.0 / kappa);
    
    // 2. Invert Beta CDF (Q_delta) for W via bisection
    real low = 0.0; real high = 1.0;
    real a = 1.0 / delta;
    for (i in 1:I) {
      real mid = (low + high) / 2.0;
      if (beta_cdf(mid| a, 2.0) < target_Q) low = mid; else high = mid;
    }
    real W = (low + high) / 2.0;
    
    // 3. Invert W = (1-G)^delta => G = 1 - W^(1/delta)
    real G = 1.0 - pow(W, 1.0 / delta);
    
    // 4. Invert GPD
    return (sigma / xi) * (pow(1.0 - G, -xi) - 1.0);
  }

  // Gumbel Copula (Keep existing)
  real gumbel_copula_lpdf_vec(vector u, vector v, real theta) {
    int N = num_elements(u);
    real alpha = 1.0 / theta;
    vector[N] ln_u = -log(u); vector[N] ln_v = -log(v);
    vector[N] t = pow(ln_u, theta) + pow(ln_v, theta);
    vector[N] r = pow(t, alpha);
    return sum(-r + (theta - 1) * (log(ln_u) + log(ln_v)) + (alpha - 2) * log(t) + log(r + theta - 1) - (log(u) + log(v)));
  }

  real gumbel_hfunc(real u, real v, real theta) {
    real ln_u = -log(u); real ln_v = -log(v);
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
  real<lower=0> kappa;
  real<lower=0> delta;       
  real<lower=0> sigma;
  real xi;
  real<lower=0> thetam1;
}

transformed parameters {
  real theta = thetam1 + 1.0;
}

model {
  kappa ~ lognormal(2, 1);
  delta ~ lognormal(2, 1);
  sigma ~ exponential(0.1);
  xi ~ normal(0, 0.3);
  thetam1 ~ gamma(2, 1);

  if (prior_check == 0) {
    target += egpd4_lpdf_vec(x, kappa, delta, sigma, xi);
    vector[T] u = exp(egpd4_lcdf_vec(x, kappa, delta, sigma, xi));
    target += gumbel_copula_lpdf_vec(u[2:T], u[1:(T-1)], theta);
  }
}

generated quantities {
  vector[T] x_rep;
  if (run_ppc == 1) {
    x_rep[1] = egpd4_quantile(uniform_rng(0.01, 0.99), kappa, delta, sigma, xi, I);
    for (t in 2:T) {
      real u_prev = exp(egpd4_lcdf_vec(x_rep[(t-1):(t-1)], kappa, delta, sigma, xi)[1]);
      real w = uniform_rng(0.01, 0.99);
      real low_c = 1e-10; real high_c = 1.0 - 1e-10;
      for (i in 1:I) {
        real mid = (low_c + high_c) / 2.0;
        if (gumbel_hfunc(mid, u_prev, theta) < w) low_c = mid; else high_c = mid;
      }
      x_rep[t] = egpd4_quantile((low_c + high_c) / 2.0, kappa, delta, sigma, xi, I);
    }
  }
}
