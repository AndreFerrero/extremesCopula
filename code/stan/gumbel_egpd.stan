functions {

  // =========================
  // GPD
  // =========================

  real gpd_lpdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-5)
      return -y - log(sigma);
    if (1 + xi * y <= 0)
      return negative_infinity();
    return -(1 + 1/xi) * log1p(xi * y) - log(sigma);
  }

  real gpd_lcdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (abs(xi) < 1e-5)
      return log1m_exp(-y);
    if (1 + xi * y <= 0)
      return negative_infinity();
    return log1m_exp(-(1/xi) * log1p(xi * y));
  }

  // =========================
  // Vectorized EGPD
  // =========================

  real egpd_lpdf(real x, real mu, real kappa, real sigma, real xi) {
    real log_G = gpd_lcdf(x - mu | sigma, xi);
    real log_g = gpd_lpdf(x - mu | sigma, xi);
    return log(kappa) + (kappa - 1) * log_G + log_g;
  }

  real egpd_lpdf_vec(vector x,
                real mu,
                real kappa,
                real sigma,
                real xi) {

    int N = num_elements(x);
    vector[N] log_G;
    vector[N] log_g;

    for (n in 1:N) {
      log_G[n] = gpd_lcdf(x[n] - mu | sigma, xi);
      log_g[n] = gpd_lpdf(x[n] - mu | sigma, xi);
    }

    return sum(log(kappa) + (kappa - 1) .* log_G + log_g);
  }

  real egpd_lcdf(real x,
                real mu,
                real kappa,
                real sigma,
                real xi) {

    return kappa * gpd_lcdf(x - mu | sigma, xi);
  }

  vector egpd_lcdf_vec(vector x,
                       real mu,
                       real kappa,
                       real sigma,
                       real xi) {

    int N = num_elements(x);
    vector[N] out;

    for (n in 1:N)
      out[n] = kappa * gpd_lcdf(x[n] - mu | sigma, xi);

    return out;
  }

  real egpd_quantile(real u,
                     real mu,
                     real kappa,
                     real sigma,
                     real xi) {

    real p_gpd = pow(u, 1.0 / kappa);
    real excess;

    if (abs(xi) < 1e-5)
      excess = -sigma * log1m(p_gpd);
    else
      excess = (sigma / xi) * (pow(1.0 - p_gpd, -xi) - 1.0);

    return mu + excess;
  }

  real egpd_rng(real mu,
                real kappa,
                real sigma,
                real xi) {

    return egpd_quantile(uniform_rng(0, 1), mu, kappa, sigma, xi);
  }

  // =========================
  // Vectorized Gumbel Copula
  // =========================
  real gumbel_copula_lpdf(real u, real v, real theta) {
    if (u <= 1e-10 || u >= 1.0-1e-10 || v <= 1e-10 || v >= 1.0-1e-10) return negative_infinity();
    real alpha = 1.0 / theta;
    real ln_u = -log(u);
    real ln_v = -log(v);
    real t = pow(ln_u, theta) + pow(ln_v, theta);
    real r = pow(t, alpha);
    return -r + (theta - 1) * (log(ln_u) + log(ln_v)) + (alpha - 2) * log(t) + log(r + theta - 1) - (log(u) + log(v));
  }

  real gumbel_copula_lpdf_vec(vector u,
                          vector v,
                          real theta) {

    int N = num_elements(u);

    real alpha = 1.0 / theta;

    vector[N] ln_u = -log(u);
    vector[N] ln_v = -log(v);

    vector[N] t = pow(ln_u, theta) + pow(ln_v, theta);
    vector[N] r = pow(t, alpha);

    return sum(
        -r
        + (theta - 1) .* (log(ln_u) + log(ln_v))
        + (alpha - 2) .* log(t)
        + log(r + theta - 1)
        - (log(u) + log(v))
    );
  }

  real gumbel_hfunc(real u,
                    real v,
                    real theta) {

    real ln_u = -log(u);
    real ln_v = -log(v);
    real t = pow(ln_u, theta) + pow(ln_v, theta);

    return exp(-pow(t, 1.0/theta))
           * pow(ln_v, theta - 1.0)
           * pow(t, 1.0/theta - 1.0)
           / v;
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

  // =========================
  // Priors
  // =========================

  mu ~ normal(5, 2.5);
  kappa ~ lognormal(2, 1);
  sigma ~ exponential(0.1);
  xi ~ gamma(2, 10);
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

    target += gumbel_copula_lpdf_vec(
                  u[2:T],
                  u[1:(T-1)],
                  theta
              );
  }
}
generated quantities {
  vector[T] x_rep;
  vector[T] log_lik; // Added for LOO-CV model comparison

  if (run_ppc == 1) {
    // Generate initial state
    x_rep[1] = egpd_rng(mu, kappa, sigma, xi);
    
    for (t in 2:T) {
      real v_prev = exp(egpd_lcdf(x_rep[t-1] | mu, kappa, sigma, xi));
      real w = uniform_rng(0, 1);
      real low = 1e-5; 
      real high = 1 - 1e-5;
      
      // Bisection solver for Copula inversion (Rosenblatt)
      for (i in 1:I) {
        real mid = (low + high) / 2.0;
        if (gumbel_hfunc(mid, v_prev, theta) < w) low = mid; else high = mid;
      }
      real u_next = (low + high) / 2.0;
      
      // Map uniform back to physical scale
      x_rep[t] = egpd_quantile(u_next, mu, kappa, sigma, xi);
    }
  }

  // Calculate pointwise log-likelihood for model selection (LOO-CV)
  log_lik[1] = egpd_lpdf(x[1] | mu, kappa, sigma, xi);
  for (t in 2:T) {
    real u = exp(egpd_lcdf(x[t] | mu, kappa, sigma, xi));
    real v = exp(egpd_lcdf(x[t-1] | mu, kappa, sigma, xi));
    log_lik[t] = egpd_lpdf(x[t] | mu, kappa, sigma, xi) + gumbel_copula_lpdf(u | v, theta);
  }
}