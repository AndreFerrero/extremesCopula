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

  // --- EGPD Functions
  // Log-CDF: F(x) = [G(x - mu)]^kappa
  real egpd_lcdf(real x, real mu, real kappa, real sigma, real xi) {
    return gpd_lcdf(x - mu | sigma, xi) * kappa;
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

  // Log-PDF: f(x) = kappa * [G(x - mu)]^(kappa-1) * g(x - mu)
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

  real egpd_rng(real mu, real kappa, real sigma, real xi) {
    return egpd_quantile(uniform_rng(0, 1), mu, kappa, sigma, xi);
  }

  // --- Gaussian Copula Density ---
  real gaussian_copula_lpdf(real u, real v, real rho) {
    real x = inv_Phi(u);
    real y = inv_Phi(v);
    real rho2 = square(rho);
    return -0.5 * log1m(rho2) - (rho2 * (square(x) + square(y)) - 2 * rho * x * y) / (2 * (1 - rho2));
  }

  real gaussian_copula_lpdf_vec(vector u, vector v, real rho) {
  
    int N = num_elements(u);
    
    vector[N] x = inv_Phi(u);
    vector[N] y = inv_Phi(v);
    
    real rho2 = square(rho);
    real log_det_term = -0.5 * log1m(rho2);
    
    vector[N] quad =
        rho2 * (square(x) + square(y))
        - 2.0 * rho .* (x .* y);
    
    return N * log_det_term
        - sum(quad) / (2.0 * (1.0 - rho2));
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
  real<lower=-0.99, upper=0.99> rho; 
}

model {
  mu ~ normal(5, 2.5); 
    kappa ~ lognormal(2, 1);
    sigma ~ exponential(0.1);
    xi ~ gamma(2, 10);
    rho ~ normal(0, 0.5);

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

    target += gaussian_copula_lpdf_vec(
                  u[2:T],
                  u[1:(T-1)],
                  rho
              );
  }
}

generated quantities {
  vector[T] x_rep;
  vector[T] log_lik;

  if (run_ppc == 1) {
    x_rep[1] = egpd_rng(mu, kappa, sigma, xi);
    for (t in 2:T) {
      real v_prev = exp(egpd_lcdf(x_rep[t-1] | mu, kappa, sigma, xi));
      real w = uniform_rng(0, 1);
      // Fast Analytical Inverse Rosenblatt for Gaussian
      real u_next = Phi(rho * inv_Phi(v_prev) + sqrt(1 - square(rho)) * inv_Phi(w));
      x_rep[t] = egpd_quantile(u_next, mu, kappa, sigma, xi);
    }
  }

  log_lik[1] = egpd_lpdf(x[1] | mu, kappa, sigma, xi);
  for (t in 2:T) {
    real u = exp(egpd_lcdf(x[t] | mu, kappa, sigma, xi));
    real v = exp(egpd_lcdf(x[t-1] | mu, kappa, sigma, xi));
    log_lik[t] = egpd_lpdf(x[t] | mu, kappa, sigma, xi) + gaussian_copula_lpdf(u | v, rho);
  }
}
