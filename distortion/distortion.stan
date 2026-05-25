functions {
  /**
   * Vectorized GEV Log-PDF
   */
  vector gev_vector_logpdf(vector x, real mu, real sigma, real xi) {
    int N = rows(x);
    vector[N] z = (x - mu) / sigma;
    vector[N] res;
    
    if (abs(xi) < 1e-8) {
      res = -log(sigma) - z - exp(-z);
    } else {
      vector[N] support = 1 + xi * z;
      // Stan will handle the reject/inf if any element of support <= 0
      res = -log(sigma) - (1/xi + 1) * log(support) - pow(support, -1/xi);
    }
    return res;
  }

  /**
   * Vectorized Distorted Joe Log-PDF
   * Returns a real sum to be added to target
   */
  real distorted_joe_all_lpdf(vector x, real mu, real sigma, real xi, real theta) {
    int K = rows(x);
    vector[K] z = (x - mu) / sigma;
    vector[K] V;     // -log(H)
    vector[K] log_h; // log(gev_pdf)

    if (abs(xi) < 1e-8) {
      V = exp(-z);
      log_h = -log(sigma) - z - V;
    } else {
      vector[K] support = 1 + xi * z;
      V = pow(support, -1/xi);
      log_h = -log(sigma) - (1/xi + 1) * log(support) - V;
    }

    vector[K] V_theta = pow(V, theta);
    
    // Vectorized Joe log(-psi_prime)
    // log(-psi_prime) = -log(theta) + (1/theta - 1) * log(1 - exp(-V^theta)) - V^theta
    vector[K] log_neg_psi_prime = -log(theta) + (1/theta - 1) * log1m_exp(-V_theta) - V_theta;
    
    // log(g) = log(theta) + (theta - 1)*log(V) + log_h + V + log_neg_psi_prime
    return sum(log(theta) + (theta - 1) * log(V) + log_h + V + log_neg_psi_prime);
  }
}

data {
  int<lower=0> K;             // Number of block maxima
  vector[K] M;                // Observed maxima
  int<lower=1, upper=2> family; // 1 = Gumbel, 2 = Joe
  
  // For Two-Stage: 1 to use fixed theta from DMLE, 0 to estimate
  int<lower=0, upper=1> use_fixed_theta;
  real<lower=1> theta_fixed;
}

parameters {
  real mu;
  real<lower=0> sigma;
  real xi; // Shape constrained to realistic EVT range
  real<lower=1> theta_param;  // Dependency parameter
}

transformed parameters {
  real<lower=1> theta;
  if (use_fixed_theta == 1) {
    theta = theta_fixed;
  } else {
    theta = theta_param;
  }
}

model {
  // Priors (Standardized)
  mu ~ normal(mean(M), sd(M) * 2);
  sigma ~ cauchy(0, 5);
  xi ~ normal(0.1, 0.2);
  
  if (!use_fixed_theta) {
    theta_param ~ gamma(2, 0.5);
  }

  // Vectorized Likelihood calls
  if (family == 1) {
    // family == 1: Gumbel distortion is just a GEV
    target += sum(gev_vector_logpdf(M, mu, sigma, xi));
  } else if (family == 2) {
    // family == 2: Joe distortion (non-GEV limit)
    target += distorted_joe_all_lpdf(M | mu, sigma, xi, theta);
  }
}

generated quantities {
  real p100 = 0.99;
  real return_level_100;
  
  if (family == 1) {
    return_level_100 = mu + (sigma / xi) * (pow(-log(p100), -xi) - 1);
  } else {
    real t = -log1m(pow(1 - p100, theta)); // psi_inv_joe using log1m for stability
    real p_adj = exp(-pow(t, 1/theta));     // Adjusted probability
    return_level_100 = mu + (sigma / xi) * (pow(-log(p_adj), -xi) - 1);
  }
}
