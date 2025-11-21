functions {
  // Diagonal likelihood for a single block
  real diag_loglik_gumbel(real u, int n, real theta) {
    real inv_val = pow(-log(u), theta);  // ψ⁻¹(u)
    real psi_p   = (-1/theta) * pow(n * inv_val, 1/theta - 1) * exp(-pow(n * inv_val, 1/theta));  // ψ'(n ψ⁻¹(u))
    real inv_p   = -theta / u * pow(-log(u), theta - 1);  // derivative of ψ⁻¹(u)
    return log(n * psi_p * inv_p);
  }

  // Total log-likelihood over all blocks
  real total_diag_loglik(vector Y, int n, real theta) {
    int N = rows(Y);
    real ll = 0;
    for (i in 1:N) {
      ll += diag_loglik_gumbel(Y[i], n, theta);
    }
    return ll;
  }
}

data {
  int<lower=1> N;            // number of blocks
  int<lower=1> n;            // block size
  vector<lower=0, upper=1>[N] Y;  // maxima of each block
}

parameters {
  real<lower=0> theta_shift;  // theta = 1 + theta_shift > 1
}

transformed parameters {
  real<lower=1> theta = 1 + theta_shift;
}

model {
  // Prior: shifted gamma, favors smaller theta
  theta_shift ~ gamma(2, 0.25);

  // Likelihood
  target += total_diag_loglik(Y, n, theta);
}
