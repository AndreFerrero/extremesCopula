functions {
  // 1. Standard GPD Log-PDF
  real gpd_lpdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (xi == 0) {
      return -y - log(sigma);
    } else {
      if (1 + xi * y <= 0) reject("x outside GPD support");
      return -(1 + 1/xi) * log1p(xi * y) - log(sigma);
    }
  }

  // 2. Standard GPD Log-CDF
  real gpd_lcdf(real x, real sigma, real xi) {
    real y = x / sigma;
    if (xi == 0) {
      return log1m(exp(-y));
    } else {
      return log1m(pow(1 + xi * y, -1/xi));
    }
  }

  // 3. EGPD Log-PDF: f(x) = kappa * G(x)^(kappa-1) * g(x)
  real egpd_lpdf(real x, real kappa, real sigma, real xi) {
    real log_G = gpd_lcdf(x | sigma, xi);
    real log_g = gpd_lpdf(x | sigma, xi);
    return log(kappa) + (kappa - 1) * log_G + log_g;
  }

  // 4. Bivariate Gumbel Copula Log-Density (d=2)
  // Derived from Hofert et al. (2012)
  real gumbel_copula_lpdf(real u, real v, real theta) {
    if (u <= 0 || u >= 1 || v <= 0 || v >= 1) return -1e10;
    
    real w_u = -log(u);
    real w_v = -log(v);
    real r = pow(pow(w_u, theta) + pow(w_v, theta), 1.0/theta);
    
    real term_top = log(r + theta - 1.0);
    real log_c = -r + (theta - 1.0) * (log(w_u) + log(w_v)) 
                 + (1.0 - 2.0 * theta) * log(r) 
                 + term_top - (log(u) + log(v));
    return log_c;
  }
}

data {
  int<lower=2> T;
  vector<lower=0>[T] x;
  
}

parameters {
  real<lower=0.01> kappa;
  real<lower=0.01> sigma;
  real<lower=-0.5, upper=0.5> xi; 
  real<lower=0> thetam1;        
}

transformed parameters {
    real theta = 1.0 + thetam1;
}

model {
  // Weakly informative: covers multiple orders of magnitude
  kappa ~ lognormal(0, 2); 
  sigma ~ lognormal(0, 2);
  
  // Very wide for EVT: allows for finite and infinite variance
  xi ~ normal(0, 0.5);
  
  thetam1 ~ gamma(2, 0.7); 

  // Likelihood
  // Initial observation margin
  target += egpd_lpdf(x[1] | kappa, sigma, xi);
  
  // Markov transitions
  for (t in 2:T) {
    real u = exp(gpd_lcdf(x[t]   | sigma, xi) * kappa);
    real v = exp(gpd_lcdf(x[t-1] | sigma, xi) * kappa);
    
    // Joint = Margins + Copula
    target += egpd_lpdf(x[t] | kappa, sigma, xi);
    target += gumbel_copula_lpdf(u | v, theta);
  }
}