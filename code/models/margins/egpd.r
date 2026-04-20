margin_egpd <- list(
  name = "egpd",
  name_param = c("mu", "sigma", "xi", "kappa"),
  
  G_dist = function(u, param) u^param["kappa"],
  G_inv  = function(u, param) u^(1 / param["kappa"]),
  g_dist = function(u, param) {
    param["kappa"] * u^(param["kappa"] - 1)
  },

  cdf = function(x, param) {
    x_std <- (x - param["mu"]) / param["sigma"]
    u <- evd::pgpd(x_std, shape = param["xi"])
    margin_egpd$G_dist(u, param)
  },

  lpdf = function(x, param) {

    x_std <- (x - param["mu"]) / param["sigma"]

    u <- evd::pgpd(x_std, shape = param["xi"])

    log(margin_egpd$g_dist(u, param)) +
      log(evd::dgpd(x_std, shape = param["xi"])) -
      log(param["sigma"])
  },

  quantile = function(p, param) {

    q_std <- evd::qgpd(
      margin_egpd$G_inv(p, param),
      shape = param["xi"]
    )

    param["mu"] + param["sigma"] * q_std
  },

  sample = function(n, param) {
    U <- runif(n)
    margin_egpd$quantile(U, param)
  }
)

margin_egpd2 <- list(
  name = "egpd2",
  name_param = c("mu", "p", "kappa1", "kappa2", "sigma", "xi"),
  
  # G(v) = p*v^k1 + (1-p)*v^k2
  G_dist = function(u, param) {
    param["p"] * u^param["kappa1"] + (1 - param["p"]) * u^param["kappa2"]
  },
  
  # g(v) = p*k1*v^(k1-1) + (1-p)*k2*v^(k2-1)
  g_dist = function(u, param) {
    param["p"] * param["kappa1"] * u^(param["kappa1"] - 1) + 
    (1 - param["p"]) * param["kappa2"] * u^(param["kappa2"] - 1)
  },

  # Numerical inversion of G(v)
  G_inv = function(prob, param) {
    # Helper to find root for a single value
    invert_single <- function(p_val) {
      if (p_val <= 0) return(0)
      if (p_val >= 1) return(1)
      f <- function(v) margin_egpd2$G_dist(v, param) - p_val
      # G(v) is monotonic on [0,1], so uniroot is very stable here
      uniroot(f, interval = c(0, 1), tol = 1e-10)$root
    }
    # Vectorize the inversion
    vapply(prob, invert_single, numeric(1))
  },

  cdf = function(x, param) {
    # Handle values below mu
    x_std <- (x - param["mu"]) / param["sigma"]
    u <- evd::pgpd(x_std, shape = param["xi"])
    margin_egpd2$G_dist(u, param)
  },

  lpdf = function(x, param) {
    x_std <- (x - param["mu"]) / param["sigma"]
    u <- evd::pgpd(x_std, shape = param["xi"])
    
    # log(f(x)) = log(g(GPD_cdf)) + log(GPD_pdf) - log(sigma)
    log(margin_egpd2$g_dist(u, param)) +
      evd::dgpd(x_std, shape = param["xi"], log = TRUE) -
      log(param["sigma"])
  },

  quantile = function(p, param) {
    # 1. Invert the mixture to get the GPD-level probability
    v <- margin_egpd2$G_inv(p, param)
    
    # 2. Map through GPD quantile and scale/shift
    q_std <- evd::qgpd(v, shape = param["xi"])
    param["mu"] + param["sigma"] * q_std
  },

  sample = function(n, param) {
    U <- runif(n)
    margin_egpd2$quantile(U, param)
  }
)

margin_egpd4 <- list(
  name = "egpd4",
  name_param = c("mu", "kappa", "delta", "sigma", "xi"),
  
  # G(v) = [1 - Q_delta((1-v)^delta)]^(kappa/2)
  G_dist = function(u, param) {
    W <- (1 - u)^param["delta"]
    # Q_delta is Beta CDF with params (1/delta, 2)
    Q <- pbeta(W, shape1 = 1/param["delta"], shape2 = 2)
    (1 - Q)^(param["kappa"] / 2)
  },
  
  # G_inv is closed form thanks to qbeta
  G_inv = function(p, param) {
    target_Q <- 1 - p^(2 / param["kappa"])
    W <- qbeta(target_Q, shape1 = 1/param["delta"], shape2 = 2)
    1 - W^(1 / param["delta"])
  },

  cdf = function(x, param) {
    x_std <- (x - param["mu"]) / param["sigma"]
    u <- evd::pgpd(x_std, shape = param["xi"])
    margin_egpd4$G_dist(u, param)
  },

  lpdf = function(x, param) {
    x_std <- (x - param["mu"]) / param["sigma"]
    G_gpd <- evd::pgpd(x_std, shape = param["xi"])
    g_gpd <- evd::dgpd(x_std, shape = param["xi"])
    
    W <- (1 - G_gpd)^param["delta"]
    Q <- pbeta(W, shape1 = 1/param["delta"], shape2 = 2)
    q_beta <- dbeta(W, shape1 = 1/param["delta"], shape2 = 2)
    
    # log density calculation
    term1 <- log(param["kappa"] / 2) + (param["kappa"]/2 - 1) * log(1 - Q)
    term2 <- log(q_beta) + (param["delta"] - 1) * log(1 - G_gpd) + log(param["delta"])
    
    term1 + term2 + log(g_gpd) - log(param["sigma"])
  },

  quantile = function(p, param) {
    v <- margin_egpd4$G_inv(p, param)
    q_std <- evd::qgpd(v, shape = param["xi"])
    param["mu"] + param["sigma"] * q_std
  },

  sample = function(n, param) {
    U <- runif(n)
    margin_egpd4$quantile(U, param)
  }
)