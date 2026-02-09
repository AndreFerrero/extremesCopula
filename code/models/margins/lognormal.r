margin_lognormal <- list(

  name = "lognormal",

  quantile = function(u, param) {
    qlnorm(u, param["mu"], param["sigma"])
  },

  lpdf = function(x, param) {
    dlnorm(x, param["mu"], param["sigma"], log = TRUE)
  },

  cdf = function(x, param) {
    plnorm(x, param["mu"], param["sigma"])
  },

  log_prior = function(param, prior_mu = c(0, 10), prior_sigma = c(0, 1)) {
    if (param["sigma"] <= 0) return(-Inf)
    dnorm(param["mu"], prior_mu[1], prior_mu[2], log = TRUE) +
      dlnorm(param["sigma"], prior_sigma[1], prior_sigma[2], log = TRUE)
  }
)
