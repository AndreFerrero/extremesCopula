margin_egp <- list(
  name = "egp",
  name_param = c("mu", "kappa", "sigma", "xi"),
  
  G_dist = function(u, param) u^param["kappa"],
  G_inv  = function(u, param) u^(1 / param["kappa"]),
  g_dist = function(u, param) {
    param["kappa"] * u^(param["kappa"] - 1)
  },

  cdf = function(x, param) {
    x_std <- (x - param["mu"]) / param["sigma"]
    u <- evd::pgpd(x_std, shape = param["xi"])
    margin_egp$G_dist(u, param)
  },

  lpdf = function(x, param) {

    x_std <- (x - param["mu"]) / param["sigma"]

    u <- evd::pgpd(x_std, shape = param["xi"])

    log(margin_egp$g_dist(u, param)) +
      log(evd::dgpd(x_std, shape = param["xi"])) -
      log(param["sigma"])
  },

  quantile = function(p, param) {

    q_std <- evd::qgpd(
      margin_egp$G_inv(p, param),
      shape = param["xi"]
    )

    param["mu"] + param["sigma"] * q_std
  },

  sample = function(n, param) {
    U <- runif(n)
    margin_egp$quantile(U, param)
  }
)
