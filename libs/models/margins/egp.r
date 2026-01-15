margin_egp <- list(
  name = "egp",

  G_dist = function(u, param) u^param["kappa"],
  G_inv  = function(u, param) u^(1 / param["kappa"]),
  g_dist = function(u, param) {
    param["kappa"] * u^(param["kappa"] - 1)
  },

  cdf = function(x, param) {
    u <- evd::pgpd(x / param["sigma"], shape = param["xi"])
    margin_egp$G_dist(u, param)
  },

  lpdf = function(x, param) {
    u <- evd::pgpd(x / param["sigma"], shape = param["xi"])

    log(margin_egp$g_dist(u, param)) +
      log(evd::dgpd(x / param["sigma"], shape = param["xi"])) -
      log(param["sigma"])
  },

  quantile = function(p, param) {
    param["sigma"] * evd::qgpd(
      margin_egp$G_inv(p, param),
      shape = param["xi"]
    )
  },

  sample = function(n, param) {
    U <- runif(n)
    margin_egp$quantile(U, param)
  }
)
