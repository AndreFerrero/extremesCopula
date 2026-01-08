build_maxima_distribution <- function(
  margin,
  copula,
  param_map
) {

  force(margin)
  force(copula)
  force(param_map)

  function(params, y, n) {

    # --- Extract parameters ---
    margin_param  <- params[param_map$margin]
    copula_param  <- params[param_map$copula]

    # --- Marginal CDF ---
    u <- margin$cdf(y, margin_param)
    u <- pmin(pmax(u, 1e-15), 1 - 1e-15)

    # --- Copula diagonal ---
    copula$diag(u, copula_param, n = n)
  }
}
