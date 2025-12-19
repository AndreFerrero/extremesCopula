#' Marginal log-posterior for the Independence Case
#' 
#' @description 
#' This version ignores the copula density. It is used to estimate the 
#' margin parameters as if the data were i.i.d.
build_logposterior_indep <- function(margin, param_map, data,
                                           inverse_transform = NULL,
                                           log_jacobian = NULL,
                                           margin_prior = NULL) {
  
  function(param_init) {
    
    # 1. Apply inverse transform
    # param_init usually contains [mu, log_sigma] 
    param <- if (!is.null(inverse_transform)) inverse_transform(param_init) else param_init
    
    # Extract only margin parameters (theta is ignored)
    param_m <- param[param_map$margin]
    
    # --- Prior ---
    # Only evaluate the marginal prior
    logprior_m <- do.call(margin$log_prior, c(list(param_m), margin_prior))
    
    if (!is.finite(logprior_m)) return(-Inf)
    
    # --- Likelihood ---
    # In the independent case, Log-Likelihood = sum( log f(x_i) )
    # There is no copula contribution because log(c(u)) = log(1) = 0
    loglik <- sum(margin$log_density(data, param_m))
    
    # --- Jacobian adjustment ---
    logjac <- if (!is.null(log_jacobian)) log_jacobian(param_init) else 0
    
    return(loglik + logprior_m + logjac)
  }
}