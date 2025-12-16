#' Gaussian random-walk proposal
#'
#' @param Sigma Initial proposal covariance matrix (p x p)
#'
#' @return A proposal object with:
#'   - init(param): initializes internal proposal state
#'   - propose(param, state): proposes new parameter vector
#' 
source("libs/packages.R")

proposal_gaussian_rw <- function(Sigma) {

  list(

    #' Initialize proposal state
    #'
    #' @param param Initial parameter vector
    #' @return List storing proposal state (e.g. covariance)
    init_state = function(param) {
      list(
        Sigma0 = Sigma,
        Sigma  = Sigma,
        scale  = 1
      )
    },

    #' Propose a new parameter value
    #'
    #' @param param Current parameter vector
    #' @param state Current proposal state
    #'
    #' @return List with:
    #'   - param: proposed parameter
    #'   - log_q_ratio: log q(param | param') - log q(param' | param)
    propose = function(param, state) {

      eps <- MASS::mvrnorm(
        n = 1,
        mu = rep(0, length(param)),
        Sigma = state$Sigma
      )

      param_prop <- param + eps

      list(
        param = param_prop,
        log_q_ratio = 0  # symmetric proposal
      )
    }
  )
}
