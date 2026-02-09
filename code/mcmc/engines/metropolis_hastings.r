#' Single Metropolis–Hastings step
#'
#' @param param Current parameter vector
#' @param logpost Current log-target value
#' @param log_target Function computing log pi(param)
#' @param proposal Proposal object
#' @param prop_state Current proposal state
#'
#' @return List with updated param, logpost, accept indicator
#' 

mh_step <- function(
  param,
  logpost,
  log_target,
  proposal,
  prop_state
) {

  # Propose new parameter
  prop <- proposal$propose(param, prop_state)

  param_prop <- prop$param
  log_q_ratio <- prop$log_q_ratio

  # Evaluate target at proposal
  logp_prop <- log_target(param_prop)

  # Log acceptance probability
  log_alpha <- logp_prop - logpost + log_q_ratio

  if (log(runif(1)) < log_alpha) {
    list(
      param = param_prop,
      logpost = logp_prop,
      accept = TRUE
    )
  } else {
    list(
      param = param,
      logpost = logpost,
      accept = FALSE
    )
  }
}
