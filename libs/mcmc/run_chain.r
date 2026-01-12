run_chain <- function(
  log_target,
  init,
  n_iter,
  proposal,
  engine_step = mh_step,
  adapt = adapt_none()
) {

  p <- length(init)

  param   <- init
  logpost <- log_target(param)

  samples <- matrix(NA, n_iter, p)
  accept  <- logical(n_iter)

  # Initialize proposal state
  prop_state <- proposal$init_state(param)

  # Adaptation bookkeeping
  adapting <- ifelse(adapt$type == "none", FALSE, TRUE)
  burnin   <- NA_integer_

  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

  for (i in seq_len(n_iter)) {

    # MH step
    step <- engine_step(
      param      = param,
      logpost    = logpost,
      log_target = log_target,
      proposal   = proposal,
      prop_state = prop_state
    )

    param        <- step$param
    logpost      <- step$logpost
    accept[i]    <- step$accept
    samples[i, ] <- param

    # Adapt only if still adapting
    if (adapting) {

      res <- adapt$update(prop_state, param, accept[i], i)
      prop_state <- res$state

      if (isTRUE(res$ready)) {
        burnin   <- i
        adapting <- FALSE
        adapt    <- adapt_none()
      }
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)
  colnames(samples) <- names(init)

  # If adaptation never stopped, treat whole chain as burn-in
  if (is.na(burnin)) {
    warning("Adaptation did not converge; entire chain treated as burn-in.")
    burnin <- n_iter
  }

  list(
    samples      = samples,
    burnin       = burnin,
    accept_rate  = mean(accept),
    conv_state   = prop_state
  )
}
