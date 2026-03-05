monitor_theoretical_theta <- function(copula_param, n_sim = 5000, n_track = 20) {
  
  # 1. Setup parameters
  alpha <- 1 / copula_param
  if (is.na(copula_param) || copula_param <= 1.0001) return(list(theta = 1))
  
  # 2. Setup monitoring storage
  # Matrix to store the product evolution for n_track chains
  path_history <- matrix(NA, nrow = n_track, ncol = 100)
  
  # Generate U (the random levels we must stay below to "escape")
  U <- runif(n_sim)
  
  max_prod     <- rep(0, n_sim)
  current_prod <- rep(1, n_sim)
  break_step   <- 100 # Default if no early break
  
  # 3. The Simulation Loop
  for(i in 1:100) {
    V <- runif(n_sim)
    A <- (V^(1/(alpha - 1)) - 1)^(-alpha)
    
    # Update the chain
    current_prod <- current_prod * A
    
    # Update the running maximum (supremum)
    max_prod <- pmax(max_prod, current_prod)
    
    # Record history for the first n_track simulations
    path_history[, i] <- current_prod[1:n_track]
    
    # Monitor the break condition
    if (all(current_prod < 1e-9)) {
      break_step <- i
      break
    }
  }
  
  # 4. Return diagnostic list
  return(list(
    theta = mean(max_prod <= U),
    break_at = break_step,
    paths = path_history[, 1:break_step],
    U_subset = U[1:n_track]
  ))
}

diag <- monitor_theoretical_theta(2, n_sim = 10000)