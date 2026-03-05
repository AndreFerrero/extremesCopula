
source("code/packages.r")
library(exdex)


spm(gumbel_sim$x, 1000)
spm(gaussian_sim$x, 1000)

?choose_b
b_vals <- c(2,3,4,5,6,8,9,10,12,15,16,18,20,24,30,32,36,40,45,48,54,60)
b_res <- choose_b(gumbel_sim$x, b_vals)
plot(b_res, ylim = c(0, 1))

# PPC on empirical extremal index
b_vals <- seq(10, 200, 10)
b_choice <- exdex::choose_b(gumbel_sim$x, b_vals)
plot(b_choice)

theta_gumbel <- spm(gumbel_sim$x, 50)
conf_theta_gumbel <- confint(theta_gumbel, interval_type = "lik")
plot(conf_theta_gumbel)

(gumbel_theta_ppc_20 <- ppc_stat(gumbel_sim$x, gumbel_fit$ppc,  stat = "stat_theta_b20"))
(gumbel_theta_ppc_50 <- ppc_stat(gumbel_sim$x, gumbel_fit$ppc,  stat = "stat_theta_b50"))
(gumbel_theta_ppc_100 <- ppc_stat(gumbel_sim$x, gumbel_fit$ppc,  stat = "stat_theta_b100"))

get_theoretical_theta_gumbel <- function(copula_param, n_sim = 5000) {
  
  # 1. Map Stan parameter to Beirlant/Smith parameter (alpha)
  alpha <- 1 / copula_param
  
  # 2. Handle the independence case 
  if (is.na(copula_param) || copula_param <= 1.0001) return(1.0)
  
  # 3. Generate U ~ Uniform(0,1)
  U <- runif(n_sim)
  
  max_prod <- rep(0, n_sim)
  current_prod <- rep(1, n_sim)
  
  # 4. Simulate the Tail Chain multipliers A
  # Logic from Beirlant (2004) Example 10.21
  for(i in 1:500) {
    V <- runif(n_sim)
    
    # Inverse CDF: a = (V^(1/(alpha - 1)) - 1)^(-alpha)
    A <- (V^(1/(alpha - 1)) - 1)^(-alpha)
    
    current_prod <- current_prod * A
    max_prod <- pmax(max_prod, current_prod)
    
    # Optimization: stop if the chains have decayed
    if (all(current_prod < 1e-9)) break
  }
  
  return(mean(max_prod <= U))
}

draws_vector <- as.numeric(as.matrix(gumbel_fit$copula_draws))

# Run the simulation for each draw
theta_ei_posterior <- sapply(draws_vector, function(a) {
  get_theoretical_theta_gumbel(a, n_sim = 10)
})
mean(theta_ei_posterior)

# Final Check: Look at the distribution
hist(theta_ei_posterior, 
     breaks = 30, 
     col = "steelblue", 
     border = "white",
     main = "Integrated Posterior of Extremal Index",
     xlab = expression(theta))
     
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

diag <- monitor_theoretical_theta(3, n_sim = 10000)

# 1. See when the 'downward drift' effectively killed the chains
print(paste("Simulation stopped at step:", diag$break_at))

# 2. Visualize the Sampled Tail Chains
# We plot the first 5 paths on a log scale to see the negative drift
matplot(t(diag$paths[1:20, ]), type = "l", lty = 1, log = "y",
        main = "Evolution of Tail Chain (Log Scale)",
        xlab = "Step (i)", ylab = "Product (A1 * ... * Ai)")
abline(h = 1e-9, col = "red", lty = 2) # The break threshold


monitor_individual_deaths <- function(copula_param, l = 100, n_sim = 5000) {
  alpha <- 1 / copula_param
  
  max_prod <- rep(0, n_sim)
  current_prod <- rep(1, n_sim)
  U <- runif(n_sim)
  
  # Track the step where each individual simulation "dies"
  death_steps <- rep(100, n_sim)
  still_alive <- rep(TRUE, n_sim)
  
  for(i in 1:l) {
    V <- runif(n_sim)
    A <- (V^(1/(alpha - 1)) - 1)^(-alpha)
    
    current_prod <- current_prod * A
    max_prod <- pmax(max_prod, current_prod)
    
    # Check which ones just died at this step
    just_died <- still_alive & (current_prod < 1e-9)
    death_steps[just_died] <- i
    still_alive[just_died] <- FALSE
    
    if (!any(still_alive)) break 
  }
  
  return(list(
    theta = mean(max_prod <= U),
    durations = death_steps
  ))
}

res <- monitor_individual_deaths(12, l = 5000)
hist(res$durations, main="Steps until cluster memory vanishes", xlab="Steps (Days)", breaks = 30)

# Posterior distribution of Extremal Index using tail chain theory
hist(gumbel_fit$draws[, "extremal_index"], breaks = 30)
median(gumbel_fit$draws[,"extremal_index"])
mean(gumbel_fit$draws[,"extremal_index"])
