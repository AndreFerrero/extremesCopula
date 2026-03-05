# Install if needed
# install.packages("copula")

library(copula)

# Function that maps Joe theta -> equivalent Gumbel theta (via Kendall's tau)
joe_to_gumbel <- function(theta_joe_values, dim = 2) {
  
  results <- data.frame(
    theta_joe = numeric(),
    tau = numeric(),
    theta_gumbel = numeric(),
    lambdaU_joe = numeric(),
    lambdaU_gumbel = numeric()
  )
  
  for (theta_j in theta_joe_values) {
    
    # Define Joe copula
    joe_cop <- joeCopula(param = theta_j, dim = dim)
    
    # Compute Kendall's tau numerically
    tau_val <- tau(joe_cop)
    
    # Convert tau to Gumbel parameter
    theta_g <- 1 / (1 - tau_val)
    
    # Define Gumbel copula
    gumbel_cop <- gumbelCopula(param = theta_g, dim = dim)
    
    # Upper tail dependence
    lambdaU_j <- lambda(joe_cop)[["upper"]]
    lambdaU_g <- lambda(gumbel_cop)[["upper"]]
    
    results <- rbind(results, data.frame(
      theta_joe = theta_j,
      tau = tau_val,
      theta_gumbel = theta_g,
      lambdaU_joe = lambdaU_j,
      lambdaU_gumbel = lambdaU_g
    ))
  }
  
  return(results)
}

# Example grid of Joe parameters
theta_grid <- seq(1.1, 5, by = 0.3)

comparison_table <- joe_to_gumbel(theta_grid)

print(comparison_table)
