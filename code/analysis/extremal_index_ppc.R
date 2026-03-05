library(exdex)

#' Extremal Index Stat Wrapper
#' @param x Vector of observations
#' @param b Block size (Tuning parameter)
#' @param type Version to extract (defaulting to the best choice)
calc_theta <- function(x, b, type = "BB2018b") {
  # Error handling: spm needs a minimum amount of data and exceedances
  res <- try({
    exdex::spm(x, b = b)
  }, silent = TRUE)
  
  if (inherits(res, "try-error")) {
    return(NA_real_)
  }
  
  # Return the sliding version of the chosen estimator
  return(as.numeric(res$theta_sl[type]))
}


# 2. The Function Factory
# This "locks in" the value of b for each specific function created
make_theta_stat <- function(b_val) {
  force(b_val) # Crucial: ensures the function remembers its specific b
  return(function(x) calc_theta(x, b = b_val))
}

# 3. Automatic Generation Loop
b_vals <- seq(10, 200, 10)

for(b in b_vals) {
  # Create a name like "stat_theta_b10", "stat_theta_b20", etc.
  func_name <- paste0("stat_theta_b", b)
  
  # Assign the generated function to that name in the global environment
  assign(func_name, make_theta_stat(b), envir = .GlobalEnv)
}