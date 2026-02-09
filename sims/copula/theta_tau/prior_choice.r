library(ggplot2)
library(patchwork)

# Settings
n_sim <- 100000
a <- 2
b <- 0.7
tau_func <- function(theta) 1 - 1/theta

# 1. Generate Samples
# Gamma choice: theta = Gamma(2, 0.7) + 1
stheta_gamma <- rgamma(n_sim, a, b) + 1
stau_gamma   <- tau_func(stheta_gamma)

# Uniform choice: theta = Uniform(0, 10) + 1  => theta in [1, 11]
stheta_unif  <- runif(n_sim, 0, 10) + 1
stau_unif    <- tau_func(stheta_unif)

# Combine into dataframes for ggplot
df_theta <- data.frame(
  val = c(stheta_gamma, stheta_unif),
  type = rep(c("Gamma(2, 0.7)", "Uniform(1, 11)"), each = n_sim)
)

df_tau <- data.frame(
  val = c(stau_gamma, stau_unif),
  type = rep(c("Gamma(2, 0.7)", "Uniform(1, 11)"), each = n_sim)
)

# 2. Plotting
p1 <- ggplot(df_theta, aes(x = val, fill = type)) +
  geom_density(alpha = 0.4) +
  labs(title = "Prior on Theta (Parameter Space)", x = expression(theta), y = "Density") +
  theme_minimal() + theme(legend.position = "none")

p2 <- ggplot(df_tau, aes(x = val, fill = type)) +
  geom_density(alpha = 0.4) +
  labs(title = "Implied Prior on Tau (Effect Space)", x = expression(tau), y = "Density") +
  theme_minimal() + labs(fill = "Prior Type")

gridExtra::grid.arrange(p1, p2)