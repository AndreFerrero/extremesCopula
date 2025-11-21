library(rstan)

# ----------------------------
# Simulation parameters
# ----------------------------
set.seed(123)
nk <- 200                  # block size (dimension of copula)
theta_true <- 6           # true Gumbel copula parameter
K <- 30                   # number of blocks

plots_folder <- here("sims", "estim", "block", "exp_stan_diag", "results", "plots")
exp_folder <- here("sims", "estim", "block", "exp_stan_diag")

# ----------------------------
# Functions
# ----------------------------
simulate_blocks <- function(K, cop) {
  lapply(1:K, function(b) as.numeric(rCopula(1, cop)))
}

psi_inv_gumbel <- function(u, theta) {
  (-log(u))^theta
}

psi_prime_gumbel <- function(t, theta) {
  (-1/theta) * t^(1/theta - 1) * exp(-t^(1/theta))
}

psi_inv_prime_gumbel <- function(u, theta) {
  -theta / u * (-log(u))^(theta - 1)
}

diag_prime_gumbel <- function(u, n, theta) {
  inv_val <- psi_inv_gumbel(u, theta)
  psi_p <- psi_prime_gumbel(n * inv_val, theta)
  inv_p <- psi_inv_prime_gumbel(u, theta)
  log(n * psi_p * inv_p)
}

# ----------------------------
# Simulate data
# ----------------------------
cop <- gumbelCopula(theta_true, dim = nk)
blocks <- simulate_blocks(K, cop)
Ymax <- sapply(blocks, max)

# ----------------------------
# Prepare data for Stan
# ----------------------------
stan_data <- list(
  N = K,
  Y = Ymax,
  n = nk
)

# ----------------------------
# Fit Bayesian model in Stan
# ----------------------------
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_fit <- stan(
  file = here(exp_folder, "diag_stan.stan"),
  data = stan_data,
  iter = 16000,
  thin = 4,
  chains = 4,
  control = list(adapt_delta = 0.99)
)

library(bayesplot)

# Extract posterior samples
posterior_samples <- as.data.frame(rstan::extract(stan_fit, "theta"))
colnames(posterior_samples) <- "theta"

# Compute 95% credible interval
ci <- quantile(posterior_samples$theta, probs = c(0.025, 0.975))

# Posterior density plot
# ggplot(posterior_samples, aes(x = theta)) +
#   geom_density(fill = "skyblue", alpha = 0.5, color = "blue", size = 1) +
#   geom_vline(xintercept = theta_true, color = "red", linetype = "dashed", size = 1) +
#   geom_vline(xintercept = ci[1], color = "darkblue", linetype = "dotted", size = 1) +
#   geom_vline(xintercept = ci[2], color = "darkblue", linetype = "dotted", size = 1) +
#   annotate("text", x = theta_true, y = 0.5 * max(density(posterior_samples$theta)$y),
#            label = "True θ", color = "red", hjust = -0.1) +
#   annotate("text", x = ci[1], y = 0.4 * max(density(posterior_samples$theta)$y),
#            label = "2.5%", color = "darkblue", hjust = 1.1) +
#   annotate("text", x = ci[2], y = 0.4 * max(density(posterior_samples$theta)$y),
#            label = "97.5%", color = "darkblue", hjust = -0.1) +
#   labs(title = "Posterior density of Gumbel copula parameter θ",
#        x = expression(theta),
#        y = "Density") +
#   theme_minimal(base_size = 14)

# --- Prior density (shifted Gamma) ---
shape <- 2
rate  <- 0.25
theta_seq <- seq(1, 30, length.out = 1000)
theta_shift_seq <- theta_seq - 1
prior_density <- dgamma(theta_shift_seq, shape = shape, rate = rate)
prior_density[theta_seq < 1] <- 0
prior_df <- data.frame(theta = theta_seq, density = prior_density)

# --- Combined plot ---
ggplot() +
  # Posterior
  geom_density(data = posterior_samples, aes(x = theta, y = ..density..),
               fill = "skyblue", alpha = 0.5, color = "blue") +
  # Prior
  geom_line(data = prior_df, aes(x = theta, y = density),
            color = "darkgreen") +
  geom_area(data = prior_df, aes(x = theta, y = density),
            fill = "lightgreen", alpha = 0.3) +
  # True theta and credible interval
  geom_vline(xintercept = theta_true, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = ci[1], color = "darkblue", linetype = "dotted", size = 1) +
  geom_vline(xintercept = ci[2], color = "darkblue", linetype = "dotted", size = 1) +
  annotate("text", x = theta_true, y = 0.5 * max(density(posterior_samples$theta)$y),
           label = "True θ", color = "red", hjust = -0.1) +
  annotate("text", x = ci[1], y = 0.4 * max(density(posterior_samples$theta)$y),
           label = "2.5%", color = "darkblue", hjust = 1.1) +
  annotate("text", x = ci[2], y = 0.4 * max(density(posterior_samples$theta)$y),
           label = "97.5%", color = "darkblue", hjust = -0.1) +
  labs(title = "Prior and Posterior for Gumbel copula parameter θ",
       x = expression(theta),
       y = "Density") +
  theme_minimal(base_size = 14)

ggsave(here(plots_folder, "density.png"))

# Extract posterior samples
posterior_array <- as.array(stan_fit)  # array of shape (iterations, chains, parameters)

# Traceplot for theta
trace <- mcmc_trace(posterior_array, pars = "theta") +
  ggtitle("Traceplot") +
  theme_minimal(base_size = 14)

# Autocorrelation function for theta
acf <- mcmc_acf(posterior_array, pars = "theta") +
  ggtitle("ACF with thinning = 4") +
  theme_minimal(base_size = 14)

bayes_check <- gridExtra::grid.arrange(trace, acf, ncol = 2)

ggsave(plot = bayes_check,filename =  here(plots_folder, "bayes_check.png"))

# Summary of posterior
fit_summary <- summary(stan_fit)$summary
theta_summary <- fit_summary["theta", ]

round(theta_summary, 2)
