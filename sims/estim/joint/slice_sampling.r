library(stabledist)

# Stable density at V
logp_V_theta <- function(V, theta) {
    dstable(V,
        alpha = 1 / theta,
        beta  = 1,
        gamma = (cos(pi / (2 * theta)))^theta,
        delta = 0,
        pm    = 1,
        log   = TRUE
    )
}

logp_U_Vtheta <- function(U, V, theta) {
    n <- length(U)

    n * log(V) + n * log(theta) - sum(log(U)) + (theta - 1) * sum(log(-log(U))) - V * sum((-log(U))^theta)
}

logp_mu <- function(mu) dnorm(mu, 0, 10, log = TRUE)

logp_sigma <- function(sigma) {
    if (sigma <= 0) {
        return(-Inf)
    }
    dlnorm(sigma, 0, 1, log = TRUE)
}

logp_theta <- function(theta) {
    if (theta <= 1) {
        return(-Inf)
    }
    dgamma(theta - 1, shape = 2, rate = 1, log = TRUE)
}

# -----------------------------------------------
# Generic slice sampler
# -----------------------------------------------
slice_sampler <- function(x0, logf, w = 0.2, m = 50, ...) {
    # Draw initial value defining the slice at random
    # y ~ Unif(0, f(x_0))
    # Using log f is safer
    # log(y) can be simulated using the exponential connection
    logy <- logf(x0, ...) - rexp(1)

    # Initial interval I = (L, R)
    L <- x0 - w * runif(1)
    R <- L + w

    # Step-out to increase I
    J <- floor(m * runif(1))
    K <- (m - 1) - J

    while (J > 0) {
        val <- logf(L, ...)
        if (!is.finite(val) || logy >= val) break
        L <- L - w
        J <- J - 1
    }

    while (K > 0) {
        val <- logf(R, ...)
        if (!is.finite(val) || logy >= val) break
        R <- R + w
        K <- K - 1
    }

    Lbar <- L
    Rbar <- R

    # Shrinkage until acceptance
    repeat {
        x1 <- Lbar + runif(1) * (Rbar - Lbar)

        val <- logf(x1, ...)
        if (!is.finite(val)) val <- -Inf
        if (logy <= val) {
            return(x1)
        }

        if (x1 < x0) Lbar <- x1 else Rbar <- x1
    }
}

# -----------------------------------------------
# Log posterior for V | U, theta
# -----------------------------------------------
logp_V_Utheta <- function(V, theta, U) {
    logp_V_theta(V, theta) + logp_U_Vtheta(U, V, theta)
}

# -----------------------------------------------
# Log posterior for theta | V, U
# -----------------------------------------------
logp_theta_UV <- function(theta, U, V) {
    if (theta <= 1) {
        return(-Inf)
    }

    logp_theta(theta) + logp_V_theta(V, theta) + logp_U_Vtheta(U, V, theta)
}


logp_mu_sigma_XVtheta <- function(mu, sigma, X, V, theta) {
    if (sigma <= 0) {
        return(-Inf)
    }

    U <- plnorm(X, meanlog = mu, sdlog = sigma)

    # log-likelihood of margins
    ll_marg <- sum(dlnorm(X, meanlog = mu, sdlog = sigma, log = TRUE))

    # log-likelihood of copula latent
    ll_cop <- logp_U_Vtheta(U, V, theta)

    logp_mu(mu) +
        logp_sigma(sigma) +
        ll_marg +
        ll_cop
}

# -----------------------------------------------
# Full latent-variable sampler
# -----------------------------------------------
gumbel_latent_slice <- function(U, n_iter = 5000, theta_init = 2, V_init = NULL) {
    n <- length(U)

    theta <- theta_init
    V <- if (is.null(V_init)) mean(-log(U)) else V_init

    theta_trace <- numeric(n_iter)
    V_trace <- numeric(n_iter)

    for (t in 1:n_iter) {
        # Sample V | theta, U
        V <- slice_sampler(V, logp_V_Utheta, w = 0.1, theta = theta, U = U)

        # Sample theta | V, U
        theta <- slice_sampler(theta, logp_theta_UV, w = 0.2, V = V, U = U)

        theta_trace[t] <- theta
        V_trace[t] <- V
    }

    list(theta = theta_trace, V = V_trace)
}

gumbel_latent_slice_margins <- function(X, n_iter = 5000,
                                        theta_init = 2, V_init = NULL,
                                        mu_init = mean(log(X)), sigma_init = sd(log(X))) {
    n <- length(X)

    theta <- theta_init
    V <- if (is.null(V_init)) mean(-log(plnorm(X, mu_init, sigma_init))) else V_init
    mu <- mu_init
    sigma <- sigma_init

    theta_trace <- numeric(n_iter)
    V_trace <- numeric(n_iter)
    mu_trace <- numeric(n_iter)
    sigma_trace <- numeric(n_iter)

    for (t in 1:n_iter) {
        # Sample V | theta, mu, sigma, X
        U <- plnorm(X, meanlog = mu, sdlog = sigma)
        V <- slice_sampler(V, logp_V_Utheta, w = 0.1, theta = theta, U = U)

        # Sample theta | V, mu, sigma, X
        theta <- slice_sampler(theta, logp_theta_UV, w = 0.1, V = V, U = U)

        # Sample mu | sigma, theta, V, X
        mu <- slice_sampler(mu, function(mu_) logp_mu_sigma_XVtheta(mu_, sigma, X, V, theta), w = 0.5)

        # Sample sigma | mu, theta, V, X
        sigma <- slice_sampler(sigma, function(sigma_) logp_mu_sigma_XVtheta(mu, sigma_, X, V, theta), w = 0.05)

        theta_trace[t] <- theta
        V_trace[t] <- V
        mu_trace[t] <- mu
        sigma_trace[t] <- sigma
    }

    list(theta = theta_trace, V = V_trace, mu = mu_trace, sigma = sigma_trace)
}

# Load packages
source("libs/packages.R")
# Load models
source("libs/models/copulas/gumbel.R")
source("libs/models/margins/lognormal.R")
source("libs/models/builders/simulator.R")

theta_true <- 2
mu_true <- 0
sigma_true <- 1
true_param <- c(mu = mu_true, sigma = sigma_true, theta = theta_true)

param_map <- list(margin = c("mu", "sigma"), copula = "theta")
simulator <- build_simulator(copula_gumbel, margin_lognormal, param_map)

n_sim <- 10000

U_sim <- copula_gumbel$simulate(theta_true, n_sim)

X_sim <- simulator(true_param, n_sim)


# run sampler for U
# res_U <- gumbel_latent_slice(U_sim, n_iter = 5000, theta_init = 1.5)

# run sampler for margin and copula

res_tot <- gumbel_latent_slice_margins(X_sim, n_iter = 100)

res <- res_tot

hist(res$theta, breaks = 40, main = "Posterior of theta", xlab = expression(theta))
abline(v = theta_true, col = "red", lwd = 2, lty = 2)
plot(1000:5000, res$theta[1000:5000], type = "l")
acf(res$theta[1000:5000])

plot(1000:5000, res$V[1000:5000], type = "l")

