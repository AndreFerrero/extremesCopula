?rstable
library(stabledist)
## Stable parameters for gumbel
# Margin F = Frechet(2)
n <- 1000

alpha <- 3
theta <- 1.1

gum_psi <- function(t, theta){
    exp(- t ^ (1/theta))
}

set.seed(46)

V <- rstable(
    n = 1,
    alpha = 1/theta,
    beta = 1,
    gamma = cospi(1/(2 * theta))^theta,
    delta = 0,
    pm = 1
)

E <- rexp(n, V) 

psi_E <- gum_psi(E, theta)

X <- qfrechet(psi_E)
Y <- qnorm(psi_E)

par(mfrow = c(1,3))
hist(psi_E)
plot(density(X))
plot(density(Y))
par(mfrow = c(1,1))

