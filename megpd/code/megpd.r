library(egpd)

rmegpd <- function(n, kappa, sigma, xi, delta) {

    r <- regpd(n, kappa = kappa, sigma = sigma, xi = xi)

    z <- rnorm(n, sd = delta(r))

    u1 <- exp(z) / (1 + exp(z))
    u2 <- 1 / (1 + exp(z))

    u <- cbind(u1, u2)

    x <- u * r

    colnames(x) <- c("x1", "x2")

    return(x)
}

delta_constant <- function(r) 0.8

delta_strong_upper <- function(r) {
  0.2 + 0.6 * exp(-r / 5)
}

delta_gamma <- function(r) 8 * dgamma(r, 2, 0.3) + 0.2
curve(delta_gamma, to = 50)

delta_upside <- function(r) (- 2.5) * dgamma(r, 9, 0.7) + 0.5
curve(delta_upside, to = 50)

x <- rmegpd(10000, 2, 1, 0.5, delta_gamma)

lx <- log(x)

kd <- MASS::kde2d(lx[,1], lx[,2], n = 200)

plot(lx)

# add contour lines
contour(kd, nlevels = 15, add = TRUE, drawlabels = FALSE, col = "blue", lwd = 2)

r <- apply(x, 1, sum)
hist(r)

1-pegpd(150, kappa = 2, sigma = 1, xi = 0.5)
round(1-mean(r <= 150), 5)

fitegpd(r)
