log_p <- function(x) {
  -x^4 + 2*x^2
}

slice_sample <- function(log_p, x0, n_samples, w = 1, m = 100) {
  samples <- numeric(n_samples)
  x <- x0
  
  for (i in 1:n_samples) {
    
    # Step 1: draw vertical level
    log_y <- log_p(x) + log(runif(1))
    
    # Step 2: initialize interval
    u <- runif(1, 0, w)
    L <- x - u
    R <- x + (w - u)
    
    # Step 3: stepping out
    j <- floor(runif(1, 0, m))
    k <- (m - 1) - j
    
    while (j > 0 && log_p(L) > log_y) {
      L <- L - w
      j <- j - 1
    }
    
    while (k > 0 && log_p(R) > log_y) {
      R <- R + w
      k <- k - 1
    }
    
    # Step 4: shrinkage
    repeat {
      x_new <- runif(1, L, R)
      
      if (log_p(x_new) >= log_y) {
        x <- x_new
        break
      } else {
        if (x_new < x) {
          L <- x_new
        } else {
          R <- x_new
        }
      }
    }
    
    samples[i] <- x
  }
  
  samples
}

mh_sample <- function(log_p, x0, n_samples, sigma = 0.5) {
  samples <- numeric(n_samples)
  x <- x0
  
  for (i in 1:n_samples) {
    x_prop <- rnorm(1, mean = x, sd = sigma)
    
    log_alpha <- log_p(x_prop) - log_p(x)
    
    if (log(runif(1)) < log_alpha) {
      x <- x_prop
    }
    
    samples[i] <- x
  }
  
  samples
}

set.seed(123)

n_samples <- 5000
samples_slice <- slice_sample(
  log_p = log_p,
  x0 = 0,
  n_samples = n_samples,
  w = 1
)

plot(samples_slice[1:500], type = "l",
     main = "Trace plot (first 500 iterations)",
     xlab = "Iteration", ylab = "x")

set.seed(123)

samples_mh <- mh_sample(
  log_p = log_p,
  x0 = 0,
  n_samples = n_samples,
  sigma = 2
)


# True density (unnormalized)
trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}

x_grid <- seq(-2.5, 2.5, length.out = 500)

p_true <- exp(log_p(x_grid))
p_true <- p_true / trapz(x_grid, p_true)


par(mfrow = c(1, 2))
hist(samples_slice, breaks = 60, probability = TRUE,
     col = "lightgray", border = "white",
     main = "Slice Sampling from a Bimodal Distribution",
     xlab = "x")

lines(x_grid, p_true, col = "red", lwd = 2)


# MH sampler
hist(samples_mh, breaks = 60, probability = TRUE,
     col = "lightgray", border = "white",
     main = "Metropolis–Hastings",
     xlab = "x")
lines(x_grid, p_true, col = "red", lwd = 2)


par(mfrow = c(2, 1))

plot(samples_slice[1:500], type = "l",
     main = "Slice sampler trace",
     xlab = "Iteration", ylab = "x")

plot(samples_mh[1:500], type = "l",
     main = "MH sampler trace",
     xlab = "Iteration", ylab = "x")
par(mfrow = c(1, 1))
