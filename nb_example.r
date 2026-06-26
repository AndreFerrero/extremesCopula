library("NeuralEstimators")
library("JuliaConnectoR")
library("ggplot2")
juliaEval('using NeuralEstimators, Flux')

sampler <- function(K) {
  mu    <- rnorm(K)
  sigma <- rgamma(K, 1)
  Theta <- matrix(c(mu, sigma), byrow = TRUE, ncol = K)
  return(Theta)
}

K <- 10000
theta_train <- sampler(K)
theta_val   <- sampler(K/10)
mu    <- theta_train[1, 1:4]
sigma <- theta_train[2, 1:4]
Z     <- Z_train[1:4]

simulate <- function(Theta, m) {
  apply(Theta, 2, function(theta) t(rnorm(m, theta[1], theta[2])), simplify = FALSE)
}

m <- 30
Z_train <- simulate(theta_train, m)
Z_val   <- simulate(theta_val, m)

estimator <- juliaEval('
  n = 1    # dimension of each data replicate (univariate)
  d = 2    # dimension of the parameter vector θ
  w = 128  # width of each hidden layer 
  
  # Final layer has output dimension d and enforces parameter constraints
  final_layer = Parallel(
      vcat,
      Dense(w, 1, identity),     # mean parameter: no constraints
      Dense(w, 1, softplus)      # standard-deviation parameter: strictly positive
  )

  # Construct inner and outer networks and combine into DeepSet
  psi = Chain(Dense(n, w, relu), Dense(w, w, relu))    
  phi = Chain(Dense(w, w, relu), final_layer)          
  network = DeepSet(psi, phi)

  # Wrap the neural network as a PointEstimator
  estimator = PointEstimator(network)
')

estimator <- train(
  estimator,
  theta_train = theta_train,
  theta_val   = theta_val,
  Z_train = Z_train,
  Z_val   = Z_val
  )

theta_test <- sampler(1000)
Z_test     <- simulate(theta_test, m)
assessment <- assess(estimator, theta_test, Z_test, estimator_names = "NBE")
plotestimates(assessment, parameter_labels = c("E1" = expression(mu), "E2" = expression(sigma)))

theta    <- as.matrix(c(0, 0.5))     # true parameters
Z        <- simulate(theta, m)       # pretend that this is observed data
estimate(estimator, Z)               # point estimate from the "observed data"

boot <- bootstrap(estimator, Z)

hist(boot[1,])
hist(boot[2,])



# Multivariate

# K: number of samples to draw from the prior
sampler <- function(K) {
  mu1    <- rnorm(K)
  mu2    <- rnorm(K)
  sigma  <- rgamma(K, 1)
  rho    <- runif(K, -1, 1)
  theta  <- matrix(c(mu1, mu2, sigma, rho), byrow = TRUE, ncol = K)
  return(theta)
}

# Marginal simulation from the statistical model
# theta: a matrix of parameters drawn from the prior
# m: number of conditionally independent replicates for each parameter vector
simulate <- function(Theta, m) {
  apply(Theta, 2, function(theta) {
    mu    <- c(theta[1], theta[2])
    sigma <- theta[3]
    rho   <- theta[4]
    Sigma <- sigma^2 * matrix(c(1, rho, rho, 1), 2, 2)
    Z <- MASS::mvrnorm(m, mu, Sigma)
    t(Z)
  }, simplify = FALSE)
}

# Initialise the estimator
estimator <- juliaEval('
  d = 2    # dimension of each replicate 
  w = 32   # number of neurons in each hidden layer
  
  # Layer to ensure valid estimates
  final_layer = Parallel(
      vcat,
      Dense(w, 2, identity), # mean parameters
      Dense(w, 1, softplus), # variance parameter
      Dense(w, 1, tanh)      # correlation parameter
    )

  psi = Chain(Dense(d, w, relu), Dense(w, w, relu), Dense(w, w, relu))
  phi = Chain(Dense(w, w, relu), Dense(w, w, relu), final_layer)
  deepset = DeepSet(psi, phi)
  estimator = PointEstimator(deepset)
')

theta <- sampler(3)
Z     <- simulate(theta, m = 250)
mu1   <- theta[1, 1:3]
mu2   <- theta[2, 1:3]
sigma <- theta[3, 1:3]
rho   <- theta[4, 1:3]

df <- Map(function(z, m1, m2, s, r) {
  data.frame(Z1 = z[1, ], Z2 = z[2, ], mu1 = m1, mu2 = m2, sigma = s, rho = r)
  }, Z, mu1, mu2, sigma, rho)
df <- do.call(rbind, df)

df$theta <- paste0(
  "mu1= ", round(df$mu1, 2), ", ",
  "mu2 = ", round(df$mu2, 2), ", ",
  "sigma = ", round(df$sigma, 2), ", ",
  "rho = ", round(df$rho, 2)
  )

df$theta <- as.factor(df$theta)

ggplot(df) + 
  geom_point(aes(x = Z1, y = Z2), alpha = 0.5) +
  facet_grid(~theta) +
  theme_bw()

K <- 5000
m <- 1000
theta_train <- sampler(K)
theta_val   <- sampler(K/10)
Z_train <- simulate(theta_train, m)
Z_val   <- simulate(theta_val, m)

# Train the estimator
estimator <- train(
  estimator,
  theta_train = theta_train,
  theta_val   = theta_val,
  Z_train = Z_train,
  Z_val   = Z_val, 
  epochs = 20
  )

theta_test <- sampler(1000)
Z_test     <- simulate(theta_test, m)
assessment <- assess(estimator, theta_test, Z_test, 
                     estimator_names = "NBE", 
                     parameter_names = c("mu1", "mu2", "sigma", "rho")) 
plotestimates(assessment)
