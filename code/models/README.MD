## Model Objects

Here we store model objects and useful builders

- `builders`: Builder functions that use model objects functionalities (Simulator with copula + margin choice, Log-posterior to be used in MCMC,...)
- `copulas`: Copula objects with related functionalities (simulate_u, log_density, log_prior...)
- `margins`: Margin objects with related functionalities (density, quantile function for simulator, cdf, log_prior)