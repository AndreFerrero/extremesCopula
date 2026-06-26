# file: r_benchmark.R

f <- function(x) {
  sin(x^2)
}

N <- 100000

time <- system.time({
  results <- numeric(N)

  for (i in 1:N) {
    results[i] <- integrate(f, lower = 0, upper = 10)$value
  }
})

cat("R mean result:", mean(results), "\n")
print(time)
