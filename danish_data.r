# Install: install.packages("evir")
library(evir)

# Load the built-in dataset
data(danish)

# The data is a 'numeric' vector of claim sizes (in millions of Danish Krone)
insurance_vector <- as.numeric(danish)

# Get the timestamps (optional, for your Copula-Markov model)
# Note: These are 'event' times, not hourly regularly spaced data.
insurance_times <- attr(danish, "times") 

cat("Loaded", length(insurance_vector), "insurance claims.\n")
cat("Min claim value:", min(insurance_vector), "\n")

# Histogram on log-scale to see the heavy tail
hist(log(insurance_vector), breaks=30, col="salmon", main="Log-Insurance Claims")
hist(insurance_vector, breaks=50, col="salmon", main="Insurance Claims")

plot(insurance_times, insurance_vector, type="h", main="Insurance Claims Over Time", xlab="Time", ylab="Claim Size (Millions DKK)")
