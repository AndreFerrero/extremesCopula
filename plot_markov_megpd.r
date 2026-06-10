source("markov_megpd.r")

###########################################################
# 4. PLOTS
###########################################################

par(mfrow = c(2, 2))

# 1. Delta Function
r_seq <- seq(0, 20, length.out = 200)
plot(r_seq, delta_f_strong_upper(r_seq), type="l", main="Delta(r)", xlab="r", ylab="delta")

# 2. Trace
plot.ts(final_chain, main="Trace Plot", ylab="X_t")

# 3. Lag-1 Scatter
plot(log(x_tm1), log(x_t), xlab="log(X_{t-1})", ylab="log(X_t)", pch=16, col=rgb(0,0,1,0.5))

plot(u_tm1, u_t, 
     main = "Lag-1 Scatter Plot (U_{t-1} vs U_t)", 
     xlab = "U_{t-1}", ylab = "U_t", pch = 16, col = rgb(0, 0, 1, 0.5))

par(mfrow=c(1,2))
plot(u_tm1, u_t,
     xlim = c(0.9, 1), ylim = c(0.9, 1),
     main = "Lag-1 Scatter Plot (U_{t-1} vs U_t)", 
     xlab = "U_{t-1}", ylab = "U_t", pch = 16, col = rgb(0, 0, 1, 0.5))
plot(u_tm1, u_t,
     xlim = c(0, 0.1), ylim = c(0, 0.1),
     main = "Lag-1 Scatter Plot (U_{t-1} vs U_t)", 
     xlab = "U_{t-1}", ylab = "U_t", pch = 16, col = rgb(0, 0, 1, 0.5))

# 4. Marginal Verification
par(mfrow=c(1,1))
hist(final_chain, breaks=60, freq=FALSE, col="lightgrey", main="Marginal Fit")
x_seq <- seq(0, max(final_chain), length.out = 1000)
y_true <- sapply(x_seq, dmegpd_marginal_integrated, kappa=kappa_val, xi=xi_val, delta_func=delta_f)
lines(x_seq, y_true, col="red", lwd=3)
lines(density(final_chain, from = 0), col="blue", lwd=2)
###########################################################
# 5. TAIL DEPENDENCE BY THRESHOLD
###########################################################

# Function to calculate empirical chi
calc_chi_emp <- function(u1, u2, p_levels, type = "upper") {
  n <- length(u1)
  sapply(p_levels, function(p) {
    if (type == "upper") {
      # Prob(U1 > p & U2 > p) / (1 - p)
      sum(u1 > p & u2 > p) / (n * (1 - p))
    } else {
      # Prob(U1 < p & U2 < p) / p
      sum(u1 < p & u2 < p) / (n * p)
    }
  })
}

# Define threshold levels
p_upper <- seq(0.1, 0.99, length.out = 50)
p_lower <- seq(0.01, 0.99, length.out = 50)

chi_upper <- calc_chi_emp(u_tm1, u_t, p_upper, type = "upper")
chi_lower <- calc_chi_emp(u_tm1, u_t, p_lower, type = "lower")

# Plotting
par(mfrow = c(1, 2))

# Upper Tail Plot
plot(p_upper, chi_upper, type = "b", pch = 16, col = "blue",
     ylim = c(0, 1), xlab = "Threshold p", ylab = expression(chi(p)),
     main = "Upper Tail Dependence")
abline(h = 0, lty = 2)
grid()

# Lower Tail Plot
plot(p_lower, chi_lower, type = "b", pch = 16, col = "red",
     ylim = c(0, 1), xlab = "Threshold p", ylab = expression(chi[L](p)),
     main = "Lower Tail Dependence")
abline(h = 0, lty = 2)
grid()