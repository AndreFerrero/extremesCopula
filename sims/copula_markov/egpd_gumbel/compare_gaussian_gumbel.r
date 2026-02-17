# --- 1. LIBRARIES & GLOBAL SETTINGS ---
library(rstan)
library(ggplot2)
library(bayesplot)
library(copula)
library(tidyr)
library(dplyr)
library(purrr)
library(evd)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

# --- 2. MATHEMATICAL ENGINE: EGPD & COPULA FUNCTIONS ---
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/margins/egp.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/margins/frechet.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/copulas/gumbel.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/copulas/gaussian.r")
source("C:/Users/Andrea Ferrero/extremesCopula/code/models/builders/markov/sim_copula_markov.R")

mod_sim <- simulate_copula_markov

# --- 3. GROUND TRUTH & DATA GENERATION ---

set.seed(46)

# 1. Generate a long series to see the density clearly
N_long <- 100000
rho_test <- 0.95
theta_test <- 8

sim_gauss <- simulate_copula_markov(N_long, copula_gaussian, rho_test, margin_egp, param_egp)
sim_gumbel <- simulate_copula_markov(N_long, copula_gumbel, theta_test, margin_egp, param_egp)

# 2. Plotting side-by-side in Uniform Space (U_t-1 vs U_t)
# This removes the marginal effect and shows ONLY the copula
par(mfrow=c(1,2))

plot(sim_gauss$U[-N_long], sim_gauss$U[-1], pch=16, col=rgb(0,0,1,0.1),
     main="Gaussian Copula (AI)", xlab="U[t-1]", ylab="U[t]")
abline(0,1, col="red")

plot(sim_gumbel$U[-N_long], sim_gumbel$U[-1], pch=16, col=rgb(1,0,0,0.1),
     main="Gumbel Copula (AD)", xlab="U[t-1]", ylab="U[t]")
abline(0,1, col="blue")

par(mfrow=c(1,1))

# 1. Filter only points in the top 5% corner (Correcting U/V order)
# x = previous state (Ut-1), y = current state (Ut)
thresh <- 0.99
corner_gauss  <- data.frame(prev = sim_gauss$U[-N_long], curr = sim_gauss$U[-1]) %>%
  filter(prev > thresh & curr > thresh)

corner_gumbel <- data.frame(prev = sim_gumbel$U[-N_long], curr = sim_gumbel$U[-1]) %>%
  filter(prev > thresh & curr > thresh)

# 2. Plotting
par(mfrow=c(1,2))

# --- GAUSSIAN PLOT ---
plot(corner_gauss$prev, corner_gauss$curr, xlim=c(thresh, 1), ylim=c(thresh, 1), 
     pch=16, col=rgb(0,0,1,0.5), main="Gaussian Corner (AI)",
     xlab=expression(U[t-1]), ylab=expression(U[t]))
abline(0, 1, col="grey", lty=2)

# --- GUMBEL PLOT ---
plot(corner_gumbel$prev, corner_gumbel$curr, xlim=c(thresh, 1), ylim=c(thresh, 1), 
     pch=16, col=rgb(1,0,0,0.5), main="Gumbel Corner (AD)",
     xlab=expression(U[t-1]), ylab=expression(U[t]))
abline(0, 1, col="grey", lty=2)

par(mfrow=c(1,1))
