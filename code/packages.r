# libs/packages.R
# ------------------------------
# Load all required packages
# Only loads packages, does NOT attempt installation
# This works both locally and on HPC/cluster
# ------------------------------

required_pkgs <- c(
  "copula",
  "MASS",
  "coda",
  "parallel",
  "here",
  "evd",
  "rstan",
  "bayesplot",
  "ggplot2"
)

invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

rstan_options(auto_write = TRUE)
