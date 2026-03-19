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
  "ggplot2",
  "posterior",
  "microbenchmark",
  "kableExtra"
)

invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

