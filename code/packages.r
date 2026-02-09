# libs/packages.R
# ------------------------------
# Load all required packages
# Only loads packages, does NOT attempt installation
# This works both locally and on HPC/cluster
# ------------------------------

required_pkgs <- c(
  "copula",
  "stabledist",
  "mvtnorm",
  "MASS",
  "coda",
  "parallel",
  "here",
  "evd"
)

invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))
