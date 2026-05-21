# ==============================================================================
# SIMULATION FUNCTIONS
# ==============================================================================

library(copula)

#' Create copula object
#'
#' @param theta Copula parameter
#' @param n Dimension
#' @param family Copula family ("gumbel" or "joe")
#' @return Copula object
create_copula <- function(theta, n, family) {

  switch(family,
    gumbel = gumbelCopula(theta, dim = n),
    joe = joeCopula(theta, dim = n),
    stop("Unknown copula family: ", family)
  )
}

#' Simulate copula data and transform to margins
#'
#' @param k Number of replications (rows)
#' @param n Dimension (columns)
#' @param theta Copula parameter
#' @param family Copula family
#' @param margin_fn Marginal transformation function (default: t with df=4)
#' @param seed Random seed (optional)
#' @return List with U (uniform), X (marginal), M (block maxima)
simulate_data <- function(k, n, theta, family,
                          margin_fn = function(u) qt(u, df = 4),
                          seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  cop <- create_copula(theta, n, family)
  U <- rCopula(k, cop)
  X <- margin_fn(U)
  M <- apply(X, 1, max)

  list(U = U, X = X, M = M)
}

#' Compute pseudo-observations from data matrix
#'
#' @param X Data matrix (k x n)
#' @return List with Uhat (pseudo-observations), Yhat (row maxima)
compute_pseudo_obs <- function(X) {

  k <- nrow(X)
  n <- ncol(X)

  Uhat <- matrix(
    rank(as.vector(X)) / (k * n + 1),
    nrow = k,
    ncol = n
  )

  Yhat <- apply(Uhat, 1, max)

  list(Uhat = Uhat, Yhat = Yhat)
}
