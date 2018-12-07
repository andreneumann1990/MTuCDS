library("mvtnorm")
library("doParallel")
library("doRNG")
library("copula")

# functions ####

seed <- function(seed = NULL) {
  # Sets a given or randomly generated seed.
  #
  # Args:
  #   seed: A given seed. If NULL then a randomly generated seed is chosen.
  #
  # Returns:
  #   NULL
  if (is.null(seed))
    seed <- as.numeric(paste(floor(runif(6, 0, 10)), collapse = ""))
  print(list(seed = seed))
  set.seed(seed)
}

calibrate_Bernstein_MC <-
  function(X, T, M, alpha) {
    # Calibrates the local significance level for the multiple test using the
    # Bernstein method.
    #
    # Args:
    #   X: Data matrix of size n x m.
    #   T: Test statistic function.
    #   M: Number of Monte Carlo repetitions.
    #   alpha: Global significance level.
    #
    # Returns:
    #   Vector (of length m) of local significance levels.
    n <- nrow(X); m <- ncol(X);

    U_star <- rBernsteinCopula(n*M, X)
    X_star <- matrix(nrow = n*M, ncol = m)
    for (j in 1:m)
      X_star[, j] <- qnorm(U_star[, j], sd = sd(X[, j]))
    T_star <- matrix(0, M, m)
    for (l in 1:M)
      T_star[l,] <- T(X_star[((l - 1) * n + 1):(l * n),])

    V_star <- 2 * pt(T_star, df = n - 1) - 1
    1 - emvQuantile(V_star, 1 - alpha)
  }

emvQuantile <- function(X, u) {
  # Calculates the (implicitly weighted) empirical multivariate quantile at
  # point u.
  #
  # Args:
  #   X: Data matrix of size n x m.
  #   u: Number between 0 and 1.
  #
  # Returns:
  #   Vector of length m.
  n <- nrow(X)
  r <- ceiling(n * u)

  rankMatrix = apply(X, 2, rank, ties.method = "first")
  rowRankVector = rank(apply(rankMatrix, 1, max), ties.method = "first")
  apply(X[rowRankVector <= r, ], 2, max)
}

test_p_values <- function(X, T, alpha_loc) {
  # Calculates and tests p-values for two-sided marginal t-tests.
  #
  # Args:
  #   X: Data matrix of size n x m.
  #   T: Test statistic function.
  #   alpha_loc: Vector (of length m) of local significance levels.
  #
  # Returns:
  #   Vector (logical, length m) of test results for the multiple test.
  n = nrow(X)
  2 * (1 - pt(T(X), n - 1)) < alpha_loc
}


# Cottin and Pfeifer (2014) ####

rBernsteinCopula <- function(N, X) {
  # Samples from the Bernstein copula using the method described in
  # Cottin and Pfeifer (2014).
  #
  # Args:
  #   N: Sample size.
  #   X: Data matrix of size n x m.
  #
  # Returns:
  #   Sample matrix (of size N x m) of the Bernstein copula.
  n <- nrow(X); m <- ncol(X);
  K <- n
  M <- apply(X, 2, rank, ties.method = "first") - 1
  U <- matrix(nrow = N, ncol = m)
  for (i in 1:N) {
    k <- M[sample(1:n, 1), ]
    for (j in 1:m)
      U[i, j] = rbeta(1, k[j] + 1, K - k[j])
  }
  U
}

rBernsteinCopula_U <- function(N, U) {
  # Samples from the Bernstein copula under known and continous marginal
  # distribution functions.
  #
  # Args:
  #   N: Sample size.
  #   U: Sample matrix (of size n x m) of the true copula.
  #
  # Returns:
  #   Sample matrix (of size N x m) of the Bernstein copula.
  n <- nrow(U); m <- ncol(U);
  K <- n
  M <- apply(U, 2, function(u) floor(u * K))
  U <- matrix(nrow = N, ncol = m)
  for (i in 1:N) {
    k <- M[sample(1:n, 1), ]
    for (j in 1:m)
      U[i, j] = rbeta(1, k[j] + 1, K - k[j])
  }
  U
}