library("mvtnorm")
library("doParallel")
library("doRNG")
library("copula")

source("functions.R")

# main ####

BC_SimulationStudy <- function() {
  # Performs the simulation study.
  #
  # Variables:
  #   m: Number of null hypotheses.
  #   m0: Vector (arbitrary size) of numbers of true null hyptheses, each
  #       between 0 and m.
  #   n: Vector (arbitrary size) of sample sizes.
  #   C: List (arbitrary size) of copula functions.
  #   L: Number of repeated tests.
  #   M: Vector (arbitrary size) of Monte Carlo repetitions.
  #   alpha: Global significance level.
  #
  # Returns:
  #   Null.
  seed(26232) #26232

  m <- 20
  m0 <- c(0, 5, 10, 15, 20)
  n <- c(20, 60, 100)
  C <- list(frankCopula(param = 2, dim = m),
            frankCopula(param = 14, dim = m),
            gumbelCopula(param = 2, dim = m),
            gumbelCopula(param = 4, dim = m),
            claytonCopula(param = 1, dim = m),
            claytonCopula(param = 6, dim = m),
            joeCopula(param = 2, dim = m),
            joeCopula(param = 7, dim = m),
            tCopula(param = 0.4, dim = m),
            tCopula(param = 0.9, dim = m),
            normalCopula(param = 0.4, dim = m),
            normalCopula(param = 0.9, dim = m),
            indepCopula(dim = m))
  L <- 1000
  M <- c(200, 600, 1000)
  alpha <- 0.05

  BC_SimulationStudy <- compare_Methods(m, m0, n, C, L, M, alpha)
  BC_SimulationStudy <- evaluate_Data(BC_SimulationStudy)
  save(BC_SimulationStudy,
       file = paste(Sys.Date(), "-BC_SimulationStudy.RData", sep = ""))
  BC_SimulationStudy
}


# functions ####

compare_Methods <- function(m, m0, n, C, L, M, alpha) {
  # Compares the Bernstein method with Bonferroni and Sidak method.
  #
  # Args:
  #   m: Number of null hypotheses.
  #   m0: Vector (arbitrary size) of numbers of true null hyptheses, each
  #       between 0 and m.
  #   n: Vector (arbitrary size) of sample sizes.
  #   C: List (arbitrary size) of copula functions.
  #   L: Number of repeated tests.
  #   M: Vector (arbitrary size) of Monte Carlo repetitions.
  #   alpha: Global significance level.
  #
  # Returns:
  #   List (size depends) of the simulation results for each setting.
  BC_SimulationStudy <- list()
  registerDoParallel(cores = max(detectCores() - 1, 1))
  for (C_ in C) {
    print(C_)
    for (m0_ in m0) {
      theta <- c(numeric(m0_), numeric(m - m0_) + 0.4)
      for (n_ in n) {
        t <- Sys.time()
        T <- function(x, mu = 0)
          sqrt(n_) * abs(apply(x, 2, mean) - mu) / apply(x, 2, sd)
        for (M_ in M) {
          tests <-
            foreach(l = 1:L, .export = lsf.str(.GlobalEnv),
                    .packages = c("mvtnorm", "copula")) %dorng% {
              U <- rCopula(n_, C_)
              X <- matrix(nrow = n_, ncol = m)
              for (j in 1:m)
                X[, j] <- qnorm(U[, j], theta[j], sd = 1)
              alpha_loc <- calibrate_Bernstein_MC(X, T, M_, alpha)
              list(
                Bonferroni = test_p_values(X, T, alpha/m),
                Sidak  = test_p_values(X, T, 1 - (1 - alpha) ^ (1 / m)),
                Bernstein = test_p_values(X, T, alpha_loc)
              )
            }
          tests <- t(as.data.frame(tests))
          L_ <- nrow(tests)

          if (m0_ > 0 )
            FWER <- list(
              Bonferroni = mean(apply(tests[seq(1, L_, 3), 1:m0_], 1, sum) > 0),
              Sidak = mean(apply(tests[seq(2, L_, 3), 1:m0_], 1, sum) > 0),
              Bernstein = mean(apply(tests[seq(3, L_, 3), 1:m0_], 1, sum) > 0)
            )
          else
            FWER <- list(Bonferroni = 0, Sidak = 0, Bernstein = 0)
          if (m0_ < m)
            power <- list(
              Bonferroni = mean(apply(tests[seq(1, L_, 3), (m0_ + 1):m], 1,
                                      mean)),
              Sidak = mean(apply(tests[seq(2, L_, 3), (m0_ + 1):m], 1, mean)),
              Bernstein = mean(apply(tests[seq(3, L_, 3), (m0_ + 1):m], 1,
                                     mean))
            )
          else
            power <- list(Bonferroni = 1, Sidak = 1, Bernstein = 1)

          BC_SimulationStudy <- c(BC_SimulationStudy,
                                  list(list(C = C_, pi0 = m0_/m, M = M_, n = n_,
                                            FWER = FWER, power = power)))
          save(BC_SimulationStudy,
               file = paste(Sys.Date(), "-BC_SimulationStudy.RData", sep = ""))

          print(Sys.time() - t)
        }
      }
    }
  }
  closeAllConnections()
  BC_SimulationStudy
}

evaluate_Data <- function(BC_SimulationStudy) {
  # Changes the structure of the simulation results to a data frame.
  # Additionally, removes unnecessary entries of the copulas.
  #
  # Args:
  #   BC_SimulationStudy: List of the simulation results.
  #
  # Returns:
  #   Data frame of the simulation results.
  results <- BC_SimulationStudy
  for (i in 1:length(BC_SimulationStudy)) {
    param <- paste(BC_SimulationStudy[[i]]$C@parameters, collapse = "-")
    results[[i]]$C <- paste(class(BC_SimulationStudy[[i]]$C), "-param-",
                            param, sep = "")
  }
  nrow_ <- length(as.data.frame(results[[1]]))
  results <- t(as.data.frame(results))
  matrix(results, nrow = nrow_, dimnames = list(rownames(results)[1:nrow_],
                                                NULL))
}