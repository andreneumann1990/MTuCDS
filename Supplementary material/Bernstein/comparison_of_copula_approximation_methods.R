library("doParallel")
library("doRNG")
library("copula")

source("functions.R")

# main ####

BC_Comparison <- function() {
  # Compares the approximation precision of the Bernstein copula with the
  # empirical copula in the settings of Omelka et al. (2009).
  # Note: The results are meant to be an addition to the large simulation study
  #       in Omelka et al. (2009).
  #
  # Variables:
  #   m: Dimension of the copula.
  #   K: Number of grid points.
  #   L: Number of test repetitions.
  #   n: Sample size.
  #
  # Returns:
  #   Null.
  t <- Sys.time()
  seed(603696) #603696
  #setting of Omelka et al. (2009) Section 3:
  m <- 2
  K <- 100

  #Kendall's tau = 1 - 1/theta*(1-D_1(theta)) = 0.25, D_1(theta) = ...
  model1 <- list(theta = 2.372, C = frankCopula(param = 2.372, dim = m))
  #Kendall's tau = theta/(theta+2) = 0.75
  model2 <- list(theta = 6, C = claytonCopula(param = 6, dim = m))

  models <- list(model1 = model1, model2 = model2)
  L <- 10000
  n <- 150

  #KS-distance (max) is used on a 101 x 101 grid
  grid <- list()
  for (j in 1:m) grid[[j]] <- 0:K
  grid <- as.matrix(expand.grid(grid))

  distance <- matrix(nrow = 2, ncol = nrow(grid))
  BC_Comparison <- list()
  registerDoParallel(cores = max(detectCores() - 1, 1))
  for (model_ in 1:2) {
    C <- models[[model_]]$C
    results_ <-
      foreach(l = 1:L, .export = lsf.str(.GlobalEnv),
              .packages = c("copula", "fBasics")) %dorng% {
        U <- rCopula(n, C)
        U_Bernstein <- rBernsteinCopula_U(100*n, U)
        cop <- pCopula(grid/K, C)
        distance[1, ] <- abs(C.n(grid/K, U) - cop)
        distance[2, ] <- abs(C.n(grid/K, U_Bernstein) - cop)
        rowMaxs(distance)
      }
    results <- t(as.data.frame(results_))
    colnames(results) <- c("Empricial", "Bernstein")
    rownames(results) <- 1:L
    BC_Comparison[[model_]] <- list(C = C, distance = results)
  }
  closeAllConnections()

  plot_results(BC_Comparison)
  save(BC_Comparison, file = paste(Sys.Date(), "-BC_Comparison.RData",
                                   sep = ""))
  print(Sys.time() - t)
  BC_Comparison
}


# functions ####

plot_results <- function(BC_Comparison) {
  # Creates and saves the boxplots of the comparison results.
  #
  # Args:
  #   BC_Comparison: List of simlulation results.
  #
  # Returns:
  #   Null.
  data_model1 <- list("Empirical" = BC_Comparison[[1]]$distance[, 1],
                      "Bernstein" = BC_Comparison[[1]]$distance[, 2])
  data_model2 <- list("Empirical" = BC_Comparison[[2]]$distance[, 1],
                      "Bernstein" = BC_Comparison[[2]]$distance[, 2])
  boxplot(data_model1)
  boxplot(data_model2)

  png(paste(Sys.Date(), "-BoxPlotModel1.png", sep = ""), width = 480,
      height = 480)
  par(mar = c(5.1, 4.4, 4.1, 2.1))
  boxplot(
    data_model1, main = "Kolmogorov-Smirnov distance", yaxt = 'n',
    pch = 19, cex.axis = 1.6, cex.lab = 1.6, cex.main = 1.6
  )
  axis(side = 2, at = (1:10)/100, cex.axis = 1.6)
  dev.off()
  png(paste(Sys.Date(), "-BoxPlotModel2.png", sep = ""), width = 480,
      height = 480)
  par(mar = c(5.1, 4.4, 4.1, 2.1))
  boxplot(
    data_model2, main = "Kolmogorov-Smirnov distance", yaxt = 'n',
    pch = 19, cex.axis = 1.6, cex.lab = 1.6, cex.main = 1.6
  )
  axis(side = 2, at = (1:10)/100, cex.axis = 1.6)
  dev.off()
}