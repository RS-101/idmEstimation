# Test script for PWC implementation with constant hazards
# We'll simulate data with known constant hazards and see if we can recover them

library(Rcpp)
library(RcppArmadillo)

# Set working directory to package root
setwd("/home/rasmus-emil/github/idmEstimation")

# Source the simulate_idm function
source("R/simulate_idm.R")

# Compile and load the C++ code
cat("Compiling C++ code...\n")
sourceCpp("src/piecewise_constant_hazard.cpp")

# Source the R wrapper
source("R/piecewise_constant_hazard.R")

cat("Compilation successful!\n\n")

# ============================================================================
# Simulate data with constant hazards
# ============================================================================

set.seed(123)

# True constant hazard rates
true_a12 <- 0.3
true_a13 <- 0.3
true_a23 <- 0.2

cat("=== True Hazard Rates ===\n")
cat("a12 (healthy -> ill):", true_a12, "\n")
cat("a13 (healthy -> death):", true_a13, "\n")
cat("a23 (ill -> death):", true_a23, "\n\n")

# Simulate data
cat("Simulating data with constant hazards...\n")
sim_result <- simulate_idm_constant_hazards(
  n = 10000,
  a12 = true_a12,
  a13 = true_a13,
  a23 = true_a23,
  target_pct_censoring = 0.2,
  average_number_of_visits = 8
)

# Extract the dataset
simulated_data <- sim_result$datasets$obs

cat("Data simulated successfully!\n")
cat("Sample size:", nrow(simulated_data), "\n\n")

# Check data structure
cat("=== Data Summary ===\n")
cat("Cases:\n")
cat("  Censored (healthy):", sum(simulated_data$status_dead == 0 & simulated_data$status_ill == 0), "\n")
cat("  Died (healthy):", sum(simulated_data$status_dead == 1 & simulated_data$status_ill == 0), "\n")
cat("  Ill:", sum(simulated_data$status_dead == 0 & simulated_data$status_ill == 1), "\n")
cat("  Died (ill):", sum(simulated_data$status_dead == 1 & simulated_data$status_ill == 1), "\n\n")

# Check for NA values
cat("NA values in V_ill:", sum(is.na(simulated_data$V_ill)), "\n")
cat("Range of V_0:", range(simulated_data$V_0), "\n")
cat("Range of V_healthy:", range(simulated_data$V_healthy), "\n")
cat("Range of V_ill (non-NA):", range(simulated_data$V_ill, na.rm = TRUE), "\n")
cat("Range of T_obs:", range(simulated_data$T_obs), "\n\n")

# ============================================================================
# Fit PWC model with few knots (should approximate constant hazards)
# ============================================================================

cat("=== Fitting PWC Model ===\n")
cat("Using 5 knots per transition (4 intervals)...\n\n")

tryCatch({
  fit_pwc <- fit_idm_pwc(
    data = simulated_data,
    n_knots = 5,
    verbose = TRUE
  )

  cat("\n=== Fitting Results ===\n")
  cat("Convergence:", ifelse(fit_pwc$fit$convergence == 0, "YES", "NO"), "\n")
  cat("Log-likelihood:", fit_pwc$fit$log_likelihood, "\n\n")

  # Extract estimated parameters
  lambda_12 <- fit_pwc$fit$lambda_hat$lambda_12
  lambda_13 <- fit_pwc$fit$lambda_hat$lambda_13
  lambda_23 <- fit_pwc$fit$lambda_hat$lambda_23

  cat("=== Estimated Hazard Rates ===\n")
  cat("lambda_12 (intervals):", round(lambda_12, 4), "\n")
  cat("lambda_13 (intervals):", round(lambda_13, 4), "\n")
  cat("lambda_23 (intervals):", round(lambda_23, 4), "\n\n")

  # Average estimated hazards (should be close to true values)
  avg_lambda_12 <- mean(lambda_12)
  avg_lambda_13 <- mean(lambda_13)
  avg_lambda_23 <- mean(lambda_23)

  cat("=== Average Estimated vs True ===\n")
  cat(sprintf("a12: Estimated = %.4f, True = %.4f, Error = %.4f\n",
              avg_lambda_12, true_a12, avg_lambda_12 - true_a12))
  cat(sprintf("a13: Estimated = %.4f, True = %.4f, Error = %.4f\n",
              avg_lambda_13, true_a13, avg_lambda_13 - true_a13))
  cat(sprintf("a23: Estimated = %.4f, True = %.4f, Error = %.4f\n\n",
              avg_lambda_23, true_a23, avg_lambda_23 - true_a23))

  # Plot comparison
  cat("Creating diagnostic plots...\n")
  times <- seq(0, max(simulated_data$T_obs), length.out = 100)

  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

  # Hazard functions
  plot(times, fit_pwc$hazards$a12(times), type = "s", col = "blue", lwd = 2,
       main = "Hazard 1->2", xlab = "Time", ylab = "Hazard",
       ylim = c(0, max(fit_pwc$hazards$a12(times)) * 1.1))
  abline(h = true_a12, col = "red", lty = 2, lwd = 2)
  legend("topright", c("Estimated", "True"), col = c("blue", "red"),
         lwd = 2, lty = c(1, 2))

  plot(times, fit_pwc$hazards$a13(times), type = "s", col = "blue", lwd = 2,
       main = "Hazard 1->3", xlab = "Time", ylab = "Hazard",
       ylim = c(0, max(fit_pwc$hazards$a13(times)) * 1.1))
  abline(h = true_a13, col = "red", lty = 2, lwd = 2)
  legend("topright", c("Estimated", "True"), col = c("blue", "red"),
         lwd = 2, lty = c(1, 2))

  plot(times, fit_pwc$hazards$a23(times), type = "s", col = "blue", lwd = 2,
       main = "Hazard 2->3", xlab = "Time", ylab = "Hazard",
       ylim = c(0, max(fit_pwc$hazards$a23(times)) * 1.1))
  abline(h = true_a23, col = "red", lty = 2, lwd = 2)
  legend("topright", c("Estimated", "True"), col = c("blue", "red"),
         lwd = 2, lty = c(1, 2))

  # Cumulative hazards
  true_A12 <- function(t) true_a12 * t
  true_A13 <- function(t) true_a13 * t
  true_A23 <- function(t) true_a23 * t

  plot(times, fit_pwc$hazards$A12(times), type = "l", col = "blue", lwd = 2,
       main = "Cumulative Hazard 1->2", xlab = "Time", ylab = "Cumulative Hazard")
  lines(times, true_A12(times), col = "red", lty = 2, lwd = 2)

  plot(times, fit_pwc$hazards$A13(times), type = "l", col = "blue", lwd = 2,
       main = "Cumulative Hazard 1->3", xlab = "Time", ylab = "Cumulative Hazard")
  lines(times, true_A13(times), col = "red", lty = 2, lwd = 2)

  plot(times, fit_pwc$hazards$A23(times), type = "l", col = "blue", lwd = 2,
       main = "Cumulative Hazard 2->3", xlab = "Time", ylab = "Cumulative Hazard")
  lines(times, true_A23(times), col = "red", lty = 2, lwd = 2)

  cat("\nTest completed successfully!\n")

}, error = function(e) {
  cat("\n!!! ERROR OCCURRED !!!\n")
  cat("Error message:", conditionMessage(e), "\n")
  cat("Error call:", deparse(conditionCall(e)), "\n")
  traceback()
})
