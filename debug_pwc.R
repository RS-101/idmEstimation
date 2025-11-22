# Simple debug script to trace through the likelihood calculation

library(Rcpp)
sourceCpp("src/piecewise_constant_hazard.cpp")
source("R/piecewise_constant_hazard.R")
source("R/simulate_idm.R")

set.seed(42)
sim_result <- simulate_idm_constant_hazards(
  n = 20,
  a12 = 0.3,
  a13 = 0.1,
  a23 = 0.5,
  target_pct_censoring = 0.2,
  average_number_of_visits = 8
)

data <- sim_result$datasets$obs

cat("Data structure:\n")
print(str(data))
cat("\n")

# Setup model
model_config <- setup_cpp_model_pwc(
  V_0 = data$V_0,
  V_healthy = data$V_healthy,
  V_ill = data$V_ill,
  T_obs = data$T_obs,
  status_dead = data$status_dead,
  status_ill = data$status_ill,
  n_knots = 4  # Just 3 intervals
)

cat("Knots:\n")
cat("  knots_12:", model_config$knots_12, "\n")
cat("  knots_13:", model_config$knots_13, "\n")
cat("  knots_23:", model_config$knots_23, "\n\n")

# Try with true constant hazards
true_lambda_12 <- rep(0.3, 3)
true_lambda_13 <- rep(0.1, 3)
true_lambda_23 <- rep(0.5, 3)

cat("Testing with true constant hazards:\n")
ll_true <- calc_log_likelihood_pwc(
  model_config$model_pointer,
  true_lambda_12,
  true_lambda_13,
  true_lambda_23
)
cat("Log-likelihood with true values:", ll_true, "\n\n")

# Try with zeros
zero_lambda <- rep(0.01, 3)
ll_zero <- calc_log_likelihood_pwc(
  model_config$model_pointer,
  zero_lambda,
  zero_lambda,
  zero_lambda
)
cat("Log-likelihood with small values:", ll_zero, "\n\n")

# Now fit
cat("Fitting model...\n")
fit <- max_likelihood_pwc(model_config, lambda_0 = c(0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.5, 0.5, 0.5))

cat("\nFitted values:\n")
cat("lambda_12:", fit$lambda_hat$lambda_12, "\n")
cat("lambda_13:", fit$lambda_hat$lambda_13, "\n")
cat("lambda_23:", fit$lambda_hat$lambda_23, "\n")
cat("Log-likelihood:", fit$log_likelihood, "\n")
cat("Convergence:", fit$convergence, "\n")
