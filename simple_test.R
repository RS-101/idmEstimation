# Ultra-simple test with manual calculation

library(Rcpp)
sourceCpp("src/piecewise_constant_hazard.cpp")
source("R/piecewise_constant_hazard.R")

# Create ultra-simple data: 1 person, censored healthy
simple_data <- data.frame(
  V_0 = 0,
  V_healthy = 0,
  V_ill = NA,
  T_obs = 1,
  status_dead = 0,
  status_ill = 0
)

cat("Simple data (1 person, censored healthy at T=1):\n")
print(simple_data)
cat("\n")

# Setup with 2 knots (1 interval) - should have constant hazard
knots <- c(0, 2)  # One interval [0, 2]

model_config <- setup_cpp_model_pwc(
  V_0 = simple_data$V_0,
  V_healthy = simple_data$V_healthy,
  V_ill = simple_data$V_ill,
  T_obs = simple_data$T_obs,
  status_dead = simple_data$status_dead,
  status_ill = simple_data$status_ill,
  knots_12 = knots,
  knots_13 = knots,
  knots_23 = knots
)

# Test with known constant hazards
lambda_12 <- 0.3
lambda_13 <- 0.1
lambda_23 <- 0.5

ll <- calc_log_likelihood_pwc(
  model_config$model_pointer,
  lambda_12,
  lambda_13,
  lambda_23
)

cat("Log-likelihood with constant hazards:\n")
cat("  lambda_12 =", lambda_12, "\n")
cat("  lambda_13 =", lambda_13, "\n")
cat("  lambda_23 =", lambda_23, "\n")
cat("  LL =", ll, "\n\n")

# Manual calculation
# Case 1: L_1 = exp(A_12(0) + A_13(0)) * [exp(-A_12(1) - A_13(1)) + I_23(0, 1, 1)]
# 
# A_12(0) = 0
# A_13(0) = 0
# A_12(1) = 0.3 * 1 = 0.3
# A_13(1) = 0.1 * 1 = 0.1
# 
# First term: exp(0) * exp(-0.3 - 0.1) = exp(-0.4)
#
# For I_23(0, 1, 1):
# Build grid on [0, 1]: just {0, 1} (one interval)
# For interval [0, 1]:
#   left = 0, right = 1, delta = 1
#   H_0 = A_12(0) + A_13(0) - A_23(0) = 0
#   theta_0 = -lambda_12 - lambda_13 + lambda_23 = -0.3 - 0.1 + 0.5 = 0.1
#   phi(0.1, 1) = (exp(0.1) - 1) / 0.1 = (1.10517 - 1) / 0.1 = 1.0517
# 
# I_23 = exp(-A_23(1)) * lambda_23 * exp(-H_0) * phi(theta_0, delta)
#      = exp(-0.5 * 1) * 0.5 * exp(0) * 1.0517
#      = exp(-0.5) * 0.5 * 1.0517
#      = 0.6065 * 0.5 * 1.0517
#      = 0.319
#
# L_1 = 1 * (exp(-0.4) + 0.319)
#     = 0.6703 + 0.319
#     = 0.989
#
# log(L_1) = log(0.989) = -0.011

cat("Manual calculation:\n")
A12_0 <- 0
A13_0 <- 0
A12_1 <- lambda_12 * 1
A13_1 <- lambda_13 * 1
A23_1 <- lambda_23 * 1

term1 <- exp(-A12_1 - A13_1)
cat("  term1 = exp(-A_12(1) - A_13(1)) = exp(", -A12_1 - A13_1, ") =", term1, "\n")

# I_23 integral
H_0 <- 0
theta_0 <- -lambda_12 - lambda_13 + lambda_23
phi_val <- (exp(theta_0 * 1) - 1) / theta_0
I23_val <- exp(-A23_1) * lambda_23 * exp(-H_0) * phi_val

cat("  H_0 =", H_0, "\n")
cat("  theta_0 =", theta_0, "\n")
cat("  phi(theta_0, 1) =", phi_val, "\n")
cat("  I_23 = exp(-A_23(1)) * lambda_23 * exp(-H_0) * phi =", I23_val, "\n")

L1 <- exp(A12_0 + A13_0) * (term1 + I23_val)
log_L1 <- log(L1)

cat("  L_1 = exp(0) * (", term1, " + ", I23_val, ") =", L1, "\n")
cat("  log(L_1) =", log_L1, "\n\n")

cat("Difference: calculated - manual =", ll - log_L1, "\n")
