# Test Case 2 (died healthy)

library(Rcpp)
sourceCpp("src/piecewise_constant_hazard.cpp")
source("R/piecewise_constant_hazard.R")

# Create simple data: 1 person, died at T=1, was healthy (never got ill)
simple_data <- data.frame(
  V_0 = 0,
  V_healthy = 0,
  V_ill = NA,
  T_obs = 1,
  status_dead = 1,
  status_ill = 0
)

cat("Simple data (1 person, died healthy at T=1):\n")
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

# Manual calculation for Case 2
# L_2 = exp(A_12(V_0) + A_13(V_0)) * [exp(-A_12(T) - A_13(T)) * alpha_13(T) + alpha_23(T) * I_12(V_m, T, T)]
# where V_m = V_healthy = 0, T = T_obs = 1
#
# A_12(0) = 0
# A_13(0) = 0
# A_12(1) = 0.3
# A_13(1) = 0.1
# alpha_13(1) = 0.1 (constant)
# alpha_23(1) = 0.5 (constant)
#
# First term: exp(0) * [exp(-0.3 - 0.1) * 0.1 + 0.5 * I_12(0, 1, 1)]
#           = exp(-0.4) * 0.1 + 0.5 * I_12(0, 1, 1)
#           = 0.6703 * 0.1 + 0.5 * I_12(0, 1, 1)
#           = 0.06703 + 0.5 * I_12(0, 1, 1)
#
# For I_12(0, 1, 1):
# Build grid on [0, 1]: just {0, 1} (one interval)
# For interval [0, 1]:
#   left = 0, right = 1, delta = 1
#   H_0 = A_12(0) + A_13(0) - A_23(0) = 0
#   theta_0 = -lambda_12 - lambda_13 + lambda_23 = 0.1
#   phi(0.1, 1) = (exp(0.1) - 1) / 0.1 = 1.0517
#
# I_12 = exp(-A_23(1)) * lambda_12 * exp(-H_0) * phi(theta_0, delta)
#      = exp(-0.5) * 0.3 * exp(0) * 1.0517
#      = 0.6065 * 0.3 * 1.0517
#      = 0.1914
#
# L_2 = 0.06703 + 0.5 * 0.1914
#     = 0.06703 + 0.0957
#     = 0.1627
#
# log(L_2) = log(0.1627) = -1.815

cat("Manual calculation:\n")
A12_0 <- 0
A13_0 <- 0
A12_1 <- lambda_12 * 1
A13_1 <- lambda_13 * 1
A23_1 <- lambda_23 * 1
alpha13_1 <- lambda_13
alpha23_1 <- lambda_23

term1 <- exp(-A12_1 - A13_1) * alpha13_1
cat("  term1 = exp(-A_12(1) - A_13(1)) * alpha_13(1) =", term1, "\n")

# I_12 integral from 0 to 1
H_0 <- 0
theta_0 <- -lambda_12 - lambda_13 + lambda_23
delta <- 1
phi_val <- (exp(theta_0 * delta) - 1) / theta_0
I12_val <- exp(-A23_1) * lambda_12 * exp(-H_0) * phi_val

cat("  I_12 =", I12_val, "\n")

term2 <- alpha23_1 * I12_val
cat("  term2 = alpha_23(1) * I_12 =", term2, "\n")

L2 <- exp(A12_0 + A13_0) * (term1 + term2)
log_L2 <- log(L2)

cat("  L_2 = exp(0) * (", term1, " + ", term2, ") =", L2, "\n")
cat("  log(L_2) =", log_L2, "\n\n")

cat("Difference: calculated - manual =", ll - log_L2, "\n")
