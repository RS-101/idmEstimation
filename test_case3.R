# Test Case 3 (illness)

library(Rcpp)
sourceCpp("src/piecewise_constant_hazard.cpp")
source("R/piecewise_constant_hazard.R")

# Create simple data: 1 person, got ill between V_healthy=0 and V_ill=0.5, still alive at T=1
simple_data <- data.frame(
  V_0 = 0,
  V_healthy = 0,
  V_ill = 0.5,
  T_obs = 1,
  status_dead = 0,
  status_ill = 1
)

cat("Simple data (1 person, got ill at V_ill=0.5, alive at T=1):\n")
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

# Manual calculation for Case 3
# L_3 = exp(A_12(V_0) + A_13(V_0) - A_23(T)) * I_12(V_k, V_{k+1}, T)
# where V_k = V_healthy = 0, V_{k+1} = V_ill = 0.5, T = T_obs = 1
#
# A_12(0) = 0
# A_13(0) = 0  
# A_23(1) = 0.5 * 1 = 0.5
#
# First term: exp(0 + 0 - 0.5) = exp(-0.5) = 0.6065
#
# For I_12(0, 0.5, 1):
# Build grid on [0, 0.5]: just {0, 0.5} (one interval)
# For interval [0, 0.5]:
#   left = 0, right = 0.5, delta = 0.5
#   H_0 = A_12(0) + A_13(0) - A_23(0) = 0
#   theta_0 = -lambda_12 - lambda_13 + lambda_23 = -0.3 - 0.1 + 0.5 = 0.1
#   phi(0.1, 0.5) = (exp(0.05) - 1) / 0.1 = (1.05127 - 1) / 0.1 = 0.5127
#
# I_12 = exp(-A_23(1)) * lambda_12 * exp(-H_0) * phi(theta_0, delta)
#      = exp(-0.5) * 0.3 * exp(0) * 0.5127
#      = 0.6065 * 0.3 * 0.5127
#      = 0.0933
#
# L_3 = exp(-0.5) * 0.0933
#     = 0.6065 * 0.0933
#     = 0.0566
#
# log(L_3) = log(0.0566) = -2.87

cat("Manual calculation:\n")
A12_0 <- 0
A13_0 <- 0
A23_1 <- lambda_23 * 1

cat("  exp(A_12(0) + A_13(0) - A_23(1)) = exp(", A12_0 + A13_0 - A23_1, ") =", exp(A12_0 + A13_0 - A23_1), "\n")

# I_12 integral from 0 to 0.5
H_0 <- 0
theta_0 <- -lambda_12 - lambda_13 + lambda_23
delta <- 0.5
phi_val <- (exp(theta_0 * delta) - 1) / theta_0
I12_val <- exp(-A23_1) * lambda_12 * exp(-H_0) * phi_val

cat("  H_0 =", H_0, "\n")
cat("  theta_0 =", theta_0, "\n")
cat("  delta =", delta, "\n")
cat("  phi(theta_0, delta) =", phi_val, "\n")
cat("  I_12 = exp(-A_23(1)) * lambda_12 * exp(-H_0) * phi =", I12_val, "\n")

L3 <- exp(A12_0 + A13_0 - A23_1) * I12_val
log_L3 <- log(L3)

cat("  L_3 = exp(-0.5) * ", I12_val, " =", L3, "\n")
cat("  log(L_3) =", log_L3, "\n\n")

cat("Difference: calculated - manual =", ll - log_L3, "\n")
