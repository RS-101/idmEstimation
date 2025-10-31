# Test suite for penalized likelihood calculation
#
# These tests validate the C++ implementation of calc_case_1_log_likelihood
# by comparing against R implementations using both:
# 1. Constant/affine hazard rates (analytically tractable)
# 2. Spline-based numerical integration
#
# Case 1: Individual is healthy at last visit (T_obs) and was healthy at
#         previous visit (V_healthy)

test_that("calc_case_1_log_likelihood works with constant hazard rates", {
  # Test using affine cumulative hazards (constant hazard rates)
  # This allows analytical computation for validation

  # Set seed for reproducibility
  set.seed(42)

  # Define constant hazard rates
  alpha_12 <- 0.03  # constant hazard 1->2
  alpha_13 <- 0.02  # constant hazard 1->3
  alpha_23 <- 0.05  # constant hazard 2->3

  # Generate test data
  n <- 10
  V_healthy <- sort(runif(n, 0, 3))  # healthy visit times
  T_obs <- sort(runif(n, 3, 8))      # observation times
  V_0 <- rep(0, n)                   # initial times

  # For constant hazards, cumulative hazard is A(t) = alpha * t
  # We'll use splines that approximate this

  # Create knots for splines (need enough to represent linear function)
  knots_12 <- seq(0, 10, length.out = 10)
  knots_13 <- seq(0, 10, length.out = 10)
  knots_23 <- seq(0, 10, length.out = 10)

  # Create spline matrices
  spline_V_healthy_12 <- make_spline_mat(V_healthy, knots_12)
  spline_V_healthy_13 <- make_spline_mat(V_healthy, knots_13)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23)

  spline_V_0_12 <- make_spline_mat(V_0, knots_12)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13)

  # Create grid for integration
  grid_T_obs <- seq(min(V_healthy), max(T_obs), length.out = 250)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23)

  # Create grid for V_ill (needed for structure)
  grid_V_ill <- seq(min(V_healthy), max(T_obs), length.out = 250)
  dx_grid_V_ill <- diff(grid_V_ill)[1]

  spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots_12)
  spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots_13)
  spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots_23)

  # Fit theta parameters to represent constant hazards
  # For constant hazard alpha, A(t) = alpha * t
  # We need to find theta such that I-spline %*% theta ≈ alpha * t

  # Use least squares to fit
  fit_constant_hazard <- function(spline_mat, times, alpha) {
    target <- alpha * times
    theta <- as.vector(solve(t(spline_mat$i_spline) %*% spline_mat$i_spline +
                               diag(1e-6, ncol(spline_mat$i_spline))) %*%
                        t(spline_mat$i_spline) %*% target)
    # Ensure non-negative
    theta <- pmax(theta, 0)
    return(theta)
  }

  # Fit times spanning the range
  fit_times <- seq(0, 10, length.out = 100)
  fit_spline_12 <- make_spline_mat(fit_times, knots_12)
  fit_spline_13 <- make_spline_mat(fit_times, knots_13)
  fit_spline_23 <- make_spline_mat(fit_times, knots_23)

  theta_12 <- fit_constant_hazard(fit_spline_12, fit_times, alpha_12)
  theta_13 <- fit_constant_hazard(fit_spline_13, fit_times, alpha_13)
  theta_23 <- fit_constant_hazard(fit_spline_23, fit_times, alpha_23)

  # Create model data
  model_data <- list(
    V_healthy_i_spline_mat_12 = spline_V_healthy_12$i_spline,
    V_healthy_i_spline_mat_13 = spline_V_healthy_13$i_spline,

    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,

    V_0_i_spline_mat_12 = spline_V_0_12$i_spline,
    V_0_i_spline_mat_13 = spline_V_0_13$i_spline,

    grid_T_obs_i_spline_mat_12 = spline_grid_T_obs_12$i_spline,
    grid_T_obs_i_spline_mat_13 = spline_grid_T_obs_13$i_spline,
    grid_T_obs_i_spline_mat_23 = spline_grid_T_obs_23$i_spline,

    grid_V_ill_i_spline_mat_12 = spline_grid_V_ill_12$i_spline,
    grid_V_ill_i_spline_mat_13 = spline_grid_V_ill_13$i_spline,
    grid_V_ill_i_spline_mat_23 = spline_grid_V_ill_23$i_spline,

    T_obs_m_spline_mat_12 = spline_T_obs_12$m_spline,
    T_obs_m_spline_mat_13 = spline_T_obs_13$m_spline,
    T_obs_m_spline_mat_23 = spline_T_obs_23$m_spline,

    dx_grid_T_obs = dx_grid_T_obs,
    dx_grid_V_ill = dx_grid_V_ill,

    grid_T_obs_m_spline_mat_12 = spline_grid_T_obs_12$m_spline,
    grid_T_obs_m_spline_mat_13 = spline_grid_T_obs_13$m_spline,
    grid_T_obs_m_spline_mat_23 = spline_grid_T_obs_23$m_spline,

    grid_V_ill_m_spline_mat_12 = spline_grid_V_ill_12$m_spline,
    grid_V_ill_m_spline_mat_13 = spline_grid_V_ill_13$m_spline,
    grid_V_ill_m_spline_mat_23 = spline_grid_V_ill_23$m_spline,

    T_obs_values = T_obs,
    V_healthy_values = V_healthy
  )

  # Create pointer to model data
  md_ptr <- create_penlik_model_data(model_data)

  # Calculate log-likelihood using C++
  log_lik_cpp <- calc_case_1_log_likelihood(md_ptr, theta_12, theta_13, theta_23)

  # Calculate log-likelihood using R (matching C++ implementation)
  A12_V_0 <- spline_V_0_12$i_spline %*% theta_12
  A13_V_0 <- spline_V_0_13$i_spline %*% theta_13
  factored_out <- exp(A12_V_0 + A13_V_0)

  A12_T_obs <- spline_T_obs_12$i_spline %*% theta_12
  A13_T_obs <- spline_T_obs_13$i_spline %*% theta_13
  term_1 <- exp(-(A12_T_obs + A13_T_obs))

  factored_out_term_2 <- exp(-(spline_T_obs_23$i_spline %*% theta_23))

  # Compute integrand
  integrand <- exp(
    -(spline_grid_T_obs_12$i_spline %*% theta_12) -
    (spline_grid_T_obs_13$i_spline %*% theta_13)
  ) *
  (spline_grid_T_obs_12$m_spline %*% theta_12) *
  exp(spline_grid_T_obs_23$i_spline %*% theta_23)

  # Compute midpoints for trapezoidal rule
  mid <- (integrand[-1] + integrand[-length(integrand)]) / 2

  # Compute cumulative integral
  integral <- c(0, cumsum(mid * dx_grid_T_obs))

  # Get indices
  grid_min <- min(grid_T_obs)
  T_obs_indices <- round((T_obs - grid_min) / dx_grid_T_obs) + 1
  V_healthy_indices <- round((V_healthy - grid_min) / dx_grid_T_obs) + 1

  # Calculate term_2
  term_2 <- factored_out_term_2 * (
    integral[T_obs_indices] - integral[V_healthy_indices]
  )

  # Calculate final log-likelihood
  log_lik_r <- sum(log(factored_out * (term_1 + term_2)))

  # Test that C++ and R implementations match
  expect_equal(log_lik_cpp, log_lik_r, tolerance = 1e-10)

  # Test that result is finite
  expect_true(is.finite(log_lik_cpp))

  # Test that result is a single number
  expect_length(log_lik_cpp, 1)
})


test_that("calc_case_1_log_likelihood analytical validation for simple case", {
  # Test with a very simple case where we can verify the result analytically
  # Use a single observation with constant hazards

  set.seed(123)

  # Single observation
  V_healthy <- 1.0
  T_obs <- 3.0
  V_0 <- 0.0

  # Constant hazards
  alpha_12 <- 0.1
  alpha_13 <- 0.05
  alpha_23 <- 0.08

  # Create knots
  knots <- seq(0, 5, length.out = 15)

  # Create splines
  spline_V_healthy_12 <- make_spline_mat(V_healthy, knots)
  spline_V_healthy_13 <- make_spline_mat(V_healthy, knots)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots)

  spline_V_0_12 <- make_spline_mat(V_0, knots)
  spline_V_0_13 <- make_spline_mat(V_0, knots)

  # Fine grid for integration
  grid_T_obs <- seq(V_healthy, T_obs, length.out = 500)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots)

  # Dummy grid for V_ill
  grid_V_ill <- seq(V_healthy, T_obs, length.out = 500)
  dx_grid_V_ill <- diff(grid_V_ill)[1]

  spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots)
  spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots)
  spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots)

  # Fit constant hazards
  fit_times <- seq(0, 5, length.out = 100)
  fit_spline <- make_spline_mat(fit_times, knots)

  fit_constant <- function(spline_mat, times, alpha) {
    target <- alpha * times
    theta <- as.vector(solve(t(spline_mat$i_spline) %*% spline_mat$i_spline +
                               diag(1e-6, ncol(spline_mat$i_spline))) %*%
                        t(spline_mat$i_spline) %*% target)
    pmax(theta, 0)
  }

  theta_12 <- fit_constant(fit_spline, fit_times, alpha_12)
  theta_13 <- fit_constant(fit_spline, fit_times, alpha_13)
  theta_23 <- fit_constant(fit_spline, fit_times, alpha_23)

  # Create model data
  model_data <- list(
    V_healthy_i_spline_mat_12 = spline_V_healthy_12$i_spline,
    V_healthy_i_spline_mat_13 = spline_V_healthy_13$i_spline,

    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,

    V_0_i_spline_mat_12 = spline_V_0_12$i_spline,
    V_0_i_spline_mat_13 = spline_V_0_13$i_spline,

    grid_T_obs_i_spline_mat_12 = spline_grid_T_obs_12$i_spline,
    grid_T_obs_i_spline_mat_13 = spline_grid_T_obs_13$i_spline,
    grid_T_obs_i_spline_mat_23 = spline_grid_T_obs_23$i_spline,

    grid_V_ill_i_spline_mat_12 = spline_grid_V_ill_12$i_spline,
    grid_V_ill_i_spline_mat_13 = spline_grid_V_ill_13$i_spline,
    grid_V_ill_i_spline_mat_23 = spline_grid_V_ill_23$i_spline,

    T_obs_m_spline_mat_12 = spline_T_obs_12$m_spline,
    T_obs_m_spline_mat_13 = spline_T_obs_13$m_spline,
    T_obs_m_spline_mat_23 = spline_T_obs_23$m_spline,

    dx_grid_T_obs = dx_grid_T_obs,
    dx_grid_V_ill = dx_grid_V_ill,

    grid_T_obs_m_spline_mat_12 = spline_grid_T_obs_12$m_spline,
    grid_T_obs_m_spline_mat_13 = spline_grid_T_obs_13$m_spline,
    grid_T_obs_m_spline_mat_23 = spline_grid_T_obs_23$m_spline,

    grid_V_ill_m_spline_mat_12 = spline_grid_V_ill_12$m_spline,
    grid_V_ill_m_spline_mat_13 = spline_grid_V_ill_13$m_spline,
    grid_V_ill_m_spline_mat_23 = spline_grid_V_ill_23$m_spline,

    T_obs_values = T_obs,
    V_healthy_values = V_healthy
  )

  md_ptr <- create_penlik_model_data(model_data)
  log_lik_cpp <- calc_case_1_log_likelihood(md_ptr, theta_12, theta_13, theta_23)

  # Analytical formula for constant hazards
  # L = exp(A12(V_0) + A13(V_0)) *
  #     (exp(-(alpha_12 + alpha_13) * T) +
  #      exp(-alpha_23 * T) * integral from V_healthy to T of
  #      exp(-(alpha_12 + alpha_13) * t) * alpha_12 * exp(alpha_23 * t) dt)

  # Simplify: exp(-(alpha_12 + alpha_13) * t) * exp(alpha_23 * t) =
  #           exp(-(alpha_12 + alpha_13 - alpha_23) * t)

  k <- alpha_12 + alpha_13 - alpha_23

  # Analytical integral
  if (abs(k) > 1e-10) {
    integral_analytical <- (alpha_12 / k) *
      (exp(-k * V_healthy) - exp(-k * T_obs))
  } else {
    # If k ≈ 0, use limiting case
    integral_analytical <- alpha_12 * (T_obs - V_healthy)
  }

  A12_V_0 <- alpha_12 * V_0
  A13_V_0 <- alpha_13 * V_0
  term_0_analytical <- exp(A12_V_0 + A13_V_0)
  term_1_analytical <- exp(-(alpha_12 + alpha_13) * T_obs)
  term_2_analytical <- exp(-alpha_23 * T_obs) * integral_analytical

  log_lik_analytical <- log(term_0_analytical * (term_1_analytical + term_2_analytical))

  # Test that spline approximation is close to analytical
  # Allow larger tolerance due to spline approximation
  expect_equal(log_lik_cpp, log_lik_analytical, tolerance = 0.01)
})


# Case 2 Tests ----------------------------------------------------------------

test_that("calc_case_2_log_likelihood works with constant hazard rates", {
  # Case 2: Individual died (state 3) at T_obs, was healthy at V_healthy

  set.seed(100)

  # Constant hazard rates
  alpha_12 <- 0.04
  alpha_13 <- 0.03
  alpha_23 <- 0.06

  # Generate test data
  n <- 10
  V_healthy <- sort(runif(n, 0, 2))
  T_obs <- sort(runif(n, 2, 6))
  V_0 <- rep(0, n)

  # Create knots
  knots_12 <- seq(0, 8, length.out = 10)
  knots_13 <- seq(0, 8, length.out = 10)
  knots_23 <- seq(0, 8, length.out = 10)

  # Create splines
  spline_V_healthy_12 <- make_spline_mat(V_healthy, knots_12)
  spline_V_healthy_13 <- make_spline_mat(V_healthy, knots_13)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23)

  spline_V_0_12 <- make_spline_mat(V_0, knots_12)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13)

  # Grid for integration
  grid_T_obs <- seq(min(V_healthy), max(T_obs), length.out = 250)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23)

  grid_V_ill <- grid_T_obs
  dx_grid_V_ill <- dx_grid_T_obs
  spline_grid_V_ill_12 <- spline_grid_T_obs_12
  spline_grid_V_ill_13 <- spline_grid_T_obs_13
  spline_grid_V_ill_23 <- spline_grid_T_obs_23

  # Fit theta parameters
  fit_times <- seq(0, 8, length.out = 100)
  fit_spline_12 <- make_spline_mat(fit_times, knots_12)
  fit_spline_13 <- make_spline_mat(fit_times, knots_13)
  fit_spline_23 <- make_spline_mat(fit_times, knots_23)

  fit_constant_hazard <- function(spline_mat, times, alpha) {
    target <- alpha * times
    theta <- as.vector(solve(t(spline_mat$i_spline) %*% spline_mat$i_spline +
                               diag(1e-6, ncol(spline_mat$i_spline))) %*%
                        t(spline_mat$i_spline) %*% target)
    pmax(theta, 0)
  }

  theta_12 <- fit_constant_hazard(fit_spline_12, fit_times, alpha_12)
  theta_13 <- fit_constant_hazard(fit_spline_13, fit_times, alpha_13)
  theta_23 <- fit_constant_hazard(fit_spline_23, fit_times, alpha_23)

  # Create model data
  model_data <- list(
    V_healthy_i_spline_mat_12 = spline_V_healthy_12$i_spline,
    V_healthy_i_spline_mat_13 = spline_V_healthy_13$i_spline,

    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,

    V_0_i_spline_mat_12 = spline_V_0_12$i_spline,
    V_0_i_spline_mat_13 = spline_V_0_13$i_spline,

    grid_T_obs_i_spline_mat_12 = spline_grid_T_obs_12$i_spline,
    grid_T_obs_i_spline_mat_13 = spline_grid_T_obs_13$i_spline,
    grid_T_obs_i_spline_mat_23 = spline_grid_T_obs_23$i_spline,

    grid_V_ill_i_spline_mat_12 = spline_grid_V_ill_12$i_spline,
    grid_V_ill_i_spline_mat_13 = spline_grid_V_ill_13$i_spline,
    grid_V_ill_i_spline_mat_23 = spline_grid_V_ill_23$i_spline,

    T_obs_m_spline_mat_12 = spline_T_obs_12$m_spline,
    T_obs_m_spline_mat_13 = spline_T_obs_13$m_spline,
    T_obs_m_spline_mat_23 = spline_T_obs_23$m_spline,

    dx_grid_T_obs = dx_grid_T_obs,
    dx_grid_V_ill = dx_grid_V_ill,

    grid_T_obs_m_spline_mat_12 = spline_grid_T_obs_12$m_spline,
    grid_T_obs_m_spline_mat_13 = spline_grid_T_obs_13$m_spline,
    grid_T_obs_m_spline_mat_23 = spline_grid_T_obs_23$m_spline,

    grid_V_ill_m_spline_mat_12 = spline_grid_V_ill_12$m_spline,
    grid_V_ill_m_spline_mat_13 = spline_grid_V_ill_13$m_spline,
    grid_V_ill_m_spline_mat_23 = spline_grid_V_ill_23$m_spline,

    T_obs_values = T_obs,
    V_healthy_values = V_healthy
  )

  md_ptr <- create_penlik_model_data(model_data)
  log_lik_cpp <- calc_case_2_log_likelihood(md_ptr, theta_12, theta_13, theta_23)

  # Calculate in R for comparison
  A12_V_0 <- spline_V_0_12$i_spline %*% theta_12
  A13_V_0 <- spline_V_0_13$i_spline %*% theta_13
  factored_out <- exp(A12_V_0 + A13_V_0)

  A12_T_obs <- spline_T_obs_12$i_spline %*% theta_12
  A13_T_obs <- spline_T_obs_13$i_spline %*% theta_13
  a13_T_obs <- spline_T_obs_13$m_spline %*% theta_13

  A23_T_obs <- spline_T_obs_23$i_spline %*% theta_23
  a23_T_obs <- spline_T_obs_23$m_spline %*% theta_23

  # Term 1: direct 1->3
  term_1 <- exp(-(A12_T_obs + A13_T_obs)) * a13_T_obs

  # Term 2: via 1->2->3
  factored_out_term_2 <- exp(-A23_T_obs) * a23_T_obs

  integrand <- exp(
    -(spline_grid_T_obs_12$i_spline %*% theta_12) -
    (spline_grid_T_obs_13$i_spline %*% theta_13)
  ) *
  (spline_grid_T_obs_12$m_spline %*% theta_12) *
  exp(spline_grid_T_obs_23$i_spline %*% theta_23)

  mid <- (integrand[-1] + integrand[-length(integrand)]) / 2
  integral <- c(0, cumsum(mid * dx_grid_T_obs))

  grid_min <- min(grid_T_obs)
  T_obs_indices <- round((T_obs - grid_min) / dx_grid_T_obs) + 1
  V_healthy_indices <- round((V_healthy - grid_min) / dx_grid_T_obs) + 1

  term_2 <- factored_out_term_2 * (integral[T_obs_indices] - integral[V_healthy_indices])

  log_lik_r <- sum(log(factored_out * (term_1 + term_2)))

  # Test C++ vs R
  expect_equal(log_lik_cpp, log_lik_r, tolerance = 1e-10)
  expect_true(is.finite(log_lik_cpp))
})


# Case 3 Tests ----------------------------------------------------------------

test_that("calc_case_3_log_likelihood works with constant hazard rates", {
  # Case 3: Subject is healthy at V_k (V_healthy), ill at V_{k+1} (V_ill), and still alive at T

  set.seed(300)

  # Constant hazard rates
  alpha_12 <- 0.05
  alpha_13 <- 0.02
  alpha_23 <- 0.04

  # Generate test data
  n <- 10
  V_healthy <- sort(runif(n, 0, 2))   # V_k - when healthy
  V_ill <- sort(runif(n, 2, 4))       # V_{k+1} - when ill
  T_obs <- sort(runif(n, 4, 6))       # T - observation time (still alive)
  V_0 <- rep(0, n)                    # initial time

  # Create knots
  knots_12 <- seq(0, 7, length.out = 10)
  knots_13 <- seq(0, 7, length.out = 10)
  knots_23 <- seq(0, 7, length.out = 10)

  # Create splines for all time points
  spline_V_healthy_12 <- make_spline_mat(V_healthy, knots_12)
  spline_V_healthy_13 <- make_spline_mat(V_healthy, knots_13)

  spline_V_0_12 <- make_spline_mat(V_0, knots_12)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13)

  spline_V_ill_12 <- make_spline_mat(V_ill, knots_12)
  spline_V_ill_13 <- make_spline_mat(V_ill, knots_13)
  spline_V_ill_23 <- make_spline_mat(V_ill, knots_23)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23)

  # Grid for integration between V_healthy and V_ill
  grid_V_ill <- seq(min(V_healthy), max(V_ill), length.out = 250)
  dx_grid_V_ill <- diff(grid_V_ill)[1]

  spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots_12)
  spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots_13)
  spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots_23)

  # Dummy grid for T_obs (required by struct)
  grid_T_obs <- seq(0, 7, length.out = 100)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23)

  # Fit theta parameters
  fit_times <- seq(0, 7, length.out = 100)
  fit_spline_12 <- make_spline_mat(fit_times, knots_12)
  fit_spline_13 <- make_spline_mat(fit_times, knots_13)
  fit_spline_23 <- make_spline_mat(fit_times, knots_23)

  fit_constant_hazard <- function(spline_mat, times, alpha) {
    target <- alpha * times
    theta <- as.vector(solve(t(spline_mat$i_spline) %*% spline_mat$i_spline +
                               diag(1e-6, ncol(spline_mat$i_spline))) %*%
                        t(spline_mat$i_spline) %*% target)
    pmax(theta, 0)
  }

  theta_12 <- fit_constant_hazard(fit_spline_12, fit_times, alpha_12)
  theta_13 <- fit_constant_hazard(fit_spline_13, fit_times, alpha_13)
  theta_23 <- fit_constant_hazard(fit_spline_23, fit_times, alpha_23)

  # Create model data
  model_data <- list(
    V_healthy_i_spline_mat_12 = spline_V_healthy_12$i_spline,
    V_healthy_i_spline_mat_13 = spline_V_healthy_13$i_spline,

    V_ill_i_spline_mat_12 = spline_V_ill_12$i_spline,
    V_ill_i_spline_mat_13 = spline_V_ill_13$i_spline,
    V_ill_i_spline_mat_23 = spline_V_ill_23$i_spline,
    V_ill_values = V_ill,

    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,

    V_0_i_spline_mat_12 = spline_V_0_12$i_spline,
    V_0_i_spline_mat_13 = spline_V_0_13$i_spline,

    grid_T_obs_i_spline_mat_12 = spline_grid_T_obs_12$i_spline,
    grid_T_obs_i_spline_mat_13 = spline_grid_T_obs_13$i_spline,
    grid_T_obs_i_spline_mat_23 = spline_grid_T_obs_23$i_spline,

    grid_V_ill_i_spline_mat_12 = spline_grid_V_ill_12$i_spline,
    grid_V_ill_i_spline_mat_13 = spline_grid_V_ill_13$i_spline,
    grid_V_ill_i_spline_mat_23 = spline_grid_V_ill_23$i_spline,

    T_obs_m_spline_mat_12 = spline_T_obs_12$m_spline,
    T_obs_m_spline_mat_13 = spline_T_obs_13$m_spline,
    T_obs_m_spline_mat_23 = spline_T_obs_23$m_spline,

    dx_grid_T_obs = dx_grid_T_obs,
    dx_grid_V_ill = dx_grid_V_ill,

    grid_T_obs_m_spline_mat_12 = spline_grid_T_obs_12$m_spline,
    grid_T_obs_m_spline_mat_13 = spline_grid_T_obs_13$m_spline,
    grid_T_obs_m_spline_mat_23 = spline_grid_T_obs_23$m_spline,

    grid_V_ill_m_spline_mat_12 = spline_grid_V_ill_12$m_spline,
    grid_V_ill_m_spline_mat_13 = spline_grid_V_ill_13$m_spline,
    grid_V_ill_m_spline_mat_23 = spline_grid_V_ill_23$m_spline,

    T_obs_values = T_obs,
    V_healthy_values = V_healthy
  )

  md_ptr <- create_penlik_model_data(model_data)
  log_lik_cpp <- calc_case_3_log_likelihood(md_ptr, theta_12, theta_13, theta_23)

  # Calculate in R for comparison
  # L = exp(A12(V_0) + A13(V_0)) * exp(-A23(T)) *
  #     ∫_{V_healthy}^{V_ill} exp(-A12(u) - A13(u) + A23(u)) * α12(u) du

  A12_V_0 <- spline_V_0_12$i_spline %*% theta_12
  A13_V_0 <- spline_V_0_13$i_spline %*% theta_13
  factored_out <- exp(A12_V_0 + A13_V_0)

  A23_T_obs <- spline_T_obs_23$i_spline %*% theta_23
  exp_neg_A23_T <- exp(-A23_T_obs)

  # Integrand
  integrand <- exp(
    -(spline_grid_V_ill_12$i_spline %*% theta_12) -
    (spline_grid_V_ill_13$i_spline %*% theta_13) +
    (spline_grid_V_ill_23$i_spline %*% theta_23)
  ) * (spline_grid_V_ill_12$m_spline %*% theta_12)

  # Trapezoidal rule
  mid <- (integrand[-1] + integrand[-length(integrand)]) / 2
  integral <- c(0, cumsum(mid * dx_grid_V_ill))

  # Get indices
  grid_min <- min(grid_V_ill)
  V_healthy_indices <- round((V_healthy - grid_min) / dx_grid_V_ill) + 1
  V_ill_indices <- round((V_ill - grid_min) / dx_grid_V_ill) + 1

  integral_value <- integral[V_ill_indices] - integral[V_healthy_indices]

  likelihood_r <- factored_out * exp_neg_A23_T * integral_value
  log_lik_r <- sum(log(likelihood_r))

  # Test C++ vs R
  expect_equal(log_lik_cpp, log_lik_r, tolerance = 1e-10)
  expect_true(is.finite(log_lik_cpp))
  expect_length(log_lik_cpp, 1)
})


# Case 4 Tests ----------------------------------------------------------------

test_that("calc_case_4_log_likelihood works with constant hazard rates", {
  # Case 4: Subject is healthy at V_k (V_healthy), ill at V_{k+1} (V_ill), and dies at T

  set.seed(400)

  # Constant hazard rates
  alpha_12 <- 0.06
  alpha_13 <- 0.01
  alpha_23 <- 0.08

  # Generate test data
  n <- 8
  V_healthy <- sort(runif(n, 0, 1.5))   # V_k - when healthy
  V_ill <- sort(runif(n, 1.5, 3))       # V_{k+1} - when ill
  T_obs <- sort(runif(n, 3, 5))         # T - death time
  V_0 <- rep(0, n)                      # initial time

  # Create knots
  knots_12 <- seq(0, 6, length.out = 10)
  knots_13 <- seq(0, 6, length.out = 10)
  knots_23 <- seq(0, 6, length.out = 10)

  # Create splines
  spline_V_healthy_12 <- make_spline_mat(V_healthy, knots_12)
  spline_V_healthy_13 <- make_spline_mat(V_healthy, knots_13)

  spline_V_0_12 <- make_spline_mat(V_0, knots_12)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13)

  spline_V_ill_12 <- make_spline_mat(V_ill, knots_12)
  spline_V_ill_13 <- make_spline_mat(V_ill, knots_13)
  spline_V_ill_23 <- make_spline_mat(V_ill, knots_23)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23)

  # Grid for integration between V_healthy and V_ill
  grid_V_ill <- seq(min(V_healthy), max(V_ill), length.out = 250)
  dx_grid_V_ill <- diff(grid_V_ill)[1]

  spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots_12)
  spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots_13)
  spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots_23)

  # Dummy grid for T_obs
  grid_T_obs <- seq(0, 6, length.out = 100)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23)

  # Fit theta parameters
  fit_times <- seq(0, 6, length.out = 100)
  fit_spline_12 <- make_spline_mat(fit_times, knots_12)
  fit_spline_13 <- make_spline_mat(fit_times, knots_13)
  fit_spline_23 <- make_spline_mat(fit_times, knots_23)

  fit_constant_hazard <- function(spline_mat, times, alpha) {
    target <- alpha * times
    theta <- as.vector(solve(t(spline_mat$i_spline) %*% spline_mat$i_spline +
                               diag(1e-6, ncol(spline_mat$i_spline))) %*%
                        t(spline_mat$i_spline) %*% target)
    pmax(theta, 0)
  }

  theta_12 <- fit_constant_hazard(fit_spline_12, fit_times, alpha_12)
  theta_13 <- fit_constant_hazard(fit_spline_13, fit_times, alpha_13)
  theta_23 <- fit_constant_hazard(fit_spline_23, fit_times, alpha_23)

  # Create model data
  model_data <- list(
    V_healthy_i_spline_mat_12 = spline_V_healthy_12$i_spline,
    V_healthy_i_spline_mat_13 = spline_V_healthy_13$i_spline,

    V_ill_i_spline_mat_12 = spline_V_ill_12$i_spline,
    V_ill_i_spline_mat_13 = spline_V_ill_13$i_spline,
    V_ill_i_spline_mat_23 = spline_V_ill_23$i_spline,
    V_ill_values = V_ill,

    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,

    V_0_i_spline_mat_12 = spline_V_0_12$i_spline,
    V_0_i_spline_mat_13 = spline_V_0_13$i_spline,

    grid_T_obs_i_spline_mat_12 = spline_grid_T_obs_12$i_spline,
    grid_T_obs_i_spline_mat_13 = spline_grid_T_obs_13$i_spline,
    grid_T_obs_i_spline_mat_23 = spline_grid_T_obs_23$i_spline,

    grid_V_ill_i_spline_mat_12 = spline_grid_V_ill_12$i_spline,
    grid_V_ill_i_spline_mat_13 = spline_grid_V_ill_13$i_spline,
    grid_V_ill_i_spline_mat_23 = spline_grid_V_ill_23$i_spline,

    T_obs_m_spline_mat_12 = spline_T_obs_12$m_spline,
    T_obs_m_spline_mat_13 = spline_T_obs_13$m_spline,
    T_obs_m_spline_mat_23 = spline_T_obs_23$m_spline,

    dx_grid_T_obs = dx_grid_T_obs,
    dx_grid_V_ill = dx_grid_V_ill,

    grid_T_obs_m_spline_mat_12 = spline_grid_T_obs_12$m_spline,
    grid_T_obs_m_spline_mat_13 = spline_grid_T_obs_13$m_spline,
    grid_T_obs_m_spline_mat_23 = spline_grid_T_obs_23$m_spline,

    grid_V_ill_m_spline_mat_12 = spline_grid_V_ill_12$m_spline,
    grid_V_ill_m_spline_mat_13 = spline_grid_V_ill_13$m_spline,
    grid_V_ill_m_spline_mat_23 = spline_grid_V_ill_23$m_spline,

    T_obs_values = T_obs,
    V_healthy_values = V_healthy
  )

  md_ptr <- create_penlik_model_data(model_data)
  log_lik_cpp <- calc_case_4_log_likelihood(md_ptr, theta_12, theta_13, theta_23)

  # Calculate in R for comparison
  # L = exp(A12(V_0) + A13(V_0)) * exp(-A23(T)) * α23(T) *
  #     ∫_{V_healthy}^{V_ill} exp(-A12(u) - A13(u) + A23(u)) * α12(u) du

  A12_V_0 <- spline_V_0_12$i_spline %*% theta_12
  A13_V_0 <- spline_V_0_13$i_spline %*% theta_13
  factored_out <- exp(A12_V_0 + A13_V_0)

  A23_T_obs <- spline_T_obs_23$i_spline %*% theta_23
  a23_T_obs <- spline_T_obs_23$m_spline %*% theta_23
  exp_neg_A23_T_times_a23 <- exp(-A23_T_obs) * a23_T_obs

  # Integrand
  integrand <- exp(
    -(spline_grid_V_ill_12$i_spline %*% theta_12) -
    (spline_grid_V_ill_13$i_spline %*% theta_13) +
    (spline_grid_V_ill_23$i_spline %*% theta_23)
  ) * (spline_grid_V_ill_12$m_spline %*% theta_12)

  # Trapezoidal rule
  mid <- (integrand[-1] + integrand[-length(integrand)]) / 2
  integral <- c(0, cumsum(mid * dx_grid_V_ill))

  # Get indices
  grid_min <- min(grid_V_ill)
  V_healthy_indices <- round((V_healthy - grid_min) / dx_grid_V_ill) + 1
  V_ill_indices <- round((V_ill - grid_min) / dx_grid_V_ill) + 1

  integral_value <- integral[V_ill_indices] - integral[V_healthy_indices]

  likelihood_r <- factored_out * exp_neg_A23_T_times_a23 * integral_value
  log_lik_r <- sum(log(likelihood_r))

  # Test C++ vs R
  expect_equal(log_lik_cpp, log_lik_r, tolerance = 1e-10)
  expect_true(is.finite(log_lik_cpp))
  expect_length(log_lik_cpp, 1)
})
