library(idmEstimation)

# Generate some test data
set.seed(123)
V_healthy <- sort(runif(10, 0, 5))  # healthy visit times
T_obs <- sort(runif(10, 5, 10))     # observation times
V_0 <- rep(0, 10)                   # initial times

# Create knots for splines
knots_12 <- seq(0, 10, length.out = 5)  # knots for transition 1->2
knots_13 <- seq(0, 10, length.out = 5)  # knots for transition 1->3
knots_23 <- seq(0, 10, length.out = 5)  # knots for transition 2->3

# Create spline matrices for all time points
spline_V_healthy_12 <- make_spline_mat(V_healthy, knots_12)
spline_V_healthy_13 <- make_spline_mat(V_healthy, knots_13)

spline_T_obs_12 <- make_spline_mat(T_obs, knots_12)
spline_T_obs_13 <- make_spline_mat(T_obs, knots_13)
spline_T_obs_23 <- make_spline_mat(T_obs, knots_23)

spline_V_0_12 <- make_spline_mat(V_0, knots_12)
spline_V_0_13 <- make_spline_mat(V_0, knots_13)

# Create grid points for numerical integration
grid_T_obs <- seq(min(V_healthy), max(T_obs), length.out = 250)
dx_grid_T_obs <- diff(grid_T_obs)[1]  # assuming uniform grid

spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12)
spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13)
spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23)

# Create grid points for ill visits (not used in case 1 but needed for the structure)
grid_V_ill <- seq(min(V_healthy), max(T_obs), length.out = 250)
dx_grid_V_ill <- diff(grid_V_ill)[1]

spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots_12)
spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots_13)
spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots_23)

# Create the list of matrices for the model
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
    
    # Add actual time values
    T_obs_values = T_obs,
    V_healthy_values = V_healthy
)

# Create the pointer to the model data
md_ptr <- create_penlik_model_data(model_data)

# Create some test parameters (random values between 0 and 1)
theta_12 <- runif(ncol(spline_T_obs_12$i_spline))
theta_13 <- runif(ncol(spline_T_obs_13$i_spline))
theta_23 <- runif(ncol(spline_T_obs_23$i_spline))

# Calculate the log-likelihood using C++
log_lik_cpp <- calc_case_1_log_likelihood(md_ptr, theta_12, theta_13, theta_23)

# Calculate the same log-likelihood in R for comparison
# First calculate A12(T_obs) and A13(T_obs)
A12_T_obs <- spline_T_obs_12$i_spline %*% theta_12
A13_T_obs <- spline_T_obs_13$i_spline %*% theta_13
A23_T_obs <- spline_T_obs_23$i_spline %*% theta_23

# Calculate term1
term_1 <- exp(-A12_T_obs - A13_T_obs)

# Calculate factored_out term
factored_out <- exp(-A23_T_obs)

# Create grid and calculate integrand
grid <- grid_T_obs
A12_grid <- spline_grid_T_obs_12$i_spline %*% theta_12
A13_grid <- spline_grid_T_obs_13$i_spline %*% theta_13
A23_grid <- spline_grid_T_obs_23$i_spline %*% theta_23
a12_grid <- spline_grid_T_obs_12$m_spline %*% theta_12

integrand <- exp(-A12_grid - A13_grid) * a12_grid * exp(A23_grid)

# Calculate midpoints for trapezoidal rule
dx <- dx_grid_T_obs
mid <- (integrand[-1] + integrand[-length(integrand)]) / 2
integral_at_grid <- c(0, cumsum(dx * mid))

# Create interpolation function
integral_fun <- splinefun(grid, integral_at_grid, method = "natural")

# Calculate term2
A12_V_healthy <- spline_V_healthy_12$i_spline %*% theta_12
term_2 <- factored_out * (integral_fun(T_obs) - integral_fun(V_healthy))

# Calculate final log-likelihood
log_lik_r <- sum(log(term_1 + term_2))

# Debug printing function
print_vec <- function(name, vec) {
    cat(sprintf("%s: First 3 values: %s\n", name,
        paste(format(head(vec, 3), digits=6), collapse=", ")))
}

# R implementation matching C++ exactly
# Calculate the same quantities but using the C++ approach
A12_T_obs <- spline_T_obs_12$i_spline %*% theta_12
A13_T_obs <- spline_T_obs_13$i_spline %*% theta_13
term_1_direct <- exp(-(A12_T_obs + A13_T_obs))

print_vec("A12_T_obs", A12_T_obs)
print_vec("A13_T_obs", A13_T_obs)
print_vec("term_1_direct", term_1_direct)

factored_out_direct <- exp(-(spline_T_obs_23$i_spline %*% theta_23))
print_vec("factored_out_direct", factored_out_direct)

# Compute integrand exactly as in C++
integrand_direct <- exp(
    -(spline_grid_T_obs_12$i_spline %*% theta_12) -
    (spline_grid_T_obs_13$i_spline %*% theta_13)
) *
(spline_grid_T_obs_12$m_spline %*% theta_12) *
exp(spline_grid_T_obs_23$i_spline %*% theta_23)

# Compute midpoints exactly as in C++
mid_direct <- (integrand_direct[-1] + integrand_direct[-length(integrand_direct)]) / 2

# Compute cumulative integral exactly as in C++
integral_direct <- c(0, cumsum(mid_direct * dx_grid_T_obs))

# Get indices for T_obs and V_healthy in the grid
# This mimics the C++ .elem() operation
T_obs_indices <- round((T_obs - min(grid_T_obs)) / dx_grid_T_obs) + 1
V_healthy_indices <- round((V_healthy - min(grid_T_obs)) / dx_grid_T_obs) + 1

print_vec("T_obs_indices", T_obs_indices)
print_vec("V_healthy_indices", V_healthy_indices)
print_vec("integral_direct", integral_direct)

# Calculate term_2 using direct indexing as in C++
term_2_direct <- factored_out_direct * (
    integral_direct[T_obs_indices] - integral_direct[V_healthy_indices]
)
print_vec("term_2_direct", term_2_direct)

# Calculate final log-likelihood
log_lik_r_direct <- sum(log(term_1_direct + term_2_direct))

# Print all results for comparison
print(paste("Log-likelihood (C++):", log_lik_cpp))
print(paste("Log-likelihood (R with splines):", log_lik_r))
print(paste("Log-likelihood (R direct match):", log_lik_r_direct))
print(paste("Difference (C++ vs R splines):", abs(log_lik_cpp - log_lik_r)))
print(paste("Difference (C++ vs R direct):", abs(log_lik_cpp - log_lik_r_direct)))
print(paste("Difference (R versions):", abs(log_lik_r - log_lik_r_direct)))
