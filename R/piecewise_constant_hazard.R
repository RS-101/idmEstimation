#### Case-specific model data constructors for PWC ####

# Case 1: Healthy at T_obs, was healthy at V_healthy
setup_case_1_data_pwc <- function(V_0, V_healthy, T_obs, knots_12, knots_13, knots_23) {
  list(
    V_0_values = as.matrix(V_0),
    V_healthy_values = V_healthy,
    V_ill_values = numeric(0),  # Not used in case 1
    T_obs_values = T_obs,
    knots_12 = knots_12,
    knots_13 = knots_13,
    knots_23 = knots_23
  )
}

# Case 2: Died at T_obs, was healthy at V_healthy
setup_case_2_data_pwc <- function(V_0, V_healthy, T_obs, knots_12, knots_13, knots_23) {
  list(
    V_0_values = as.matrix(V_0),
    V_healthy_values = V_healthy,
    V_ill_values = numeric(0),  # Not used in case 2
    T_obs_values = T_obs,
    knots_12 = knots_12,
    knots_13 = knots_13,
    knots_23 = knots_23
  )
}

# Case 3: Ill at T_obs, was healthy at V_healthy, became ill at V_ill
setup_case_3_data_pwc <- function(V_0, V_healthy, V_ill, T_obs, knots_12, knots_13, knots_23) {
  list(
    V_0_values = as.matrix(V_0),
    V_healthy_values = V_healthy,
    V_ill_values = V_ill,
    T_obs_values = T_obs,
    knots_12 = knots_12,
    knots_13 = knots_13,
    knots_23 = knots_23
  )
}

# Case 4: Died at T_obs, was healthy at V_healthy, became ill at V_ill
setup_case_4_data_pwc <- function(V_0, V_healthy, V_ill, T_obs, knots_12, knots_13, knots_23) {
  list(
    V_0_values = as.matrix(V_0),
    V_healthy_values = V_healthy,
    V_ill_values = V_ill,
    T_obs_values = T_obs,
    knots_12 = knots_12,
    knots_13 = knots_13,
    knots_23 = knots_23
  )
}

#### Main setup function ####

setup_cpp_model_pwc <- function(V_0,
                                V_healthy,
                                V_ill,
                                T_obs,
                                status_dead,
                                status_ill,
                                n_knots = 7,
                                knots_12 = NULL,
                                knots_13 = NULL,
                                knots_23 = NULL) {
  
  # Set default knots if not provided
  if (is.null(knots_12)) {
    knots_12 <- seq(min(V_0, na.rm = TRUE),
                    max(V_ill, na.rm = TRUE), length.out = n_knots)
  }
  if (is.null(knots_13)) {
    knots_13 <- seq(min(V_0, na.rm = TRUE),
                    max(T_obs, na.rm = TRUE), length.out = n_knots)
  }
  if (is.null(knots_23)) {
    knots_23 <- seq(min(V_healthy, na.rm = TRUE),
                    max(T_obs, na.rm = TRUE), length.out = n_knots)
  }
  
  # Number of intervals for each transition
  n_lambda_12 <- length(knots_12) - 1
  n_lambda_13 <- length(knots_13) - 1
  n_lambda_23 <- length(knots_23) - 1
  
  # Determine case for each observation
  # Case 1: status_dead = 0 and status_ill = 0 (healthy at T_obs)
  # Case 2: status_dead = 1 and status_ill = 0 (died, was healthy)
  # Case 3: status_dead = 0 and status_ill = 1 (ill at T_obs)
  # Case 4: status_dead = 1 and status_ill = 1 (died, was ill)
  
  case_1_idx <- which(status_dead == 0 & status_ill == 0)
  case_2_idx <- which(status_dead == 1 & status_ill == 0)
  case_3_idx <- which(status_dead == 0 & status_ill == 1)
  case_4_idx <- which(status_dead == 1 & status_ill == 1)
  
  # Create model data for each case
  model_data_list <- list()
  
  if (length(case_1_idx) > 0) {
    model_data_list$case_1 <- setup_case_1_data_pwc(
      V_0[case_1_idx],
      V_healthy[case_1_idx],
      T_obs[case_1_idx],
      knots_12, knots_13, knots_23
    )
  }
  
  if (length(case_2_idx) > 0) {
    model_data_list$case_2 <- setup_case_2_data_pwc(
      V_0[case_2_idx],
      V_healthy[case_2_idx],
      T_obs[case_2_idx],
      knots_12, knots_13, knots_23
    )
  }
  
  if (length(case_3_idx) > 0) {
    model_data_list$case_3 <- setup_case_3_data_pwc(
      V_0[case_3_idx],
      V_healthy[case_3_idx],
      V_ill[case_3_idx],
      T_obs[case_3_idx],
      knots_12, knots_13, knots_23
    )
  }
  
  if (length(case_4_idx) > 0) {
    model_data_list$case_4 <- setup_case_4_data_pwc(
      V_0[case_4_idx],
      V_healthy[case_4_idx],
      V_ill[case_4_idx],
      T_obs[case_4_idx],
      knots_12, knots_13, knots_23
    )
  }
  
  list(
    model_pointer = create_pwc_model_data(model_data_list),
    knots_12 = knots_12,
    knots_13 = knots_13,
    knots_23 = knots_23,
    n_lambda_12 = n_lambda_12,
    n_lambda_13 = n_lambda_13,
    n_lambda_23 = n_lambda_23
  )
}

#### Fit function ####

max_likelihood_pwc <- function(model_config, lambda_0 = NULL) {
  n_lambda_12 <- model_config$n_lambda_12
  n_lambda_13 <- model_config$n_lambda_13
  n_lambda_23 <- model_config$n_lambda_23
  
  n_lambda <- n_lambda_12 + n_lambda_13 + n_lambda_23
  
  obj_fun <- function(long_lambda) {
    res <- calc_log_likelihood_pwc(
      model_config$model_pointer,
      long_lambda[1:n_lambda_12],
      long_lambda[(n_lambda_12 + 1):(n_lambda_12 + n_lambda_13)],
      long_lambda[(n_lambda_12 + n_lambda_13 + 1):n_lambda]
    )
    
    if (!is.finite(res)) return(1e10)
    -res
  }
  
  if (is.null(lambda_0)) {
    lambda_0 <- rep(0.1, n_lambda)
  }
  
  res <- optim(
    par = lambda_0,
    fn = obj_fun,
    method = "L-BFGS-B",
    lower = 1e-8
  )
  
  lambda_hat <- list(
    lambda_12 = res$par[1:n_lambda_12],
    lambda_13 = res$par[(n_lambda_12 + 1):(n_lambda_12 + n_lambda_13)],
    lambda_23 = res$par[(n_lambda_12 + n_lambda_13 + 1):n_lambda]
  )
  
  log_likelihood <- calc_log_likelihood_pwc(
    model_config$model_pointer,
    lambda_hat$lambda_12,
    lambda_hat$lambda_13,
    lambda_hat$lambda_23
  )
  
  list(
    lambda_hat = lambda_hat,
    model_config = model_config,
    log_likelihood = log_likelihood,
    convergence = res$convergence
  )
}

#### Hazard extraction functions ####

create_hazards_pwc <- function(model_config, fit) {
  knots_12 <- model_config$knots_12
  knots_13 <- model_config$knots_13
  knots_23 <- model_config$knots_23
  
  lambda_12 <- fit$lambda_hat$lambda_12
  lambda_13 <- fit$lambda_hat$lambda_13
  lambda_23 <- fit$lambda_hat$lambda_23
  
  # Piecewise constant hazard function
  # Returns the hazard rate at time x
  make_pwc_hazard <- function(x, knots, lambda) {
    n <- length(x)
    result <- numeric(n)
    m <- length(lambda)
    
    for (i in seq_len(n)) {
      # Find which interval x[i] falls into
      idx <- findInterval(x[i], knots, left.open = FALSE)
      # Ensure idx is within valid range [1, m]
      idx <- max(1, min(idx, m))
      result[i] <- lambda[idx]
    }
    result
  }
  
  # Cumulative hazard function
  # Returns the cumulative hazard up to time x
  make_cumulative_hazard <- function(x, knots, lambda) {
    n <- length(x)
    result <- numeric(n)
    m <- length(lambda)
    
    for (i in seq_len(n)) {
      total <- 0
      for (j in seq_len(m)) {
        left <- knots[j]
        right <- knots[j + 1]
        len <- max(0, min(x[i], right) - left)
        total <- total + lambda[j] * len
      }
      result[i] <- total
    }
    result
  }
  
  hazards <- list(
    a12 = function(x) make_pwc_hazard(x, knots_12, lambda_12),
    a13 = function(x) make_pwc_hazard(x, knots_13, lambda_13),
    a23 = function(x) make_pwc_hazard(x, knots_23, lambda_23),
    A12 = function(x) make_cumulative_hazard(x, knots_12, lambda_12),
    A13 = function(x) make_cumulative_hazard(x, knots_13, lambda_13),
    A23 = function(x) make_cumulative_hazard(x, knots_23, lambda_23)
  )
  
  # add class to hazards
  class(hazards) <- c("idm_hazards", class(hazards))
  hazards
}

#### Main fitting function ####

fit_idm_pwc <- function(data,
                        knots_12 = NULL,
                        knots_13 = NULL,
                        knots_23 = NULL,
                        n_knots = 7,
                        lambda_0 = NULL,
                        verbose = TRUE) {
  
  # verify data components
  stopifnot(all(c("V_0", "V_healthy", "V_ill", "T_obs", "status_dead", "status_ill") %in% names(data)))
  
  V_0 <- data$V_0
  V_healthy <- data$V_healthy
  V_ill <- data$V_ill
  T_obs <- data$T_obs
  status_dead <- data$status_dead
  status_ill <- data$status_ill
  
  # Setup model configuration
  model_config <- setup_cpp_model_pwc(
    V_0 = V_0, V_healthy = V_healthy, V_ill = V_ill, T_obs = T_obs,
    status_dead = status_dead, status_ill = status_ill, n_knots = n_knots,
    knots_12 = knots_12, knots_13 = knots_13, knots_23 = knots_23
  )
  
  if (verbose) {
    message(sprintf("Fitting PWC model with %d knots per transition", n_knots))
    message(sprintf("  n_lambda_12 = %d", model_config$n_lambda_12))
    message(sprintf("  n_lambda_13 = %d", model_config$n_lambda_13))
    message(sprintf("  n_lambda_23 = %d", model_config$n_lambda_23))
  }
  
  # Fit the model
  fit <- max_likelihood_pwc(model_config, lambda_0)
  
  if (verbose) {
    message(sprintf("Optimization converged: %s", ifelse(fit$convergence == 0, "YES", "NO")))
    message(sprintf("Log-likelihood: %.4f", fit$log_likelihood))
  }
  
  # Create hazard functions
  hazards <- create_hazards_pwc(model_config, fit)
  
  list(
    hazards = hazards,
    fit = fit,
    model_config = model_config
  )
}
