library(splines2)
#### Create spline hazard ####
make_spline_mat <- function(x, knots, degree = 3) {
  n_knots <- length(knots)
  i_spline_mat <- splines2::iSpline(x,
    degree = degree,
    knots = knots[-c(1, n_knots)],
    Boundary.knots = knots[c(1, n_knots)],
    intercept = TRUE,
    warn.outside = FALSE
  )

  m_spline_mat <- splines2::mSpline(x,
    degree = degree,
    knots = knots[-c(1, n_knots)],
    Boundary.knots = knots[c(1, n_knots)],
    intercept = TRUE,
    warn.outside = FALSE
  )

  spline_mat_list <- list(
    i_spline = i_spline_mat,
    m_spline = m_spline_mat
  )

  class(spline_mat_list) <- c("spline_mat_list", class(spline_mat_list))
  spline_mat_list
}

calculate_penalty_matrix <- function(knots, degree = 3) {
  n_knots <- length(knots)
  d <- diff(knots)
  g_ab <- splines2::mSpline(x = knots,
                  knots = knots[-c(1, n_knots)],
                  Boundary.knots = range(knots),
                  degree = degree,
                  derivs = 2,
                  intercept = TRUE
  )
  knots_mid <- knots[-length(knots)] + d / 2
  g_ab_mid <- splines2::mSpline(x = knots_mid,
                      knots = knots[-c(1, n_knots)],
                      Boundary.knots = range(knots),
                      degree = degree,
                      derivs = 2,
                      intercept = TRUE
  )
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) +
      4 * crossprod(d * g_ab_mid, g_ab_mid) +
      crossprod(d * g_b, g_b)) / 6
}


#### Case-specific model data constructors ####

# Case 1: Healthy at T_obs, was healthy at V_healthy
# Integration from V_healthy to T_obs
setup_case_1_data <- function(V_0, V_healthy, T_obs, knots_12, knots_13, knots_23, degree = 3) {
  # Create splines at observation points
  spline_V_0_12 <- make_spline_mat(V_0, knots_12, degree)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13, degree)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree)

  # Create integration grid from V_healthy to T_obs
  grid_T_obs <- seq(min(V_healthy), max(T_obs), length.out = 250)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12, degree)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13, degree)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23, degree)

  list(
    V_0_i_spline_mat_12 = spline_V_0_12$i_spline,
    V_0_i_spline_mat_13 = spline_V_0_13$i_spline,

    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,

    grid_T_obs_i_spline_mat_12 = spline_grid_T_obs_12$i_spline,
    grid_T_obs_i_spline_mat_13 = spline_grid_T_obs_13$i_spline,
    grid_T_obs_i_spline_mat_23 = spline_grid_T_obs_23$i_spline,
    grid_T_obs_m_spline_mat_12 = spline_grid_T_obs_12$m_spline,

    dx_grid_T_obs = dx_grid_T_obs,

    T_obs_values = T_obs,
    V_healthy_values = V_healthy
  )
}

# Case 2: Died at T_obs, was healthy at V_healthy
# Integration from V_healthy to T_obs, needs hazard rates at T_obs
setup_case_2_data <- function(V_0, V_healthy, T_obs, knots_12, knots_13, knots_23, degree = 3) {
  # Create splines at observation points
  spline_V_0_12 <- make_spline_mat(V_0, knots_12, degree)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13, degree)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree)

  # Create integration grid from V_healthy to T_obs
  grid_T_obs <- seq(min(V_healthy), max(T_obs), length.out = 250)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12, degree)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13, degree)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23, degree)

  list(
    V_0_i_spline_mat_12 = spline_V_0_12$i_spline,
    V_0_i_spline_mat_13 = spline_V_0_13$i_spline,

    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,
    T_obs_m_spline_mat_13 = spline_T_obs_13$m_spline,
    T_obs_m_spline_mat_23 = spline_T_obs_23$m_spline,

    grid_T_obs_i_spline_mat_12 = spline_grid_T_obs_12$i_spline,
    grid_T_obs_i_spline_mat_13 = spline_grid_T_obs_13$i_spline,
    grid_T_obs_i_spline_mat_23 = spline_grid_T_obs_23$i_spline,
    grid_T_obs_m_spline_mat_12 = spline_grid_T_obs_12$m_spline,

    dx_grid_T_obs = dx_grid_T_obs,

    T_obs_values = T_obs,
    V_healthy_values = V_healthy
  )
}

# Case 3: Ill at T_obs, was healthy at V_healthy, became ill at V_ill
# Integration from V_healthy to V_ill
setup_case_3_data <- function(V_0, V_healthy, V_ill, T_obs, knots_12, knots_13, knots_23, degree = 3) {
  # Create splines at observation points
  spline_V_0_12 <- make_spline_mat(V_0, knots_12, degree)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13, degree)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree)

  # Create integration grid from V_healthy to V_ill
  grid_V_ill <- seq(min(V_healthy), max(V_ill), length.out = 250)
  dx_grid_V_ill <- diff(grid_V_ill)[1]

  spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots_12, degree)
  spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots_13, degree)
  spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots_23, degree)

  list(
    V_0_i_spline_mat_12 = spline_V_0_12$i_spline,
    V_0_i_spline_mat_13 = spline_V_0_13$i_spline,

    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,

    grid_V_ill_i_spline_mat_12 = spline_grid_V_ill_12$i_spline,
    grid_V_ill_i_spline_mat_13 = spline_grid_V_ill_13$i_spline,
    grid_V_ill_i_spline_mat_23 = spline_grid_V_ill_23$i_spline,
    grid_V_ill_m_spline_mat_12 = spline_grid_V_ill_12$m_spline,

    dx_grid_V_ill = dx_grid_V_ill,

    T_obs_values = T_obs,
    V_healthy_values = V_healthy,
    V_ill_values = V_ill
  )
}

# Case 4: Died at T_obs, was healthy at V_healthy, became ill at V_ill
# Integration from V_healthy to V_ill, needs hazard rate at T_obs
setup_case_4_data <- function(V_0, V_healthy, V_ill, T_obs, knots_12, knots_13, knots_23, degree = 3) {
  # Create splines at observation points
  spline_V_0_12 <- make_spline_mat(V_0, knots_12, degree)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13, degree)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree)

  # Create integration grid from V_healthy to V_ill
  grid_V_ill <- seq(min(V_healthy), max(V_ill), length.out = 250)
  dx_grid_V_ill <- diff(grid_V_ill)[1]

  spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots_12, degree)
  spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots_13, degree)
  spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots_23, degree)

  list(
    V_0_i_spline_mat_12 = spline_V_0_12$i_spline,
    V_0_i_spline_mat_13 = spline_V_0_13$i_spline,

    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,
    T_obs_m_spline_mat_23 = spline_T_obs_23$m_spline,

    grid_V_ill_i_spline_mat_12 = spline_grid_V_ill_12$i_spline,
    grid_V_ill_i_spline_mat_13 = spline_grid_V_ill_13$i_spline,
    grid_V_ill_i_spline_mat_23 = spline_grid_V_ill_23$i_spline,
    grid_V_ill_m_spline_mat_12 = spline_grid_V_ill_12$m_spline,

    dx_grid_V_ill = dx_grid_V_ill,

    T_obs_values = T_obs,
    V_healthy_values = V_healthy,
    V_ill_values = V_ill
  )
}

#### Main setup function ####

setup_cpp_model <- function(V_0,
                            V_healthy,
                            V_ill,
                            T_obs,
                            status_dead,
                            status_ill,
                            n_knots = 7,
                            knots_12 = NULL,
                            knots_13 = NULL,
                            knots_23 = NULL,
                            degree = 3) {

  # Set default knots if not provided
  if (is.null(knots_12)) {
    knots_12 <- seq(min(V_0), max(T_obs), length.out = n_knots)
  }
  if (is.null(knots_13)) {
    knots_13 <- seq(min(V_0), max(T_obs), length.out = n_knots)
  }
  if (is.null(knots_23)) {
    knots_23 <- seq(min(V_0), max(T_obs), length.out = n_knots)
  }

  n_theta_12 <- length(knots_12) + degree - 1
  n_theta_13 <- length(knots_13) + degree - 1
  n_theta_23 <- length(knots_23) + degree - 1



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
    model_data_list$case_1 <- setup_case_1_data(
      V_0[case_1_idx],
      V_healthy[case_1_idx],
      T_obs[case_1_idx],
      knots_12, knots_13, knots_23, degree
    )
  }

  if (length(case_2_idx) > 0) {
    model_data_list$case_2 <- setup_case_2_data(
      V_0[case_2_idx],
      V_healthy[case_2_idx],
      T_obs[case_2_idx],
      knots_12, knots_13, knots_23, degree
    )
  }

  if (length(case_3_idx) > 0) {
    model_data_list$case_3 <- setup_case_3_data(
      V_0[case_3_idx],
      V_healthy[case_3_idx],
      V_ill[case_3_idx],
      T_obs[case_3_idx],
      knots_12, knots_13, knots_23, degree
    )
  }

  if (length(case_4_idx) > 0) {
    model_data_list$case_4 <- setup_case_4_data(
      V_0[case_4_idx],
      V_healthy[case_4_idx],
      V_ill[case_4_idx],
      T_obs[case_4_idx],
      knots_12, knots_13, knots_23, degree
    )
  }


  model_data_list$penalty_matrix_12 <- calculate_penalty_matrix(knots_12, degree)
  model_data_list$penalty_matrix_13 <- calculate_penalty_matrix(knots_13, degree)
  model_data_list$penalty_matrix_23 <- calculate_penalty_matrix(knots_23, degree)

  list(model_pointer = create_penlik_model_data(model_data_list),
        knots_12 = knots_12,
        knots_13 = knots_13,
        knots_23 = knots_23,
        degree = degree,
        n_theta_12 = n_theta_12,
        n_theta_13 = n_theta_13,
        n_theta_23 = n_theta_23)
}


#### Fit ####
max_pen_likelihood <- function(model_config, kappa_12, kappa_13, kappa_23) {
  n_theta_12 <- model_config$n_theta_12
  n_theta_13 <- model_config$n_theta_13
  n_theta_23 <- model_config$n_theta_23

  n_theta <- n_theta_12 + n_theta_13 + n_theta_23

  obj_fun <- function(long_theta) {
    calc_penalized_log_likelihood(
      model_config$model_pointer,
      long_theta[1:n_theta_12],
      long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
      long_theta[(n_theta_12 + n_theta_13 + 1):n_theta],
      kappa_12,
      kappa_13,
      kappa_23
    )$penalized_log_likelihood
  }

  res <- optim(
    par = rep(0.1, n_theta),
    fn = obj_fun,
    control = list(fnscale = -1, maxit = 500)
  )

  theta_hat <- list(
    theta_12 = res$par[1:n_theta_12],
    theta_13 = res$par[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
    theta_23 = res$par[(n_theta_12 + n_theta_13 + 1):n_theta]
  )

  pl_at_theta_hat <- calc_penalized_log_likelihood(
    model_config$model_pointer,
    theta_hat$theta_12,
    theta_hat$theta_13,
    theta_hat$theta_23,
    kappa_12,
    kappa_13,
    kappa_23
  )

  list(
    theta_hat = theta_hat,
    model_config = model_config,
    kappa_12 = kappa_12,
    kappa_13 = kappa_13,
    kappa_23 = kappa_23,
    log_likelihood = pl_at_theta_hat$log_likelihood,
    penalized_log_likelihood = pl_at_theta_hat$penalized_log_likelihood,
    penalty = pl_at_theta_hat$penalty,
    convergence = res$convergence
  )
}

safe_solve <- function(H_pl, H_ll, ridge_start = 0, tol = 1e-8, max_iter = 8) {
  # Ensure symmetry (Hessians should be symmetric, but numerics can drift)
  H <- (H_pl + t(H_pl)) / 2

  # Scale for the ridge so lambda is magnitude-aware
  sc <- mean(abs(diag(H)))
  if (!is.finite(sc) || sc == 0) sc <- 1

  # Try Cholesky with increasing ridge
  lambda <- ridge_start
  for (i in 0:max_iter) {
    H_reg <- H + (lambda * sc + .Machine$double.eps) * diag(nrow(H))
    R <- try(chol(H_reg), silent = TRUE)  # H_reg = t(R) %*% R
    if (!inherits(R, "try-error")) {
      # Solve H_reg X = H_ll via two triangular solves (no explicit inverse)
      Y <- forwardsolve(t(R), H_ll)   # t(R) Y = H_ll
      X <- backsolve(R, Y)            # R X = Y
      return(list(X = X, method = "chol", lambda = lambda * sc, iters = i))
    }
    lambda <- if (lambda == 0) tol else lambda * 10
  }

  # Fallback: SVD pseudo-inverse with thresholding
  sv <- svd(H)
  thr <- max(tol, max(sv$d) * 1e-12)
  d_inv <- ifelse(sv$d > thr, 1 / sv$d, 0)
  X <- sv$v %*% (d_inv * t(sv$u) %*% H_ll)
  list(X = X, method = "svd", lambda = lambda * sc, rank = sum(sv$d > thr))
}


approx_cv <- function(fit_result) {
  model_config <- fit_result$model_config
  theta_hat <- fit_result$theta_hat
  kappa_12 <- fit_result$kappa_12
  kappa_13 <- fit_result$kappa_13
  kappa_23 <- fit_result$kappa_23

  n_theta_12 <- model_config$n_theta_12
  n_theta_13 <- model_config$n_theta_13
  n_theta_23 <- model_config$n_theta_23
  n_theta <- n_theta_12 + n_theta_13 + n_theta_23

  long_theta_hat <- c(theta_hat$theta_12, theta_hat$theta_13, theta_hat$theta_23)

  calc_ll_in_long_theta <- function(long_theta) {
    calc_log_likelihood(
      model_config$model_pointer,
      long_theta[1:n_theta_12],
      long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
      long_theta[(n_theta_12 + n_theta_13 + 1):n_theta]
    )
  }

  calc_pl_in_long_theta <- function(long_theta) {
    calc_penalized_log_likelihood(
      model_config$model_pointer,
      long_theta[1:n_theta_12],
      long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
      long_theta[(n_theta_12 + n_theta_13 + 1):n_theta],
      kappa_12,
      kappa_13,
      kappa_23
    )$penalized_log_likelihood
  }

  H_pl <- numDeriv::hessian(calc_pl_in_long_theta, long_theta_hat)
  H_ll <- numDeriv::hessian(calc_ll_in_long_theta, long_theta_hat)

  res <- safe_solve(H_pl, H_ll)
  tr_val <- sum(diag(res$X))

  approx_cv_value <- fit_result$log_likelihood - tr_val

  fit_result$approx_cv_value <- approx_cv_value
  fit_result$trace_value <- tr_val
  fit_result
}

create_hazards <- function(model_config, fit) {
  hazards <- list(a12 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_12, model_config$degree)$m_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_12)
  },
  a13 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_13, model_config$degree)$m_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_13)
  },
  a23 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_23, model_config$degree)$m_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_23)
  },
  A12 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_12, model_config$degree)$i_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_12)
  },
  A13 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_13, model_config$degree)$i_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_13)
  },
  A23 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_23, model_config$degree)$i_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_23)
  })

  # add class to hazards
  class(hazards) <- c("idm_hazards", class(hazards))
  hazards
}

fit_idm <- function(data,
                    knots_12 = NULL,
                    knots_13 = NULL,
                    knots_23 = NULL,
                    degree = 3,
                    n_knots = 7,
                    kappa_values = NULL,
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
  model_config <- setup_cpp_model(
   V_0 = V_0, V_healthy = V_healthy, V_ill = V_ill, T_obs = T_obs,
    status_dead = status_dead, status_ill = status_ill, n_knots = n_knots,
    knots_12 = knots_12, knots_13 = knots_13, knots_23 = knots_23, degree = degree
  )


  if (is.null(kappa_values)) {
    kappa_grid <- expand.grid(
      kappa_12 = 10^seq(-3, 3, length.out = 7),
      kappa_13 = 10^seq(-3, 3, length.out = 7),
      kappa_23 = 10^seq(-3, 3, length.out = 7)
    )
  } else {
    kappa_grid <- as.data.frame(kappa_values)
    if (!all(c("kappa_12", "kappa_13", "kappa_23") %in% colnames(kappa_grid))) {
      stop("kappa_values must contain columns 'kappa_12', 'kappa_13', and 'kappa_23'")
    }
  }
  n <- nrow(kappa_grid)

  cvs <- numeric(n)
  fits <- vector("list", n)

  for (i in seq_len(n)) {
    kappa_12 <- kappa_grid$kappa_12[i]
    kappa_13 <- kappa_grid$kappa_13[i]
    kappa_23 <- kappa_grid$kappa_23[i]

    fit <- max_pen_likelihood(model_config, kappa_12, kappa_13, kappa_23)
    fit_cv <- approx_cv(fit)

    fits[[i]] <- fit_cv
    cvs[i] <- fit_cv$approx_cv_value

    if (verbose) {
      message(sprintf(
        "[%d/%d] kappa=(%.2g, %.2g, %.2g)  cv=%.4f",
        i, n, kappa_12, kappa_13, kappa_23, cvs[i]
      ))
    }
  }

  best <- which.max(cvs)
  final_kappas <- c(
    kappa_12 = kappa_grid$kappa_12[best],
    kappa_13 = kappa_grid$kappa_13[best],
    kappa_23 = kappa_grid$kappa_23[best]
  )
  final_fit <- fits[[best]]
  final_cv <- final_fit
  hazards <- create_hazards(model_config, final_fit)

  return(list(
        hazards = hazards,
        kappas = final_kappas,
        fit = final_fit,
        cv_value = final_cv$approx_cv_value,
        model_config = model_config
    ))
}





fit_idm_greedy <- function(data,
                    knots_12 = NULL,
                    knots_13 = NULL,
                    knots_23 = NULL,
                    initial_kappas = c(1,1,1),
                    max_iter = 2,
                    degree = 3,
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
  model_config <- setup_cpp_model(
    V_0, V_healthy, V_ill, T_obs,
    status_dead, status_ill,
    knots_12, knots_13, knots_23, degree
  )


    # Initialize log of kappas
    b2 <- log(initial_kappas)

    # For each of the 3 transitions (0→1, 0→2, 1→2)
    for(iter in 1:max_iter) {
        for(ii in 1:3) {
            if(verbose) cat(sprintf("\nIteration %d, Parameter %d\n", iter, ii))

            # Initial step size
            uh <- 1.0

            # Initialize search points
            if(b2[ii] > 2.0) {
                u <- c(
                    b2[ii] - uh,  # u[1]
                    b2[ii],       # u[2]
                    b2[ii] + uh   # u[3]
                )
            } else {
                u <- c(
                    3.0 - uh,  # u[1]
                    3.0,       # u[2]
                    3.0 + uh   # u[3]
                )
            }

            # Evaluate initial points
            fu <- numeric(3)

            # First point
            b2[ii] <- u[1]
            kappas <- exp(b2)
            fit <- max_pen_likelihood(model_config, kappas[1], kappas[2], kappas[3])
            fu[1] <- approx_cv(fit)$approx_cv_value

            # Second point
            b2[ii] <- u[2]
            kappas <- exp(b2)
            fit <- max_pen_likelihood(model_config, kappas[1], kappas[2], kappas[3])
            fu[2] <- approx_cv(fit)$approx_cv_value

            # Third point
            b2[ii] <- u[3]
            kappas <- exp(b2)
            fit <- max_pen_likelihood(model_config, kappas[1], kappas[2], kappas[3])
            fu[3] <- approx_cv(fit)$approx_cv_value

            # Iterate through different step sizes
            step_sizes <- c(1.0, 0.5, 0.25, 0.125)
            for(uh in step_sizes) {
                if(verbose) cat(sprintf("Step size: %.3f\n", uh))

                # Case 1: Moving in positive direction
                if(fu[1] < fu[2] && fu[2] < fu[3]) {
                    if(verbose) cat("Case 1: Moving positive\n")
                    for(j in 1:5) {
                        u[1] <- u[2]
                        fu[1] <- fu[2]
                        u[2] <- u[3]
                        fu[2] <- fu[3]
                        u[3] <- u[2] + uh

                        b2[ii] <- u[3]
                        kappas <- exp(b2)
                        fit <- max_pen_likelihood(model_config, kappas[1], kappas[2], kappas[3])
                        fu[3] <- approx_cv(fit)$approx_cv_value

                        if(verbose) {
                            cat(sprintf("u: %.3f, %.3f, %.3f\n", u[1], u[2], u[3]))
                            cat(sprintf("fu: %.3f, %.3f, %.3f\n", fu[1], fu[2], fu[3]))
                        }
                    }
                }

                # Case 2: Zoom between points
                else if(fu[1] < fu[2] && fu[2] > fu[3]) {
                    if(verbose) cat("Case 2: Zooming\n")
                    u[3] <- u[2]
                    fu[3] <- fu[2]
                    u[2] <- u[1] + uh

                    b2[ii] <- u[2]
                    kappas <- exp(b2)
                    fit <- max_pen_likelihood(model_config, kappas[1], kappas[2], kappas[3])
                    fu[2] <- approx_cv(fit)$approx_cv_value
                }

                # Case 3: Moving in negative direction
                else if(fu[1] > fu[2] && fu[2] > fu[3]) {
                    if(verbose) cat("Case 3: Moving negative\n")
                    u[3] <- u[2]
                    fu[3] <- fu[2]
                    u[2] <- u[1]
                    fu[2] <- fu[1]
                    u[1] <- u[3] - 2.0 * uh

                    b2[ii] <- u[1]
                    kappas <- exp(b2)
                    fit <- max_pen_likelihood(model_config, kappas[1], kappas[2], kappas[3])
                    fu[1] <- approx_cv(fit)$approx_cv_value
                }
            }

            # Update b2 with best value found
            if(fu[1] < fu[2] && fu[2] < fu[3]) {
                b2[ii] <- u[3]
            } else if(fu[1] < fu[2] && fu[2] > fu[3]) {
                b2[ii] <- u[2]
            } else if(fu[1] > fu[2] && fu[2] > fu[3]) {
                b2[ii] <- u[1]
            }

            if(verbose) cat(sprintf("Final value for parameter %d: %.3f\n", ii, b2[ii]))
        }
    }

    # Convert back to kappa values
    final_kappas <- exp(b2)

    # Final fit with optimal kappas
    final_fit <- max_pen_likelihood(model_config, final_kappas[1], final_kappas[2], final_kappas[3])
    final_cv <- approx_cv(final_fit)


    hazards <- create_hazards(model_config, final_fit)

    return(list(
        hazards = hazards,
        kappas = final_kappas,
        fit = final_fit,
        cv_value = final_cv$approx_cv_value,
        log_kappas = b2,
        model_config = model_config
    ))
}


