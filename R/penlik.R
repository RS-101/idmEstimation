library(splines2)
#### Create spline hazard ####
make_spline_mat <- function(x, knots, degree, use_bspline) {

  if(!use_bspline){

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
  } else {
    n_knots <- length(knots)
    i_spline_mat <- splines2::bSpline(x,
                                      degree = degree,
                                      knots = knots[-c(1, n_knots)],
                                      Boundary.knots = knots[c(1, n_knots)],
                                      intercept = TRUE,
                                      warn.outside = FALSE,
                                      integral = T
    )

    m_spline_mat <- splines2::bSpline(x,
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

  }

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
setup_case_1_data <- function(V_0, V_healthy, T_obs, knots_12, knots_13, knots_23, degree, use_bspline) {
  # Create splines at observation points
  spline_V_0_12 <- make_spline_mat(V_0, knots_12, degree, use_bspline)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13, degree, use_bspline)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree, use_bspline)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree, use_bspline)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree, use_bspline)

  # Create integration grid from V_healthy to T_obs
  grid_T_obs <- seq(min(V_healthy), max(T_obs), length.out = 250)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12, degree, use_bspline)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13, degree, use_bspline)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23, degree, use_bspline)

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
setup_case_2_data <- function(V_0, V_healthy, T_obs, knots_12, knots_13, knots_23, degree, use_bspline) {
  # Create splines at observation points
  spline_V_0_12 <- make_spline_mat(V_0, knots_12, degree, use_bspline)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13, degree, use_bspline)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree, use_bspline)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree, use_bspline)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree, use_bspline)

  # Create integration grid from V_healthy to T_obs
  grid_T_obs <- seq(min(V_healthy), max(T_obs), length.out = 250)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12, degree, use_bspline)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13, degree, use_bspline)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23, degree, use_bspline)

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
setup_case_3_data <- function(V_0, V_healthy, V_ill, T_obs, knots_12, knots_13, knots_23, degree, use_bspline) {
  # Create splines at observation points
  spline_V_0_12 <- make_spline_mat(V_0, knots_12, degree, use_bspline)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13, degree, use_bspline)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree, use_bspline)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree, use_bspline)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree, use_bspline)

  # Create integration grid from V_healthy to V_ill
  grid_V_ill <- seq(min(V_healthy), max(V_ill), length.out = 250)
  dx_grid_V_ill <- diff(grid_V_ill)[1]

  spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots_12, degree, use_bspline)
  spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots_13, degree, use_bspline)
  spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots_23, degree, use_bspline)

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
setup_case_4_data <- function(V_0, V_healthy, V_ill, T_obs, knots_12, knots_13, knots_23, degree, use_bspline) {
  # Create splines at observation points
  spline_V_0_12 <- make_spline_mat(V_0, knots_12, degree, use_bspline)
  spline_V_0_13 <- make_spline_mat(V_0, knots_13, degree, use_bspline)

  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree, use_bspline)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree, use_bspline)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree, use_bspline)

  # Create integration grid from V_healthy to V_ill
  grid_V_ill <- seq(min(V_healthy), max(V_ill), length.out = 250)
  dx_grid_V_ill <- diff(grid_V_ill)[1]

  spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots_12, degree, use_bspline)
  spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots_13, degree, use_bspline)
  spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots_23, degree, use_bspline)

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
                            n_knots,
                            knots_12 = NULL,
                            knots_13 = NULL,
                            knots_23 = NULL,
                            degree,
                            use_bspline) {

  # Set default knots if not provided
  if (is.null(knots_12)) {
    knots_12 <- seq(min(V_0, na.rm = T),
                    max(V_ill, na.rm = T), length.out = n_knots)
  }
  if (is.null(knots_13)) {
    knots_13 <- seq(min(V_0, na.rm = T),
                    max(T_obs, na.rm = T), length.out = n_knots)
  }
  if (is.null(knots_23)) {
    knots_23 <- seq(min(V_healthy, na.rm = T),
                    max(T_obs, na.rm = T), length.out = n_knots)
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
      knots_12, knots_13, knots_23, degree, use_bspline
    )
  }

  if (length(case_2_idx) > 0) {
    model_data_list$case_2 <- setup_case_2_data(
      V_0[case_2_idx],
      V_healthy[case_2_idx],
      T_obs[case_2_idx],
      knots_12, knots_13, knots_23, degree, use_bspline
    )
  }

  if (length(case_3_idx) > 0) {
    model_data_list$case_3 <- setup_case_3_data(
      V_0[case_3_idx],
      V_healthy[case_3_idx],
      V_ill[case_3_idx],
      T_obs[case_3_idx],
      knots_12, knots_13, knots_23, degree, use_bspline
    )
  }

  if (length(case_4_idx) > 0) {
    model_data_list$case_4 <- setup_case_4_data(
      V_0[case_4_idx],
      V_healthy[case_4_idx],
      V_ill[case_4_idx],
      T_obs[case_4_idx],
      knots_12, knots_13, knots_23, degree, use_bspline
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
        n_theta_23 = n_theta_23,
        use_bspline = use_bspline)
}


#### Fit ####
max_pen_likelihood <- function(model_config, kappa_12, kappa_13, kappa_23, long_theta_0 = NULL) {
  n_theta_12 <- model_config$n_theta_12
  n_theta_13 <- model_config$n_theta_13
  n_theta_23 <- model_config$n_theta_23

  n_theta <- n_theta_12 + n_theta_13 + n_theta_23

  obj_fun <- function(long_theta) {

    res <- calc_penalized_log_likelihood(
      model_config$model_pointer,
      long_theta[1:n_theta_12],
      long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
      long_theta[(n_theta_12 + n_theta_13 + 1):n_theta],
      kappa_12,
      kappa_13,
      kappa_23
    )$penalized_log_likelihood

    if (!is.finite(res)) return(1e10)
    -res
  }


  if(is.null(long_theta_0)) {
    long_theta_0 <- rep(0.5, n_theta)
  }

  res <- optim(
    par = long_theta_0,
    fn = obj_fun,
    method= "L-BFGS-B",
    lower = 0)

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
    log_likelihood = max(pl_at_theta_hat$log_likelihood, -1e10),
    penalized_log_likelihood = max(pl_at_theta_hat$penalized_log_likelihood, -1e10),
    penalty = pl_at_theta_hat$penalty,
    convergence = res$convergence
  )
}


max_pen_likelihood_unconstrained <- function(model_config, kappa_12, kappa_13, kappa_23, long_theta_0 = NULL) {
  n_theta_12 <- model_config$n_theta_12
  n_theta_13 <- model_config$n_theta_13
  n_theta_23 <- model_config$n_theta_23

  n_theta <- n_theta_12 + n_theta_13 + n_theta_23

  obj_fun <- function(long_theta) {
    long_theta <- exp(long_theta)

    -calc_penalized_log_likelihood(
      model_config$model_pointer,
      long_theta[1:n_theta_12],
      long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
      long_theta[(n_theta_12 + n_theta_13 + 1):n_theta],
      kappa_12,
      kappa_13,
      kappa_23
    )$penalized_log_likelihood
  }


  if(is.null(long_theta_0)) {
    long_theta_0 <- rep(-1, n_theta)
  }

  res <- optim(
    par = long_theta_0,
    fn = obj_fun,
    method= "BFGS"
  )

  res$par <- exp(res$par)

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





approx_cv <- function(pl_optim) {

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

  theta_hat <- pl_optim$theta_hat
  ll_value <- pl_optim$log_likelihood
  pl_value <- pl_optim$penalized_log_likelihood

  if(abs(ll_value) > 1e9 | abs(pl_value) > 1e9){
    return(-1e9)
  }

  model_pointer <- pl_optim$model_config$model_pointer
  n_theta_12 <- pl_optim$model_config$n_theta_12
  n_theta_13 <- pl_optim$model_config$n_theta_13
  n_theta_23 <- pl_optim$model_config$n_theta_23
  n_theta <- n_theta_12 + n_theta_13 + n_theta_23


  ll_in_long_theta <- function(long_theta){
    calc_log_likelihood(md_ptr = model_pointer,
                        theta_12 = long_theta[1:n_theta_12],
                        theta_13 = long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
                        theta_23 = long_theta[(n_theta_12 + n_theta_13 + 1):n_theta])
  }

  pl_in_long_theta <- function(long_theta){
    calc_penalized_log_likelihood(md_ptr = model_pointer,
                        theta_12 = long_theta[1:n_theta_12],
                        theta_13 = long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
                        theta_23 = long_theta[(n_theta_12 + n_theta_13 + 1):n_theta],
                        kappa_12 = pl_optim$kappa_12,
                        kappa_13 = pl_optim$kappa_13,
                        kappa_23 = pl_optim$kappa_23)$penalized_log_likelihood
  }


  long_theta_hat <- unlist(theta_hat)

  H_pl <- numDeriv::hessian(pl_in_long_theta, long_theta_hat)
  H_ll  <- numDeriv::hessian(ll_in_long_theta, long_theta_hat)

  res <- safe_solve(H_pl, H_ll) # solves H_pl X = H_ll
  tr_val <- sum(diag(res$X))

  approx_cv <- ll_value - tr_val
  approx_cv
}

create_hazards <- function(model_config, fit) {
  hazards <- list(a12 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_12, model_config$degree,
                                  model_config$use_bspline)$m_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_12)
  },
  a13 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_13, model_config$degree,
                                  model_config$use_bspline)$m_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_13)
  },
  a23 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_23, model_config$degree,
                                  model_config$use_bspline)$m_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_23)
  },
  A12 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_12, model_config$degree,
                                  model_config$use_bspline)$i_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_12)
  },
  A13 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_13, model_config$degree,
                                  model_config$use_bspline)$i_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_13)
  },
  A23 = function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_23, model_config$degree,
                                  model_config$use_bspline)$i_spline
    as.vector(spline_mat %*% fit$theta_hat$theta_23)
  })

  # add class to hazards
  class(hazards) <- c("idm_hazards", class(hazards))
  hazards
}

fit_idm <- function(data,
                    run_in_parallel = TRUE,
                    knots_12 = NULL,
                    knots_13 = NULL,
                    knots_23 = NULL,
                    degree = 3,
                    n_knots = 7,
                    kappa_values = NULL,
                    use_bspline = F,
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
    knots_12 = knots_12, knots_13 = knots_13, knots_23 = knots_23, degree = degree,
    use_bspline = use_bspline
  )


  if (is.null(kappa_values)) {
    kappa_grid <- expand.grid(
      kappa_12 = 10^seq(-2, 20, length.out = 9),
      kappa_13 = 10^seq(-2, 20, length.out = 9),
      kappa_23 = 10^seq(-2, 20, length.out = 9)
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

  if(run_in_parallel){
    res <- parallel::mclapply(seq_len(n), function(i) {
      kappa_12 <- kappa_grid$kappa_12[i]
      kappa_13 <- kappa_grid$kappa_13[i]
      kappa_23 <- kappa_grid$kappa_23[i]

      fit <- max_pen_likelihood(model_config, kappa_12, kappa_13, kappa_23)
      approx_cv_value <- approx_cv(fit)

      list(fit = fit, cv = approx_cv_value)
    }, mc.cores = max(1, parallel::detectCores() - 1))

    fits <- lapply(res, `[[`, "fit")
    cvs  <- vapply(res, function(x) x$cv, numeric(1))


  } else {
    for (i in seq_len(n)) {
      kappa_12 <- kappa_grid$kappa_12[i]
      kappa_13 <- kappa_grid$kappa_13[i]
      kappa_23 <- kappa_grid$kappa_23[i]

      fit <- max_pen_likelihood(model_config, kappa_12, kappa_13, kappa_23)
      approx_cv_value <- approx_cv(fit)

      fits[[i]] <- fit
      cvs[i] <- approx_cv_value

      if (verbose) {
        message(sprintf(
          "[%d/%d] kappa=(%.2g, %.2g, %.2g)  cv=%.4f",
          i, n, kappa_12, kappa_13, kappa_23, cvs[i]
        ))
      }
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
        cv_values = cvs,
        other_fits = fits,
        model_config = model_config
    ))
}


