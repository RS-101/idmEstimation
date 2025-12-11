library(splines2)
#### Create spline hazard ####
make_spline_mat <- function(x, knots, degree, use_bSpline) {

  kn_rn <- range(knots)
  if (min(x) < kn_rn[1] || max(x) > kn_rn[2]) {
    warning("evaluation is value outside knot range, hazard set value in closet knot")

    x[x < kn_rn[1]] = kn_rn[1]
    x[x > kn_rn[2]] = kn_rn[2]
  }

  if(!use_bSpline){

    n_knots <- length(knots)
    i_spline_mat <- splines2::iSpline(x,
      degree = degree,
      knots = knots[-c(1, n_knots)],
      Boundary.knots = knots[c(1, n_knots)],
      intercept = TRUE,
      warn.outside = FALSE,
    )

    m_spline_mat <- splines2::mSpline(x,
      degree = degree,
      knots = knots[-c(1, n_knots)],
      Boundary.knots = knots[c(1, n_knots)],
      intercept = TRUE,
      warn.outside = FALSE,
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
                                      warn.outside = TRUE,
                                      integral = T
    )

    m_spline_mat <- splines2::bSpline(x,
                                      degree = degree,
                                      knots = knots[-c(1, n_knots)],
                                      Boundary.knots = knots[c(1, n_knots)],
                                      intercept = TRUE,
                                      warn.outside = TRUE
    )

    spline_mat_list <- list(
      i_spline = i_spline_mat,
      m_spline = m_spline_mat
    )

  }

  class(spline_mat_list) <- c("spline_mat_list", class(spline_mat_list))
  spline_mat_list
}

setup_case_A_data <- function(data_list, model_config) {
  T_obs <- data_list$T_obs
  V_healthy <- data_list$V_healthy
  V_ill <- data_list$V_ill

  knots_12 <- model_config$knots_12
  knots_13 <- model_config$knots_13
  knots_23 <- model_config$knots_23
  degree <- model_config$degree
  use_bSpline <- ifelse(is.null(model_config$use_bSpline), FALSE, model_config$use_bSpline)

  # Create splines at observation points
  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree, use_bSpline)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree, use_bSpline)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree, use_bSpline)

  # Create integration grid from V_healthy to V_ill
  grid_V_ill <- seq(min(V_healthy), max(V_ill), length.out = 250)
  dx_grid_V_ill <- diff(grid_V_ill)[1]

  spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots_12, degree, use_bSpline)
  spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots_13, degree, use_bSpline)
  spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots_23, degree, use_bSpline)

  list(
    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,
    T_obs_m_spline_mat_23 = spline_T_obs_23$m_spline,

    grid_V_ill_i_spline_mat_12 = spline_grid_V_ill_12$i_spline,
    grid_V_ill_i_spline_mat_13 = spline_grid_V_ill_13$i_spline,
    grid_V_ill_i_spline_mat_23 = spline_grid_V_ill_23$i_spline,
    grid_V_ill_m_spline_mat_12 = spline_grid_V_ill_12$m_spline,

    dx_grid_V_ill = dx_grid_V_ill,

    V_healthy_values = V_healthy,
    V_ill_values = V_ill
  )
}

setup_case_B_data <- function(data_list, model_config) {
  T_obs <- data_list$T_obs
  V_healthy <- data_list$V_healthy
  V_ill <- data_list$V_ill

  knots_12 <- model_config$knots_12
  knots_13 <- model_config$knots_13
  knots_23 <- model_config$knots_23
  degree <- model_config$degree
  use_bSpline <- ifelse(is.null(model_config$use_bSpline), FALSE, model_config$use_bSpline)

  # Create splines at observation points
  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree, use_bSpline)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree, use_bSpline)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree, use_bSpline)

  # Create integration grid from V_healthy to V_ill
  grid_V_ill <- seq(min(V_healthy), max(V_ill), length.out = 250)
  dx_grid_V_ill <- diff(grid_V_ill)[1]

  spline_grid_V_ill_12 <- make_spline_mat(grid_V_ill, knots_12, degree, use_bSpline)
  spline_grid_V_ill_13 <- make_spline_mat(grid_V_ill, knots_13, degree, use_bSpline)
  spline_grid_V_ill_23 <- make_spline_mat(grid_V_ill, knots_23, degree, use_bSpline)

  list(
    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_i_spline_mat_23 = spline_T_obs_23$i_spline,

    grid_V_ill_i_spline_mat_12 = spline_grid_V_ill_12$i_spline,
    grid_V_ill_i_spline_mat_13 = spline_grid_V_ill_13$i_spline,
    grid_V_ill_i_spline_mat_23 = spline_grid_V_ill_23$i_spline,
    grid_V_ill_m_spline_mat_12 = spline_grid_V_ill_12$m_spline,

    dx_grid_V_ill = dx_grid_V_ill,

    V_healthy_values = V_healthy,
    V_ill_values = V_ill
  )
}

setup_case_C_data <- function(data_list, model_config) {

  T_obs <- data_list$T_obs

  knots_12 <- model_config$knots_12
  knots_13 <- model_config$knots_13
  degree <- model_config$degree
  use_bSpline <- ifelse(is.null(model_config$use_bSpline), FALSE, model_config$use_bSpline)

  # Create splines at observation points
  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree, use_bSpline)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree, use_bSpline)

  list(
    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline,
    T_obs_m_spline_mat_13 = spline_T_obs_13$m_spline
  )
}

setup_case_D_data <- function(data_list, model_config) {
  T_obs <- data_list$T_obs
  knots_12 <- model_config$knots_12
  knots_13 <- model_config$knots_13
  degree <- model_config$degree
  use_bSpline <- ifelse(is.null(model_config$use_bSpline), FALSE, model_config$use_bSpline)

  # Create splines at observation points
  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree, use_bSpline)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree, use_bSpline)

  list(
    T_obs_i_spline_mat_12 = spline_T_obs_12$i_spline,
    T_obs_i_spline_mat_13 = spline_T_obs_13$i_spline
  )
}

setup_case_E_data <- function(data_list, model_config) {
  V_healthy <- data_list$V_healthy
  T_obs <- data_list$T_obs
  knots_12 <- model_config$knots_12
  knots_13 <- model_config$knots_13
  knots_23 <- model_config$knots_23
  degree <- model_config$degree
  use_bSpline <- ifelse(is.null(model_config$use_bSpline), FALSE, model_config$use_bSpline)

  # Create splines at observation points
  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree, use_bSpline)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree, use_bSpline)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree, use_bSpline)

  # Create integration grid from V_healthy to T_obs
  grid_T_obs <- seq(min(V_healthy), max(T_obs), length.out = 250)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12, degree, use_bSpline)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13, degree, use_bSpline)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23, degree, use_bSpline)

  list(
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

setup_case_F_data <- function(data_list, model_config) {
  V_healthy <- data_list$V_healthy
  T_obs <- data_list$T_obs
  knots_12 <- model_config$knots_12
  knots_13 <- model_config$knots_13
  knots_23 <- model_config$knots_23
  degree <- model_config$degree
  use_bSpline <- ifelse(is.null(model_config$use_bSpline), FALSE, model_config$use_bSpline)

  # Create splines at observation points
  spline_T_obs_12 <- make_spline_mat(T_obs, knots_12, degree, use_bSpline)
  spline_T_obs_13 <- make_spline_mat(T_obs, knots_13, degree, use_bSpline)
  spline_T_obs_23 <- make_spline_mat(T_obs, knots_23, degree, use_bSpline)

  # Create integration grid from V_healthy to T_obs
  grid_T_obs <- seq(min(V_healthy), max(T_obs), length.out = 250)
  dx_grid_T_obs <- diff(grid_T_obs)[1]

  spline_grid_T_obs_12 <- make_spline_mat(grid_T_obs, knots_12, degree, use_bSpline)
  spline_grid_T_obs_13 <- make_spline_mat(grid_T_obs, knots_13, degree, use_bSpline)
  spline_grid_T_obs_23 <- make_spline_mat(grid_T_obs, knots_23, degree, use_bSpline)

  list(
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



#### Main setup function ####

#' Set Up C++ Model Data Structures
#'
#' Internal function that prepares spline basis matrices and observation data
#' for likelihood computation in C++.
#'
#' @param data_object List of case-specific data from \code{\link{create_case_data}}.
#' @param model_config List with knots, degree, and other model specifications.
#'
#' @return External pointer to C++ model data structure.
#' @keywords internal
setup_cpp_model <- function(data_object,
                            model_config) {

  model_data_list <- list()

  if (!is.null(data_object$case_A)) {
    model_data_list$case_A <- setup_case_A_data(data_object$case_A, model_config)
  }
  if (!is.null(data_object$case_B)) {
    model_data_list$case_B <- setup_case_B_data(data_object$case_B, model_config)
  }
  if (!is.null(data_object$case_C)) {
    model_data_list$case_C <- setup_case_C_data(data_object$case_C, model_config)
  }
  if (!is.null(data_object$case_D)) {
    model_data_list$case_D <- setup_case_D_data(data_object$case_D, model_config)
  }
  if (!is.null(data_object$case_E)) {
    model_data_list$case_E <- setup_case_E_data(data_object$case_E, model_config)
  }
  if (!is.null(data_object$case_F)) {
    model_data_list$case_F <- setup_case_F_data(data_object$case_F, model_config)
  }

  # return pointer to C++ model data
  create_full_data(model_data_list)
}



#' Create Hazard and Distribution Function Estimators
#'
#' Internal function that constructs hazard, cumulative hazard, and distribution
#' function estimators from fitted spline coefficients.
#'
#' @param model_config List with knots, degree, and basis specifications.
#' @param theta_hat_list List with estimated coefficients for each transition.
#'
#' @return List of class \code{"idm_estimators"} with hazard, cumulative hazard,
#'   and distribution functions.
#' @keywords internal
create_estimators <- function(model_config, theta_hat_list) {

  use_bSpline <- ifelse(is.null(model_config$use_bSpline), FALSE, model_config$use_bSpline)
  a12 <- function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_12, model_config$degree,
                                  use_bSpline)$m_spline
    as.vector(spline_mat %*% theta_hat_list$theta_12)
  }
  a13 <- function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_13, model_config$degree,
                                  use_bSpline)$m_spline
    as.vector(spline_mat %*% theta_hat_list$theta_13)
  }
  a23 <- function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_23, model_config$degree,
                                  use_bSpline)$m_spline
    as.vector(spline_mat %*% theta_hat_list$theta_23)
  }
  A12 <- function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_12, model_config$degree,
                                  use_bSpline)$i_spline
    as.vector(spline_mat %*% theta_hat_list$theta_12)
  }
  A13 <- function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_13, model_config$degree,
                                  use_bSpline)$i_spline
    as.vector(spline_mat %*% theta_hat_list$theta_13)
  }
  A23 <- function(x) {
    spline_mat <- make_spline_mat(x, model_config$knots_23, model_config$degree,
                                  use_bSpline)$i_spline
    as.vector(spline_mat %*% theta_hat_list$theta_23)
  }
  # range for approximation grid
  t_max_12 <- max(model_config$knots_12)
  t_max_13 <- max(model_config$knots_13)
  t_max <- max(t_max_12, t_max_13)

  # allow user to override with a larger range if present
  if (!is.null(model_config$t_max)) {
    t_max <- max(t_max, model_config$t_max) + 100
  }

  # number of grid points (can be overridden via model_config$n_grid)
  n_grid <- if (!is.null(model_config$n_grid)) model_config$n_grid else 2000L

  # time grid
  t_grid <- seq(0, t_max, length.out = n_grid)

  # integrands on the grid
  integrand12 <- exp(-A12(t_grid) - A13(t_grid)) * a12(t_grid)
  integrand13 <- exp(-A12(t_grid) - A13(t_grid)) * a13(t_grid)

  # cumulative trapezoidal integration
  dt <- diff(t_grid)
  trap12 <- dt * (head(integrand12, -1) + tail(integrand12, -1)) / 2
  trap13 <- dt * (head(integrand13, -1) + tail(integrand13, -1)) / 2

  F12_vals <- c(0, cumsum(trap12))
  F13_vals <- c(0, cumsum(trap13))

  # standalone, vectorized approximations of the integrals
  F12 <- approxfun(t_grid, F12_vals, rule = 2)  # constant extrapolation
  F13 <- approxfun(t_grid, F13_vals, rule = 2)

  P22 = function(t, entry_time = 0) {
    ifelse(t >= entry_time, exp(-A23(t) + A23(entry_time)), NA_real_)
  }

  estimators <- list(
    hazard_functions = list(
      a12 = a12,
      a13 = a13,
      a23 = a23
    ),
    cum_hazard_functions = list(
      A12 = A12,
      A13 = A13,
      A23 = A23
    ),
    distribution_functions = list(
      F12 = F12,
      F13 = F13,
      P22 = P22
    )
  )

  # add class to hazards
  class(estimators) <- c("idm_estimators", class(estimators))
  estimators
}
