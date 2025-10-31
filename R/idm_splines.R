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

  model_data <- list(
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

  create_penlik_model_data(model_data)
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

  model_data <- list(
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
  create_penlik_model_data(model_data)
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

  model_data <- list(
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
  create_penlik_model_data(model_data)
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

  model_data <- list(
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
  create_penlik_model_data(model_data)
}

#### Main setup function ####

setup_cpp_model <- function(V_0,
                            V_healthy,
                            V_ill,
                            T_obs,
                            status_dead,
                            status_ill,
                            knots_12 = NULL,
                            knots_13 = NULL,
                            knots_23 = NULL,
                            degree = 3) {
  
  # Set default knots if not provided
  if (is.null(knots_12)) {
    knots_12 <- seq(min(V_0), max(T_obs), length.out = 10)
  }
  if (is.null(knots_13)) {
    knots_13 <- seq(min(V_0), max(T_obs), length.out = 10)
  }
  if (is.null(knots_23)) {
    knots_23 <- seq(min(V_0), max(T_obs), length.out = 10)
  }
  
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
  model_pointers_list <- list()
  
  if (length(case_1_idx) > 0) {
    model_pointers_list$case_1 <- setup_case_1_data(
      V_0[case_1_idx],
      V_healthy[case_1_idx],
      T_obs[case_1_idx],
      knots_12, knots_13, knots_23, degree
    )
  }
  
  if (length(case_2_idx) > 0) {
    model_pointers_list$case_2 <- setup_case_2_data(
      V_0[case_2_idx],
      V_healthy[case_2_idx],
      T_obs[case_2_idx],
      knots_12, knots_13, knots_23, degree
    )
  }
  
  if (length(case_3_idx) > 0) {
    model_pointers_list$case_3 <- setup_case_3_data(
      V_0[case_3_idx],
      V_healthy[case_3_idx],
      V_ill[case_3_idx],
      T_obs[case_3_idx],
      knots_12, knots_13, knots_23, degree
    )
  }
  
  if (length(case_4_idx) > 0) {
    model_pointers_list$case_4 <- setup_case_4_data(
      V_0[case_4_idx],
      V_healthy[case_4_idx],
      V_ill[case_4_idx],
      T_obs[case_4_idx],
      knots_12, knots_13, knots_23, degree
    )
  }
  
  model_pointers_list
}

