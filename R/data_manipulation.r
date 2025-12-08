create_case_data <- function(obs_data) {

    # verify data components
  stopifnot(all(c("V_0", "V_healthy", "V_ill", "T_obs", "status_dead", "status_ill") %in% names(obs_data)))

  V_healthy <- obs_data$V_healthy
  V_ill <- obs_data$V_ill
  T_obs <- obs_data$T_obs
  status_dead <- obs_data$status_dead
  status_ill <- obs_data$status_ill

  # Determine case for each observation
  # Case A: status_dead = 1 and status_ill = 1
  # Case B: status_dead = 0 and status_ill = 1
  # Case C: status_dead = 1 and status_ill = 0 and V_healthy == T_obs
  # Case D: status_dead = 0 and status_ill = 0 and V_healthy == T_obs
  # Case E: status_dead = 1 and status_ill = 0 and V_healthy < T_obs
  # Case F: status_dead = 0 and status_ill = 0 and V_healthy < T_obs

  case_A_idx <- which(status_dead == 1 & status_ill == 1) # needs V_healthy, V_ill, T_obs
  case_B_idx <- which(status_dead == 0 & status_ill == 1) # needs V_healthy, V_ill, T_obs
  case_C_idx <- which(status_dead == 1 & status_ill == 0 & abs(V_healthy- T_obs) < 1e-8) # needs T_obs
  case_D_idx <- which(status_dead == 0 & status_ill == 0 & abs(V_healthy- T_obs) < 1e-8) # needs T_obs
  case_E_idx <- which(status_dead == 1 & status_ill == 0 & (T_cutoff - V_healthy) >= 1e-8) # needs V_healthy, T_obs
  case_F_idx <- which(status_dead == 0 & status_ill == 0 & (T_cutoff - V_healthy) >= 1e-8) # needs V_healthy, T_obs


  case_data <- list()

  if (length(case_A_idx) > 0) {
    case_data[["case_A"]] <- list(
      V_healthy = V_healthy[case_A_idx],
      V_ill = V_ill[case_A_idx],
      T_obs = T_obs[case_A_idx])
  }
  if (length(case_B_idx) > 0) {
    case_data[["case_B"]] <- list(
      V_healthy = V_healthy[case_B_idx],
      V_ill = V_ill[case_B_idx],
      T_obs = T_obs[case_B_idx])
  }
  if (length(case_C_idx) > 0) {
    case_data[["case_C"]] <- list(
      T_obs = T_obs[case_C_idx])
  }
  if (length(case_D_idx) > 0) {
    case_data[["case_D"]] <- list(
      T_obs = T_obs[case_D_idx])
  }
  if (length(case_E_idx) > 0) {
    case_data[["case_E"]] <- list(
      V_healthy = V_healthy[case_E_idx],
      T_obs = T_obs[case_E_idx])
  }
  if (length(case_F_idx) > 0) {
    case_data[["case_F"]] <- list(
      V_healthy = V_healthy[case_F_idx],
      T_obs = T_obs[case_F_idx])
  }

  case_data
}




summarise_data_object <- function(data_object) {
  summary_data_object <- list()
  for (case_name in names(data_object)) {
    case_info <- data_object[[case_name]]
    n_cases <- length(case_info$T_obs)
    summary_data_object[[case_name]] <- list(
      n_cases = n_cases,
      V_healthy_range = if ("V_healthy" %in% names(case_info)) range(case_info$V_healthy) else NULL,
      V_ill_range = if ("V_ill" %in% names(case_info)) range(case_info$V_ill) else NULL,
      T_obs_range = range(case_info$T_obs)
    )
  }



  # Last V_ill form A and B and last T obs from E or F
  summary_data_object$last_possible_12 <- if(
      "case_A" %in% names(data_object) || 
      "case_B" %in% names(data_object) || 
      "case_E" %in% names(data_object) || 
      "case_F" %in% names(data_object)) {
    max(c(
      if ("case_A" %in% names(data_object)) max(data_object$case_A$V_ill, na.rm = TRUE) else -Inf,
      if ("case_B" %in% names(data_object)) max(data_object$case_B$V_ill, na.rm = TRUE) else -Inf,
      if ("case_E" %in% names(data_object)) max(data_object$case_E$T_obs, na.rm = TRUE) else -Inf,
      if ("case_F" %in% names(data_object)) max(data_object$case_F$T_obs, na.rm = TRUE) else -Inf
    ), na.rm = TRUE)
  } else {
    NA
  }



  # Last T_obs form C or E
  summary_data_object$last_possible_13 <- if(
      "case_C" %in% names(data_object) || 
      "case_E" %in% names(data_object)) {
    max(c(
      if ("case_C" %in% names(data_object)) max(data_object$case_C$T_obs, na.rm = TRUE) else -Inf,
      if ("case_E" %in% names(data_object)) max(data_object$case_E$T_obs, na.rm = TRUE) else -Inf
    ), na.rm = TRUE)
  } else {
    NA
  }


  # First V_healthy from A and B or E or F
  summary_data_object$first_possible_23 <- if(
      "case_A" %in% names(data_object) ||
      "case_B" %in% names(data_object) ||
      "case_E" %in% names(data_object) ||
      "case_F" %in% names(data_object)) {
    max(c(
      if ("case_A" %in% names(data_object)) min(data_object$case_A$V_ill, na.rm = TRUE) else -Inf,
      if ("case_B" %in% names(data_object)) min(data_object$case_B$V_ill, na.rm = TRUE) else -Inf,
      if ("case_E" %in% names(data_object)) min(data_object$case_E$V_healthy, na.rm = TRUE) else -Inf,
      if ("case_F" %in% names(data_object)) min(data_object$case_F$V_healthy, na.rm = TRUE) else -Inf
    ), na.rm = TRUE)
  } else {
    NA
  }

  # Last T_obs in A, B, F or E
  summary_data_object$last_possible_23 <- if(
      "case_A" %in% names(data_object) || 
      "case_B" %in% names(data_object) || 
      "case_F" %in% names(data_object) || 
      "case_E" %in% names(data_object)) {
    max(c(
      if ("case_A" %in% names(data_object)) max(data_object$case_A$T_obs, na.rm = TRUE) else -Inf,
      if ("case_B" %in% names(data_object)) max(data_object$case_B$T_obs, na.rm = TRUE) else -Inf,
      if ("case_F" %in% names(data_object)) max(data_object$case_F$T_obs, na.rm = TRUE) else -Inf,
      if ("case_E" %in% names(data_object)) max(data_object$case_E$T_obs, na.rm = TRUE) else -Inf
    ), na.rm = TRUE)
  } else {
    NA
  }
  
  summary_data_object
}
