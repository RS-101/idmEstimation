
#' Summarize Simulated Illness-Death Data
#'
#' Computes comprehensive summary statistics for simulated illness-death data,
#' comparing observed patterns with true underlying event times.
#'
#' @param obs_data Data frame with observed (censored) illness-death data containing:
#'   \describe{
#'     \item{status_ill}{Illness observation indicator (1 = observed, 0 = not)}
#'     \item{status_dead}{Death indicator (1 = dead, 0 = censored)}
#'     \item{V_healthy}{Last visit when observed healthy}
#'     \item{T_obs}{Observation or censoring time}
#'     \item{V_0}{Baseline time}
#'   }
#' @param exact_data Data frame with true (uncensored) event times from
#'   \code{\link{simulate_exact_idm}}, containing \code{T_ill}, \code{T_death},
#'   and optionally \code{path}.
#' @param cens_data Optional data frame with censoring mechanism details.
#'   Default is \code{NULL}.
#'
#' @return An object of class \code{"idm_summary"} containing:
#'   \describe{
#'     \item{n}{Total number of observations}
#'     \item{death}{List with death-related counts and proportions}
#'     \item{illness}{List with illness observation statistics}
#'     \item{censoring}{List with censoring counts}
#'     \item{missing_transition}{List with missing transition statistics}
#'     \item{visit_timing}{List with visit interval statistics}
#'     \item{tables}{Contingency tables of observation patterns}
#'   }
#'
#' @details
#' This function compares observed data patterns with true event times to quantify:
#' \itemize{
#'   \item Proportion of illnesses that were observed vs. unobserved
#'   \item Proportion of subjects with missing transition information
#'   \item Death rates and patterns (via illness vs. direct)
#'   \item Visit timing and interval characteristics
#' }
#'
#' Missing transitions occur when a subject was last seen healthy but later
#' dies or is censored, creating uncertainty about whether illness occurred.
#'
#' @examples
#' # Simulate data
#' set.seed(456)
#' sim_data <- simulate_idm_constant_hazards(n = 200)
#'
#' # Typically accessed via summary() method on idm_exact_object:
#' summary(sim_data)
#'
#' # Can also be called directly if you have both exact and observed data:
#' # summary_stats <- summarise_simulated_data(
#' #   obs_data = sim_data$data,
#' #   exact_data = sim_data$exact_idm
#' # )
#'
#' @seealso \code{\link{summarise_obs_data}}, \code{\link{simulate_idm}},
#'   \code{\link{summary.idm_object}}
#' @keywords internal
summarise_simulated_data <- function(obs_data, exact_data, cens_data = NULL) {
  # Validate required columns
  need_cols <- c("status_ill", "status_dead", "V_healthy", "T_obs")
  missing_cols <- setdiff(need_cols, names(obs_data))
  if (length(missing_cols)) {
    stop("Missing columns in observed data: ", paste(missing_cols, collapse = ", "))
  }

  n <- nrow(obs_data)
  observed_illness <- as.logical(obs_data$status_ill)
  observed_death <- as.logical(obs_data$status_dead)
  unknown_transition <- !observed_illness & (obs_data$V_healthy < obs_data$T_obs)

  # Basic counts
  n_any_death <- sum(observed_death, na.rm = TRUE)
  n_death_via_illness <- sum(observed_death & observed_illness, na.rm = TRUE)
  n_death_exact <- sum(observed_death & !unknown_transition & !observed_illness, na.rm = TRUE)
  n_death_missing_transition <- sum(observed_death & unknown_transition & !observed_illness, na.rm = TRUE)
  n_censoring_exact <- sum(!observed_death & !unknown_transition, na.rm = TRUE)
  n_censoring_missing_transition <- sum(!observed_death & unknown_transition & !observed_illness, na.rm = TRUE)
  n_ill_observed <- sum(observed_illness, na.rm = TRUE)
  n_no_ill_obs <- n - n_ill_observed
  n_missing <- n_censoring_missing_transition + n_death_missing_transition

  # Percentages
  p_death <- n_any_death / n
  p_illness_observed <- n_ill_observed / n
  p_missing_transition_given_no_obs <- if (n_no_ill_obs > 0) n_missing / n_no_ill_obs else NA_real_
  p_missing_transition_overall <- n_missing / n

  # From exact data: percentage with illness through state 2
  illness_happened_before_T_obs <- is.finite(exact_data$T_ill) & (exact_data$T_ill <= obs_data$T_obs)
  n_true_illness_before_cutoff <- sum(illness_happened_before_T_obs, na.rm = TRUE)
  n_unobserved_illness <- sum(illness_happened_before_T_obs & !observed_illness, na.rm = TRUE)
  p_illness_unobserved <- if (n_true_illness_before_cutoff > 0) {
    n_unobserved_illness / n_true_illness_before_cutoff
  } else NA_real_

  # Visit timing statistics (only for those with observed illness)
  if (n_ill_observed > 0) {
    ill_idx <- which(observed_illness)
    visit_intervals <- obs_data$V_ill[ill_idx] - obs_data$V_healthy[ill_idx]
    mean_interval_between_visits <- mean(visit_intervals, na.rm = TRUE)
    sd_interval_between_visits <- sd(visit_intervals, na.rm = TRUE)
  } else {
    mean_interval_between_visits <- NA_real_
    sd_interval_between_visits <- NA_real_
  }

  # Time between last visit and cutoff for non-observed illness
  if (n_no_ill_obs > 0) {
    no_ill_idx <- which(!observed_illness)
    time_to_cutoff <- obs_data$T_obs[no_ill_idx] - obs_data$V_healthy[no_ill_idx]
    mean_time_last_visit_to_cutoff <- mean(time_to_cutoff, na.rm = TRUE)
    sd_time_last_visit_to_cutoff <- sd(time_to_cutoff, na.rm = TRUE)
  } else {
    mean_time_last_visit_to_cutoff <- NA_real_
    sd_time_last_visit_to_cutoff <- NA_real_
  }

  # True path classification
  true_path <- rep(NA_character_, n)
  if (!is.null(exact_data$path)) {
    true_path <- as.character(exact_data$path)
  }

  # Build result
  summary_result <- list(
    n = n,

    # Death statistics
    death = list(
      n_any_death = n_any_death,
      n_death_via_illness = n_death_via_illness,
      n_death_exact = n_death_exact,
      n_death_missing_transition = n_death_missing_transition,
      p_death = p_death
    ),

    # Illness statistics
    illness = list(
      n_ill_observed = n_ill_observed,
      n_no_ill_obs = n_no_ill_obs,
      n_true_illness_before_cutoff = n_true_illness_before_cutoff,
      n_unobserved_illness = n_unobserved_illness,
      p_illness_observed = p_illness_observed,
      p_illness_unobserved = p_illness_unobserved
    ),

    # Censoring statistics
    censoring = list(
      n_censoring_exact = n_censoring_exact,
      n_censoring_missing_transition = n_censoring_missing_transition
    ),

    # Missing transition statistics
    missing_transition = list(
      n_missing = n_missing,
      p_missing_given_no_obs = p_missing_transition_given_no_obs,
      p_missing_overall = p_missing_transition_overall
    ),

    # Visit timing
    visit_timing = list(
      mean_interval_between_visits = mean_interval_between_visits,
      sd_interval_between_visits = sd_interval_between_visits,
      mean_time_last_visit_to_cutoff = mean_time_last_visit_to_cutoff,
      sd_time_last_visit_to_cutoff = sd_time_last_visit_to_cutoff
    ),

    # Cross-tabulation
    tables = list(
      observed_status = table(
        Illness = observed_illness,
        Death = observed_death,
        useNA = "ifany"
      ),
      true_vs_observed = if (!all(is.na(true_path))) {
        table(
          TruePath = true_path,
          ObservedIllness = observed_illness,
          useNA = "ifany"
        )
      } else NULL
    )
  )

  class(summary_result) <- c("idm_summary", class(summary_result))
  summary_result
}

#' Print Method for IDM Summary Objects
#'
#' Prints a formatted summary of illness-death model simulation statistics
#' comparing observed data with true event times.
#'
#' @param x An object of class \code{"idm_summary"} from
#'   \code{\link{summarise_simulated_data}} or \code{\link{summary.idm_object}}.
#' @param ... Additional arguments (not currently used).
#'
#' @return Invisibly returns \code{x}.
#'
#' @details
#' Prints a comprehensive summary including:
#' \itemize{
#'   \item Sample size
#'   \item Death statistics (total, via illness, direct, with missing transitions)
#'   \item Illness observation statistics (including unobserved illness)
#'   \item Missing transition proportions
#'   \item Visit timing statistics
#'   \item Cross-tabulations of observation patterns
#' }
#'
#' This print method is automatically called when printing the result of
#' \code{summary()} on an \code{"idm_exact_object"}.
#'
#' @seealso \code{\link{summary.idm_object}}, \code{\link{print.idm_obs_summary}}
#' @export
print.idm_summary <- function(x, ...) {
  cat("=== IDM Simulation Summary ===\n")
  cat(sprintf("Sample size: %d\n\n", x$n))

  cat("--- Death Statistics ---\n")
  cat(sprintf("  Total deaths: %d (%.1f%%)\n",
              x$death$n_any_death, 100 * x$death$p_death))
  cat(sprintf("  Deaths via illness: %d\n", x$death$n_death_via_illness))
  cat(sprintf("  Direct deaths (exact): %d\n", x$death$n_death_exact))
  cat(sprintf("  Deaths with missing transition: %d\n\n",
              x$death$n_death_missing_transition))

  cat("--- Illness Statistics ---\n")
  cat(sprintf("  Illness observed: %d (%.1f%%)\n",
              x$illness$n_ill_observed, 100 * x$illness$p_illness_observed))
  cat(sprintf("  True illness before cutoff: %d\n",
              x$illness$n_true_illness_before_cutoff))
  cat(sprintf("  Unobserved illness: %d (%.1f%% of true illness)\n\n",
              x$illness$n_unobserved_illness,
              100 * ifelse(is.na(x$illness$p_illness_unobserved), 0,
                           x$illness$p_illness_unobserved)))

  cat("--- Missing Transitions ---\n")
  cat(sprintf("  Total missing transitions: %d\n", x$missing_transition$n_missing))
  cat(sprintf("  Among no observed illness: %.1f%%\n",
              100 * x$missing_transition$p_missing_given_no_obs))
  cat(sprintf("  Overall: %.1f%%\n\n",
              100 * x$missing_transition$p_missing_overall))

  cat("--- Visit Timing ---\n")
  if (!is.na(x$visit_timing$mean_interval_between_visits)) {
    cat(sprintf("  Mean interval between visits (observed illness): %.2f (SD: %.2f)\n",
                x$visit_timing$mean_interval_between_visits,
                x$visit_timing$sd_interval_between_visits))
  }
  if (!is.na(x$visit_timing$mean_time_last_visit_to_cutoff)) {
    cat(sprintf("  Mean time last visit to cutoff (no obs illness): %.2f (SD: %.2f)\n",
                x$visit_timing$mean_time_last_visit_to_cutoff,
                x$visit_timing$sd_time_last_visit_to_cutoff))
  }

  cat("\n--- Observed Status Table ---\n")
  print(x$tables$observed_status)

  if (!is.null(x$tables$true_vs_observed)) {
    cat("\n--- True Path vs Observed Illness ---\n")
    print(x$tables$true_vs_observed)
  }

  invisible(x)
}

#' Summarize Observed Illness-Death Data
#'
#' Computes basic summary statistics for observed illness-death data without
#' access to true underlying event times.
#'
#' @param obs_data Data frame with observed illness-death data containing:
#'   \describe{
#'     \item{status_ill}{Illness observation indicator (1 = observed, 0 = not)}
#'     \item{status_dead}{Death indicator (1 = dead, 0 = censored)}
#'     \item{V_healthy}{Last visit when observed healthy}
#'     \item{T_obs}{Observation or censoring time}
#'     \item{V_0}{Baseline time}
#'   }
#'
#' @return An object of class \code{"idm_obs_summary"} containing:
#'   \describe{
#'     \item{n}{Total number of observations}
#'     \item{n_deaths}{Number of observed deaths}
#'     \item{n_illness_observed}{Number of subjects with observed illness}
#'     \item{n_censored}{Number of censored subjects}
#'     \item{p_death}{Proportion of deaths}
#'     \item{p_illness}{Proportion with observed illness}
#'     \item{status_table}{Contingency table of illness and death status}
#'   }
#'
#' @details
#' This function provides summary statistics based only on observed data,
#' without comparing to true event times. Use \code{\link{summarise_simulated_data}}
#' when both observed and exact data are available.
#'
#' @examples
#' # Simulate data
#' set.seed(789)
#' sim_data <- simulate_idm_constant_hazards(n = 200)
#'
#' # Summarize observed data only
#' obs_summary <- summarise_obs_data(sim_data$data)
#' print(obs_summary)
#'
#' @seealso \code{\link{summarise_simulated_data}}
#' @keywords internal
summarise_obs_data <- function(obs_data) {
  # Validate required columns
  need_cols <- c("status_ill", "status_dead")
  missing_cols <- setdiff(need_cols, names(obs_data))
  if (length(missing_cols)) {
    stop("Missing columns in observed data: ", paste(missing_cols, collapse = ", "))
  }

  n <- nrow(obs_data)
  observed_illness <- as.logical(obs_data$status_ill)
  observed_death <- as.logical(obs_data$status_dead)

  n_deaths <- sum(observed_death, na.rm = TRUE)
  n_illness_observed <- sum(observed_illness, na.rm = TRUE)
  n_censored <- n - n_deaths

  p_death <- n_deaths / n
  p_illness <- n_illness_observed / n

  status_table <- table(
    Illness = observed_illness,
    Death = observed_death,
    useNA = "ifany"
  )

  result <- list(
    n = n,
    n_deaths = n_deaths,
    n_illness_observed = n_illness_observed,
    n_censored = n_censored,
    p_death = p_death,
    p_illness = p_illness,
    status_table = status_table
  )

  class(result) <- c("idm_obs_summary", class(result))
  result
}

#' Print Method for Observed IDM Summary
#'
#' Prints a formatted summary of observed illness-death data.
#'
#' @param x An object of class \code{"idm_obs_summary"} from
#'   \code{\link{summarise_obs_data}}.
#' @param ... Additional arguments (not currently used).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.idm_obs_summary <- function(x, ...) {
  cat("=== IDM Observed Data Summary ===\n")
  cat(sprintf("Sample size: %d\n\n", x$n))

  cat("--- Observed Events ---\n")
  cat(sprintf("  Deaths: %d (%.1f%%)\n", x$n_deaths, 100 * x$p_death))
  cat(sprintf("  Illness observed: %d (%.1f%%)\n", x$n_illness_observed, 100 * x$p_illness))
  cat(sprintf("  Censored: %d (%.1f%%)\n\n", x$n_censored, 100 * (1 - x$p_death)))

  cat("--- Status Table ---\n")
  print(x$status_table)

  invisible(x)
}

#' Summary Method for IDM Objects
#'
#' Produces summary statistics for illness-death model objects. For objects
#' with class \code{"idm_exact_object"}, compares observed data with true
#' event times. For regular \code{"idm_object"}, summarizes only observed data.
#'
#' @param object An object of class \code{"idm_object"} from simulation functions
#'   like \code{\link{simulate_idm_constant_hazards}} or fitted models.
#' @param ... Additional arguments (not currently used).
#'
#' @return An object of class \code{"idm_summary"} if exact data is available,
#'   otherwise \code{"idm_obs_summary"}.
#'
#' @details
#' For simulated data with exact event times (\code{"idm_exact_object"}), this
#' computes comprehensive statistics comparing observed patterns with true
#' events using \code{\link{summarise_simulated_data}}.
#'
#' For fitted models or observed data only (\code{"idm_object"}), this provides
#' basic summary statistics using \code{\link{summarise_obs_data}}.
#'
#' @examples
#' # Simulated data with exact event times
#' set.seed(123)
#' sim_data <- simulate_idm_constant_hazards(n = 200)
#' summary(sim_data)  # Full comparison with exact data
#'
#' @seealso \code{\link{summarise_simulated_data}}, \code{\link{summarise_obs_data}}
#' @export
summary.idm_object <- function(object, ...) {
  if (inherits(object, "idm_exact_object") && !is.null(object$exact_idm)) {
    # Full summary with exact data comparison
    summarise_simulated_data(
      obs_data = object$data,
      exact_data = object$exact_idm,
      cens_data = if (!is.null(object$model_config$observation_scheme)) {
        object$model_config$observation_scheme
      } else NULL
    )
  } else {
    # Observed data summary only
    summarise_obs_data(object$data)
  }
}

