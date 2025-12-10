# case classification:
  # Case A: illness and death observed
  # Case B: illness observed and censored while ill
  # Case C: death observed at state 1 without illness V_healthy == T_obs
  # Case D: censored at state 1 without illness V_healthy == T_obs
  # Case E: not observed and death observed either at state 1 or 2
  # Case F: not observed and censored either at state 1 or 2

#' Simulate Exact Illness-Death Model Data
#'
#' Simulates exact (uncensored) event times from an illness-death model with
#' arbitrary time-dependent hazard functions. Uses inverse transform sampling
#' with numerical integration to generate transition times.
#'
#' @param n Integer. Number of subjects to simulate.
#' @param a12 Function of one argument (absolute time \code{t}) giving the hazard
#'   for the 1→2 (healthy to illness) transition.
#' @param a13 Function of one argument (absolute time \code{t}) giving the hazard
#'   for the 1→3 (healthy to death) transition.
#' @param a23 Function of one argument (absolute time \code{t}) giving the hazard
#'   for the 2→3 (illness to death) transition. This is a calendar-time hazard.
#' @param t0 Numeric. Starting time (entry time) for all subjects. Default is 0.
#' @param tmax Numeric. Maximum follow-up time. Events occurring after \code{tmax}
#'   are treated as never occurring (set to \code{Inf}). Default is \code{Inf}.
#' @param init_step Numeric. Initial step size for the bracketing search in
#'   inversion sampling. Default is 1.
#' @param int_abs_tol Numeric. Absolute tolerance for numerical integration
#'   (passed to \code{\link[stats]{integrate}}). Default is 1e-8.
#' @param root_tol Numeric. Tolerance for root finding (passed to
#'   \code{\link[stats]{uniroot}}). Default is 1e-8.
#' @param max_doublings Integer. Maximum number of bracket doublings allowed
#'   in the search for event times. Default is 60.
#'
#' @return A data frame with class \code{"exact_idm"} containing:
#'   \describe{
#'     \item{id}{Subject ID (1 to \code{n}).}
#'     \item{T_ill}{Time of illness (transition 1→2). \code{Inf} if direct death occurred.}
#'     \item{T_death}{Time of death (transition to state 3).}
#'     \item{path}{Factor indicating the path: \code{"1->3"} (direct death) or
#'       \code{"1->2->3"} (illness then death).}
#'   }
#'
#' @details
#' The function simulates competing risks from state 1 (healthy) using calendar-time
#' hazards \code{a12(t)} and \code{a13(t)}. If illness occurs first, the death time
#' from state 2 is drawn using the calendar-time hazard \code{a23(t)} starting from
#' the illness time.
#'
#' Event times are generated via inverse transform sampling: for a hazard \code{h(t)},
#' the event time \code{T} satisfies \eqn{\int_{t_0}^T h(u) du = E}, where \code{E ~ Exp(1)}.
#' This integral equation is solved using numerical integration and root finding.
#'
#' @examples
#' # Constant hazards
#' set.seed(123)
#' exact_data <- simulate_exact_idm(
#'   n = 100,
#'   a12 = function(t) rep(0.01, length(t)),
#'   a13 = function(t) rep(0.005, length(t)),
#'   a23 = function(t) rep(0.02, length(t))
#' )
#' head(exact_data)
#' table(exact_data$path)
#'
#' @seealso \code{\link{add_censoring_type_1}}, \code{\link{simulate_idm}}
#'
#' @export
simulate_exact_idm <- function(
  n,
  a12, a13, a23, # all are functions of ABSOLUTE calendar time t
  t0 = 0,
  tmax = Inf, # if finite: events after tmax are treated as never occurring (Inf)
  init_step = 1, # initial bracket step for inversion search
  int_abs_tol = 1e-8, # absolute tolerance for integrate()
  root_tol = 1e-8, # tolerance for uniroot()
  max_doublings = 60 # safety cap for growing the search bracket
) {
  # --- Safeguards ----------
  if (length(formals(a12)) != 1L || length(formals(a13)) != 1L) {
    stop("a12 and a13 must be functions of one argument: absolute time t.")
  }
  if (length(formals(a23)) != 1L) {
    stop("a23 must be a function of one argument: absolute time t (calendar-time hazard).")
  }

  # --- Helper: cumulative hazard from t_start to t_end ---------
  cumhaz <- function(h, t_start, t_end) {
    if (t_end <= t_start) {
      return(0)
    }
    res <- stats::integrate(function(u) h(u),
      lower = t_start, upper = t_end,
      stop.on.error = TRUE, abs.tol = int_abs_tol
    )
    as.numeric(res$value)
  }

  # --- Helper: draw absolute event time by inverting the cumulative hazard ----
  # Given hazard h(t) and starting clock t_start, draw T ≥ t_start with
  #  P(T > t | T ≥ t_start) = exp(-∫_{t_start}^t h(u) du).
  draw_time_from_hazard <- function(h, t_start) {
    target <- stats::rexp(1) # -log(U) ~ Exp(1)

    # grow an upper bound ub until ∫ h >= target, or we hit tmax
    lb <- t_start
    step <- init_step
    ub <- lb + step
    ch <- cumhaz(h, lb, ub)

    n_doubles <- 0
    while (is.finite(ch) && ch < target && ub < tmax && n_doubles < max_doublings) {
      lb <- ub
      step <- step * 2
      ub <- lb + step
      ch <- ch + cumhaz(h, lb, ub) # accumulate ∫ on the new segment
      n_doubles <- n_doubles + 1
    }

    # if we couldn’t reach target before tmax, the event never occurs
    if (!is.finite(ch) || (ub >= tmax && ch < target)) {
      return(Inf)
    }

    # root-find F(t) = ∫_{t_start}^{t} h(u) du - target = 0 on [lb0, ub0]
    F <- function(t) cumhaz(h, t_start, t) - target
    lb0 <- t_start
    ub0 <- ub
    # F(lb0) < 0 always (target > 0); ensure F(ub0) >= 0
    val_ub <- F(ub0)
    while (val_ub < 0 && ub0 < tmax && n_doubles < max_doublings + 20) {
      lb0 <- ub0
      ub0 <- ub0 + step
      val_ub <- F(ub0)
      n_doubles <- n_doubles + 1
    }
    if (val_ub < 0) {
      return(Inf)
    }

    uniroot(F, lower = lb0, upper = ub0, tol = root_tol)$root # absolute time
  }

  # --- Main simulation loop ---------------------------------------------------
  out <- vector("list", n)
  for (i in seq_len(n)) {
    # Draw absolute times from state 1 (competing risks, calendar time)
    # Time-change theorem: if E ~ Exp(1), then T = inf{ t ≥ t0 : ∫_{t0}^t a_{1k}(u) du ≥ E }
    t12 <- draw_time_from_hazard(a12, t0)
    t13 <- draw_time_from_hazard(a13, t0)

    if (t12 <= t13) {
      # Illness occurs first at entry2 = t12

      # From state 2, with calendar-time hazard a23(t):
      # Conditional on entry at u, T | U=u has survival exp(-∫_u^T a23(v) dv).
      # We sample the ABSOLUTE death time:
      death_time <- draw_time_from_hazard(a23, t12)

      out[[i]] <- list(
        id = i,
        T_ill = t12,
        T_death = death_time,
        path = "1->2->3"
      )
    } else {
      # Direct death from state 1 at t13
      out[[i]] <- list(
        id = i,
        T_ill = Inf,
        T_death = t13,
        path = "1->3"
      )
    }
  }

  df <- do.call(rbind, lapply(out, as.data.frame))
  rownames(df) <- NULL
  df$path <- factor(df$path, levels = c("1->3", "1->2->3"))
  class(df) <- c(class(df), "exact_idm")
  df
}


#' Add Interval Censoring to Exact Illness-Death Data
#'
#' Each subject has a visiting scheme generated by a Gamma distribution with shape
#'  and scale parameters derived from the mean and standard deviation of time
#'  between visits. Right censoring is sampled with a binomial with
#'  probability_of_right_censoring. For right censored subjects, the censoring time
#'  is drawn from a uniform distribution between 0 and T_obs.
#'
#' @param exact_idm Data frame with class \code{"exact_idm"} from
#'   \code{\link{simulate_exact_idm}}, containing true event times.
#' @param mean_time_between_visits Numeric. Average time between visits
#' @param sd_time_between_visits Numeric. sd of time between visits
#' @param max_right_censoring_time Numeric. Maximum time for right censoring by uniform variable on 0 to max_right_censoring_time
#'
#' @return A list with two components:
#'   \describe{
#'     \item{obs}{Data frame with observed (censored) data containing:
#'       \itemize{
#'         \item \code{id}: Subject ID
#'         \item \code{V_healthy}: Last visit when observed healthy
#'         \item \code{V_ill}: First visit when illness observed (NA if never observed)
#'         \item \code{T_obs}: Observation/censoring time
#'         \item \code{status_dead}: Indicator for observed death (1) or censoring (0)
#'         \item \code{status_ill}: Indicator for observed illness (1) or not (0)
#'         \item \code{case}: Factor classifying observation pattern (A-F)
#'       }}
#'     \item{cens_mechanism}{Data frame with censoring mechanism details including
#'       observation schedules and censoring times.}
#'   }
#'
#' @details
#'
#' Six observation cases are defined:
#' \itemize{
#'   \item Case A: Illness and death both observed
#'   \item Case B: Illness observed, censored while ill
#'   \item Case C: Death observed, healthy at last visit before death
#'   \item Case D: Censored, healthy at last visit
#'   \item Case E: Death observed, unknown if healthy or ill (missing transition)
#'   \item Case F: Censored, unknown if healthy or ill (missing transition)
#' }
#'
#' @export
add_censoring_type_1 <- function(exact_idm,
                                 mean_time_between_visits,
                                 sd_time_between_visits,
                                 max_right_censoring_time,
                                 prob_censoring_at_last_visit
) {

  stopifnot("exact_idm" %in% class(exact_idm))
  # Extract core times
  id <- exact_idm$id
  time_to_illness <- as.numeric(exact_idm$T_ill)
  time_to_death <- as.numeric(exact_idm$T_death)
  n <- length(time_to_illness)

  # Random right censoring
  censoring_time <- stats::runif(n, 0, max_right_censoring_time)
  status_dead <- ifelse(time_to_death <= censoring_time, 1, 0)

  T_obs <- pmin(censoring_time, time_to_death)
  max_number_of_visits <- ceiling(max(T_obs) / mean_time_between_visits) + 10L

  # Create observation schedule with visits
  obs_schedule <- matrix(NA_real_, nrow = n, ncol = max_number_of_visits) # preallocate large matrix

  # first visit at time 0
  obs_schedule[, 1] <- 0

  # we sample from a gamma distribution. reparameterized by shape and scale
  shape <- (mean_time_between_visits^2) / (sd_time_between_visits^2)
  scale <- (sd_time_between_visits^2) / mean_time_between_visits

  # sample all varialbes and take the cumulative sum row-wise
  inter_visit_times <- matrix(stats::rgamma(n * (max_number_of_visits - 1), shape = shape, scale = scale),
         nrow = n, ncol = (max_number_of_visits - 1)
  )

  obs_schedule[, 2:max_number_of_visits] <- t(apply(inter_visit_times, 1, cumsum))


  # ---- Determine last healthy visit before illness (or before cutoff if earlier)
  # We use T_end = min(illness time, cutoff) row-wise
  T_end <- pmin(time_to_illness, T_obs)
  # Count visits strictly before T_end, row-wise
  idx_healthy <- rowSums(sweep(obs_schedule, 1, T_end, "<"))

  V_healthy <- rep(0, n)
  sel <- idx_healthy > 0
  V_healthy[sel] <- obs_schedule[cbind(which(sel), idx_healthy[sel])]

  # ---- First visit after the last healthy visit (candidate "ill" visit)
  idx_ill <- pmin(idx_healthy + 1L, max_number_of_visits)
  V_ill_candidate <- obs_schedule[cbind(seq_len(n), idx_ill)]

  # Illness occurs before cutoff?
  ill_before_cutoff <- is.finite(time_to_illness) & (time_to_illness <= T_obs)

  # Candidate ill-visit is valid only if it happens on/before cutoff and there is a next visit
  has_next_visit_before_cutoff <- ill_before_cutoff &
    (idx_healthy < max_number_of_visits) &
    (V_ill_candidate <= T_obs + 1e-12)

  V_ill <- ifelse(has_next_visit_before_cutoff, V_ill_candidate, NA_real_)
  status_ill <- !is.na(V_ill)


  # For censored subjects, with probability prob_censoring_at_last_visit,
  # set T_obs to last visit time if earlier than current T_obs
  idx_censored <- ifelse(status_dead == 0 & rbinom(n, 1, prob_censoring_at_last_visit) == 1, 1, 0)

  T_obs <- ifelse(
    idx_censored == 1 & status_ill == 0,
    V_healthy,
    T_obs
  )



  status <- character(n)
  # Case A illness and death observed
  # Case B illness observed and censored while ill
  # Case C death observed at state 1 without illness V_healthy == T_obs
  # Case D censored at state 1 without illness V_healthy == T_obs
  # Case E not observed and death observed either at state 1 or 2
  # Case F not observed and censored either at state 1 or 2

  status[status_ill & status_dead] <- "A"
  status[status_ill & !status_dead] <- "B"
  status[!status_ill & status_dead & abs(V_healthy - T_obs) < 1e-8] <- "C"
  status[!status_ill & !status_dead & abs(V_healthy - T_obs) < 1e-8] <- "D"
  status[!status_ill & status_dead & (T_obs- V_healthy) >= 1e-8] <- "E"
  status[!status_ill & !status_dead & (T_obs - V_healthy) >= 1e-8] <- "F"

  # Return
  df_obs_idm <- data.frame(
    id = id,
    V_healthy = V_healthy,
    V_ill = V_ill,
    T_obs = T_obs,
    status_dead = as.numeric(status_dead),
    status_ill = as.numeric(status_ill),
    case = factor(
      status,
      levels = c("A", "B", "C", "D", "E", "F")
    )
  )

  # Store the (possibly wide) schedule as a list-column to avoid unintended column expansion
  df_observation_scheme <- data.frame(
    id = id,
    obs_schedule = obs_schedule
  )

  list(obs = df_obs_idm, cens_mechanism = df_observation_scheme)
}

#' Add Interval Censoring to Exact Illness-Death Data
#'
#' Applies type 2 censoring
#'
#' At each visit, subject have a probability of dropping out of the study,
#'  such that they are no longer followed detemind by probability_of_dropout.
#'  The dropout probability increases at each visit by a factor of
#'  increment_in_dropout_prop.  The time between visits is generated from a
#'  gamma distribution with shape and scale parameters derived from the mean
#'  and standard deviation of time between visits.

#' @param exact_idm Data frame with class \code{"exact_idm"} from
#'   \code{\link{simulate_exact_idm}}, containing true event times.
#' @param probability_of_dropout Numeric. Probability of dropout at each visit
#' @param mean_time_between_visits Numeric. Average time between visits
#' @param sd_time_between_visits Numeric. sd of time between visits
#' @param increment_in_dropout_prop Numeric. Relative increment in probability of
#'  dropout at each visit
#' @return A list with two components:
#'   \describe{
#'     \item{obs}{Data frame with observed (censored) data containing:
#'       \itemize{
#'         \item \code{id}: Subject ID
#'         \item \code{V_healthy}: Last visit when observed healthy
#'         \item \code{V_ill}: First visit when illness observed (NA if never observed)
#'         \item \code{T_obs}: Observation/censoring time
#'         \item \code{status_dead}: Indicator for observed death (1) or censoring (0)
#'         \item \code{status_ill}: Indicator for observed illness (1) or not (0)
#'         \item \code{case}: Factor classifying observation pattern (A-F)
#'       }}
#'     \item{cens_mechanism}{Data frame with censoring mechanism details including
#'       observation schedules and censoring times.}
#'   }
#'
#' @details
#' The function creates an observation schedule with visits spaced to achieve
#' approximately \code{average_number_of_visits} over each subject's follow-up.
#' Visits are jittered with small random noise to mimic realistic data.
#'
#' Censoring times are drawn from Uniform(0, 3 * death_time) for each subject.
#' Illness is only observed if it occurs between two consecutive visits within
#' the observation period.
#'
#' Six observation cases are defined:
#' \itemize{
#'   \item Case A: Illness and death both observed
#'   \item Case B: Illness observed, censored while ill
#'   \item Case C: Death observed, healthy at last visit before death
#'   \item Case D: Censored, healthy at last visit
#'   \item Case E: Death observed, unknown if healthy or ill (missing transition)
#'   \item Case F: Censored, unknown if healthy or ill (missing transition)
#' }
#'
#' @export
add_censoring_type_2 <- function(exact_idm,
                                 probability_of_dropout,
                                 mean_time_between_visits,
                                 sd_time_between_visits,
                                 increment_in_dropout_prop,
                                 max_right_censoring_time,
                                 prob_censoring_at_last_visit) {

  stopifnot("exact_idm" %in% class(exact_idm))
  stopifnot(mean_time_between_visits > 0,
            sd_time_between_visits > 0,
            probability_of_dropout >= 0,
            probability_of_dropout <= 1,
            prob_censoring_at_last_visit >= 0,
            prob_censoring_at_last_visit <= 1)

  ## Extract core times
  id <- exact_idm$id
  time_to_illness <- as.numeric(exact_idm$T_ill)
  time_to_death <- as.numeric(exact_idm$T_death)
  n <- length(time_to_illness)

  ## Random right censoring
  ## status_dead: 1 = death observed, 0 = right censored

  censoring_time <- stats::runif(n, 0, max_right_censoring_time)
  status_dead <- ifelse(time_to_death <= censoring_time, 1, 0)

  ## Observation time: if death observed -> time_to_death,
  ## otherwise uniform on (0, time_to_death)
  T_obs <- pmin(censoring_time, time_to_death)

  ## Max number of visits (sufficiently large upper bound)
  max_number_of_visits <- ceiling(max(time_to_death) / mean_time_between_visits) + 10L

  ## Observation schedule matrix
  obs_schedule <- matrix(NA_real_, nrow = n, ncol = max_number_of_visits)
  obs_schedule[, 1] <- 0  # first visit at time 0

  ## Gamma parameters for inter-visit times
  shape <- (mean_time_between_visits^2) / (sd_time_between_visits^2)
  scale <- (sd_time_between_visits^2) / mean_time_between_visits

  ## Sample inter-visit times and cumulate row-wise
  inter_visit_times <- matrix(
    stats::rgamma(n * (max_number_of_visits - 1),
                  shape = shape, scale = scale),
    nrow = n, ncol = (max_number_of_visits - 1)
  )

  obs_schedule[, 2:max_number_of_visits] <- t(apply(inter_visit_times, 1, cumsum))

  ## Apply dropout at each visit
  for (visit in 2:max_number_of_visits) {
    dropout_indicator <- stats::rbinom(n, 1, probability_of_dropout)
    dropped_out <- which(dropout_indicator == 1)
    if (length(dropped_out) > 0) {
      obs_schedule[dropped_out, visit:max_number_of_visits] <- NA_real_
    }
    ## Increment dropout probability, but cap at 1
    probability_of_dropout <- min(1, probability_of_dropout * (1 + increment_in_dropout_prop))
  }


  ## Use cutoff = T_obs (already equals time_to_death for those who die)
  T_end <- pmin(time_to_illness, T_obs)

  ## Count visits strictly before T_end, ignoring NAs from dropout
  idx_healthy <- rowSums(sweep(obs_schedule, 1, T_end, "<"), na.rm = TRUE)

  V_healthy <- numeric(n)
  sel <- idx_healthy > 0
  if (any(sel)) {
    V_healthy[sel] <- obs_schedule[cbind(which(sel), idx_healthy[sel])]
  }

  ## First visit after last healthy visit (candidate illness-visit)
  idx_ill <- pmin(idx_healthy + 1L, max_number_of_visits)
  V_ill_candidate <- obs_schedule[cbind(seq_len(n), idx_ill)]

  ## Illness occurs before the observation cutoff?
  ill_before_cutoff <- is.finite(time_to_illness) & (time_to_illness <= T_obs + 1e-12)

  ## Candidate ill-visit valid if it is on/before cutoff and there is a next visit
  has_next_visit_before_cutoff <- ill_before_cutoff &
    (idx_healthy < max_number_of_visits) &
    (V_ill_candidate <= T_obs + 1e-12)

  V_ill <- ifelse(has_next_visit_before_cutoff, V_ill_candidate, NA_real_)

  status_ill <- !is.na(V_ill)

  ## (We only act on censored subjects; probability = prob_censoring_at_last_visit)

  idx_censored <- ifelse(
    status_dead == 0 & status_ill == 0 &
      stats::rbinom(n, 1, prob_censoring_at_last_visit) == 1,
    1,
    0
  )

  T_obs <- ifelse(
    idx_censored == 1,
    V_healthy,   # last healthy visit time (always <= original T_obs)
    T_obs
  )

  status <- character(n)
  # Case A: illness and death observed
  status[status_ill & status_dead == 1] <- "A"
  # Case B: illness observed and censored while ill
  status[status_ill & status_dead == 0] <- "B"
  # Case C: death observed, last visit healthy (V_healthy == T_obs)
  status[!status_ill & status_dead == 1 & abs(V_healthy - T_obs) < 1e-8] <- "C"
  # Case D: censored, last visit healthy
  status[!status_ill & status_dead == 0 & abs(V_healthy - T_obs) < 1e-8] <- "D"
  # Case E: death observed, last visit not exactly at T_obs (missing transition)
  status[!status_ill & status_dead == 1 & (T_obs - V_healthy) >= 1e-8] <- "E"
  # Case F: censored, last visit not exactly at T_obs (missing transition)
  status[!status_ill & status_dead == 0 & (T_obs - V_healthy) >= 1e-8] <- "F"

  ## ---- Output data frames ----

  df_obs_idm <- data.frame(
    id = id,
    V_healthy = V_healthy,
    V_ill = V_ill,
    T_obs = T_obs,
    status_dead = as.numeric(status_dead),
    status_ill = as.numeric(status_ill),
    case = factor(status, levels = c("A", "B", "C", "D", "E", "F"))
  )

  ## If you really want a list-column with each subject's visit times, you can do:
  ## obs_schedule_list <- split(as.data.frame(t(obs_schedule)), seq_len(n))
  ## but here I'll keep a wide matrix as-is.
  df_observation_scheme <- data.frame(
    id = id,
    obs_schedule = obs_schedule
  )

  list(
    obs = df_obs_idm,
    cens_mechanism = df_observation_scheme
  )
}


add_censoring_frydman <- function(exact_idm, scenario = 1L) {
  # The illness–death model with constant intensities, 0.0008, 0.0002, and 0.0016 of 1 → 2, 1 → 3, and 2 → 3 transitions,
  # respectively, under four different scenarios where each observation had a potential maximum follow-up time of approximately 4 years.
  # Under these parameters, this number of simulations is sufficient to estimate F(s) with 5% accuracy.
  #
  # In scenarios 1 and 2, the probability of unknown nonfatal event status was 0.5% for the first post baseline assessment, doubling at every assessment thereafter. In scenarios 3 and 4, the probability of unknown nonfatal event status was 3% for the first postbaseline assessment, increasing twofold at every assessment thereafter.
  #
  # In scenarios 1 and 3, the nonfatal event was assessed every 6 months, ±N(0, 100) days.
  # In scenarios 2 and 4, the nonfatal event was assessed every 2 months for the first 6 months, every 6 months through 2 years, and at years 3 and 4, ±N(0, 100) days.
  #
  # Nonfatal event status was unknown for a mean of 31%, 41%, 36%, and 48% of the observations in scenarios 1–4, respectively.
  # Vital status was known for all observations at the end of the follow-up period.


  stopifnot("exact_idm" %in% class(exact_idm))

  n <- nrow(exact_idm)

  id <- exact_idm$id
  time_to_illness <- as.numeric(exact_idm$T_ill)
  time_to_death <- as.numeric(exact_idm$T_death)

  eof <- 4 * 365 #+ rnorm(n, 0, sd = sqrt(100))

  T_obs <- pmin(time_to_death, eof)
  status_dead <- T_obs == time_to_death

  last_visit <- numeric(n)

  illness_this_visit <- logical(n)
  illness_last_visit <- logical(n)
  status_illness <- logical(n)
  V_healthy <- numeric(n)
  V_ill <- rep(NA_real_, n)
  is_smaller <- T
  visit_number <- 0
  assessment_stopped <- logical(n)
  is_dead_or_censored <- logical(n)

  unknown_prop <- if (scenario %in% c(1, 2)) 0.0025 else if (scenario %in% c(3, 4)) 0.015

  status_illness_missing <- logical(n)

  while (is_smaller) {
    visit_number <- visit_number + 1
    #   Scenarios 1 & 3: Assessments every 6 months ± N(0,100) days
    #   Scenarios 2 & 4: was assessed every 2 months for the first 6 months, every 6 months through 2 years,
    this_visit <-
      if (scenario %in% c(1, 3)) {
        last_visit + pmax(6 * 30.4 + rnorm(n, 0, sqrt(100)), 1)
      } else if (scenario %in% c(2, 4)) {
        last_visit + pmax(ifelse(last_visit <= 6 * 30.4, 2 * 30.4, 6 * 30.4) + ifelse(last_visit >= 0, rnorm(n, 0, 100), 0), 1)
      }

    # Increasing two fold / doubling
    unknown_prop <- unknown_prop * 2


    # Assessment stop migth need to be change to stop_assement_this_visit
    stop_assessment_this_visit <- (runif(n) <= unknown_prop)
    assessment_stopped <- stop_assessment_this_visit | assessment_stopped



    is_dead_or_censored_this_visit <- last_visit < T_obs  & T_obs <= this_visit
    is_dead_or_censored <- is_dead_or_censored_this_visit | is_dead_or_censored

    # If % are wrong we might need to change how death is handled
    # Vectorized with ifelse
    mask <- !stop_assessment_this_visit & !is_dead_or_censored
    illness_this_visit <- mask & (last_visit < time_to_illness) & (time_to_illness <= this_visit)

    status_illness <- status_illness | illness_this_visit

    V_healthy <- ifelse(illness_this_visit, last_visit, V_healthy)
    V_ill     <- ifelse(illness_this_visit, this_visit,  V_ill)


    V_healthy <- ifelse(!is_dead_or_censored & !status_illness & !assessment_stopped, this_visit, V_healthy)

    V_healthy <- ifelse(is_dead_or_censored_this_visit &
                          !status_illness &
                          !stop_assessment_this_visit, T_obs, V_healthy)


    last_visit <- this_visit

    is_smaller <- min(this_visit) <= max(T_obs)
  }


  # Standardized status classification
  has_interval <- !is.na(V_ill)

  status <- character(n)
  status[has_interval & status_dead] <- "A"
  status[has_interval & !status_dead] <- "B"
  status[!has_interval & status_dead & abs(V_healthy - T_obs) < 1e-8] <- "C"
  status[!has_interval & !status_dead & abs(V_healthy - T_obs) < 1e-8] <- "D"
  status[!has_interval & status_dead & (T_obs - V_healthy) >= 1e-8] <- "E"
  status[!has_interval & !status_dead & (T_obs - V_healthy) >= 1e-8] <- "F"

  # Return
  df_obs_idm <- data.frame(
    id = id,
    V_healthy = V_healthy,
    V_ill = V_ill,
    T_obs = T_obs,
    status_dead = as.numeric(status_dead),
    status_ill = as.numeric(has_interval),
    case = factor(status, levels = c("A", "B", "C", "D", "E", "F"))
  )

  list(obs = df_obs_idm, cens_mechanism = NULL)
}

add_censoring_joly <- function(exact_idm) {
  stopifnot("exact_idm" %in% class(exact_idm))

  n <- nrow(exact_idm)

  id <- exact_idm$id
  time_to_illness <- as.numeric(exact_idm$T_ill)

  time_to_death <- as.numeric(exact_idm$T_death)

  eof <- 2 + 50 * runif(n)

  T_obs <- pmin(time_to_death, eof)
  status_dead <- T_obs == time_to_death

  last_visit <- numeric(n)

  illness_this_visit <- logical(n)
  illness_last_visit <- logical(n)
  status_illness <- logical(n)
  V_healthy <- numeric(n)
  V_ill <- rep(NA_real_, n)
  is_smaller <- T
  while (is_smaller) {
    this_visit <- last_visit + 1 + 3 * runif(n)


    is_dead <- this_visit > T_obs

    illness_this_visit <- time_to_illness < this_visit & !is_dead

    illness_just_observed <- !illness_last_visit & illness_this_visit

    status_illness <- status_illness | illness_just_observed

    V_healthy <- ifelse(status_illness | is_dead, V_healthy, this_visit)
    V_ill <- ifelse(illness_just_observed & !is_dead, this_visit, V_ill)

    illness_last_visit <- illness_this_visit
    last_visit <- this_visit

    is_smaller <- min(this_visit) <= max(T_obs)
  }


  # Standardized status classification
  has_interval <- !is.na(V_ill)

  status <- character(n)
  # Case A: illness and death observed
  # Case B: illness observed and censored while ill
  # Case C: death observed at state 1 without illness V_healthy == T_obs
  # Case D: censored at state 1 without illness V_healthy == T_obs
  # Case E: not observed and death observed either at state 1 or 2
  # Case F: not observed and censored either at state 1 or 2

  status[has_interval & status_dead] <- "A"
  status[has_interval & !status_dead] <- "B"
  status[!has_interval & status_dead & abs(V_healthy - T_obs) < 1e-8] <- "C"
  status[!has_interval & !status_dead & abs(V_healthy - T_obs) < 1e-8] <- "D"
  status[!has_interval & status_dead & (T_obs - V_healthy) >= 1e-8] <- "E"
  status[!has_interval & !status_dead & (T_obs - V_healthy) >= 1e-8] <- "F"

  # Return
  df_obs_idm <- data.frame(
    id = id,
    V_healthy = V_healthy,
    V_ill = V_ill,
    T_obs = T_obs,
    status_dead = as.numeric(status_dead),
    status_ill = as.numeric(has_interval),
    case = factor(status, levels = c("A", "B", "C", "D", "E", "F"))
  )

  list(obs = df_obs_idm, cens_mechanism = NULL)
}


#### General Wrapper ####

#' Simulate Illness-Death Model Data with Censoring
#'
#' General wrapper function that applies a censoring mechanism to exact
#' illness-death data and computes summary statistics.
#'
#' @param exact_idm Data frame with class \code{"exact_idm"} containing true event times.
#' @param censoring_fn Function that applies censoring to \code{exact_idm}.
#'   Must return a list with components \code{obs} (observed data) and
#'   \code{cens_mechanism} (censoring details).
#' @param ... Additional arguments passed to \code{censoring_fn}.
#'
#' @return A list with class \code{"simulated_idm"} containing:
#'   \describe{
#'     \item{obs}{Observed (censored) data frame}
#'     \item{cens_mechanism}{Censoring mechanism details}
#'     \item{summary}{Summary statistics from \code{\link{summarise_simulated_data}}}
#'   }
#'
#' @seealso \code{\link{simulate_exact_idm}}, \code{\link{add_censoring_type_1}}
#' @keywords internal
simulate_idm <- function(exact_idm, censoring_fn, ...) {
  # Apply censoring function
  censored_data <- censoring_fn(exact_idm, ...)

  # Extract observed and censoring mechanism data
  obs_data <- censored_data$obs
  cens_data <- censored_data$cens_mechanism

  # Compute summary statistics
  summary_stats <- summarise_simulated_data(obs_data, exact_idm, cens_data)

  # Return structured result
  result <- list(
    obs = obs_data,
    cens_mechanism = cens_data,
    summary = summary_stats
  )

  class(result) <- c("simulated_idm", class(result))
  result
}

#### Specific Simulation Endpoints ####

#' Simulate Illness-Death Data with Constant Hazards
#'
#' Generates illness-death model data with constant transition hazards and
#' interval-censored illness times. Returns a complete \code{idm_object} with
#' data and true estimators.
#'
#' @param n Integer. Number of subjects to simulate. Default is 300.
#' @param a12 Numeric. Constant hazard for 1→2 (healthy to illness) transition. Default is 0.0008.
#' @param a13 Numeric. Constant hazard for 1→3 (healthy to death) transition. Default is 0.0002.
#' @param a23 Numeric. Constant hazard for 2→3 (illness to death) transition. Default is 0.0016.
#' @param average_number_of_visits Numeric. Average number of observation visits. Default is 10.
#'
#' @return An object of class \code{c("idm_object", "idm_exact_object")} containing:
#'   \describe{
#'     \item{data}{Simulated observed data frame}
#'     \item{model_type}{Character string describing the simulation model}
#'     \item{model_config}{List with simulation parameters and summaries}
#'     \item{estimators}{True hazard, cumulative hazard, and distribution functions}
#'     \item{exact_idm}{Data frame with true (uncensored) event times}
#'   }
#'
#' @details
#' These default hazard values are commonly used in illness-death model literature
#' and provide realistic event rates for many applications.
#'
#' @examples
#' set.seed(789)
#' sim_data <- simulate_idm_constant_hazards(n = 200, average_number_of_visits = 8)
#' head(sim_data$data)
#'
#' # View summary comparing observed vs exact data
#' summary(sim_data)
#'
#' # Compare with estimated model
#' fit <- fit_pc_model(sim_data$data, n_knots = 5)
#' plot(fit, sim_data)
#'
#' @seealso \code{\link{simulate_idm_weibull}}, \code{\link{fit_pc_model}}
#' @export
simulate_idm_constant_hazards <- function(
  n = 300,
  a12 = 0.0008,
  a13 = 0.0002,
  a23 = 0.0016,
  mean_time_between_visits = 200,
  sd_time_between_visits = 20,
  max_right_censoring_time = 10000,
  prob_censoring_at_last_visit = 0.2
) {

  estimators <- create_constant_hazard(a12, a13, a23)

  exact_idm <- simulate_exact_idm(
    n = n,
    a12 = estimators$hazard_functions$a12,
    a13 = estimators$hazard_functions$a13,
    a23 = estimators$hazard_functions$a23
  )

  # Apply censoring and summarize
  result <- simulate_idm(
    exact_idm = exact_idm,
    censoring_fn = add_censoring_type_1,
    mean_time_between_visits,
    sd_time_between_visits,
    max_right_censoring_time,
    prob_censoring_at_last_visit
  )

  data <- result$obs

  model_config <- list(
    censoring_mechanism = "constant_hazards",
    observation_scheme = result$cens_mechanism,
    censoring_summary = result$summary,
    a12 = a12,
    a13 = a13,
    a23 = a23,
    mean_time_between_visits = mean_time_between_visits,
    sd_time_between_visits = sd_time_between_visits,
    max_right_censoring_time= max_right_censoring_time,
    prob_censoring_at_last_visit = prob_censoring_at_last_visit
  )

  retval <- list(
    data = data,
    model_type = paste0("simulation_constant_hazards_n", n),
    model_config = model_config,
    estimators = estimators,
    exact_idm = exact_idm
  )

  class(retval) <- c("idm_object", "idm_exact_object", class(retval))
  retval
}



#' Create Weibull hazards
#
#
#' @param shape12 Numeric. Shape parameter for Weibull hazard from state 1 to 2.
#' @param scale12 Numeric. Scale parameter for Weibull hazard from state 1 to 2.
#' @param shape13 Numeric. Shape parameter for Weibull hazard from state 1 to 3.
#' @param scale13 Numeric. Scale parameter for Weibull hazard from state 1 to 3.
#' @param shape23 Numeric. Shape parameter for Weibull hazard from state 2 to 3.
#' @param scale23 Numeric. Scale parameter for Weibull hazard from state 2 to 3.
#
#' @return idm_estimators object.
#'
#' @export

create_weibull_hazard <- function(shape12 = 1.5, scale12 = 1,
                                  shape13 = 1.5, scale13 = 1,
                                  shape23 = 1.5, scale23 = 1,
                                  t_max  = NULL,
                                  n_grid = 3000L) {
  # Weibull hazard: (k/λ) * (t/λ)^(k-1), vectorized and safe at t=0
  h_weibull <- function(t, shape, scale) {
    t <- pmax(as.numeric(t), .Machine$double.eps)
    (shape / scale) * (t / scale)^(shape - 1)
  }

  a12 <- function(t) h_weibull(t, shape12, scale12)
  a13 <- function(t) h_weibull(t, shape13, scale13)
  a23 <- function(t) h_weibull(t, shape23, scale23)

  A12 <- function(t) (t / scale12)^shape12
  A13 <- function(t) (t / scale13)^shape13
  A23 <- function(t) (t / scale23)^shape23

  ## ---------- Fast F12 and F13 via precomputed cumulative integrals ----------

  # If t_max not supplied, set it to a high quantile of the Weibull distributions
  if (is.null(t_max)) {
    t_max <- max(
      stats::qweibull(0.999, shape = shape12, scale = scale12),
      stats::qweibull(0.999, shape = shape13, scale = scale13),
      stats::qweibull(0.999, shape = shape23, scale = scale23)
    )
  }

  # Time grid on [0, t_max]
  t_grid <- seq(0, t_max, length.out = n_grid)

  # Integrands on the grid
  integrand12 <- exp(-A12(t_grid) - A13(t_grid)) * a12(t_grid)
  integrand13 <- exp(-A12(t_grid) - A13(t_grid)) * a13(t_grid)

  # Cumulative trapezoidal integration (corrected)
  dt     <- diff(t_grid)
  f12_lo <- head(integrand12, -1)
  f12_hi <- tail(integrand12, -1)
  f13_lo <- head(integrand13, -1)
  f13_hi <- tail(integrand13, -1)

  trap12 <- dt * (f12_lo + f12_hi) / 2
  trap13 <- dt * (f13_lo + f13_hi) / 2

  F12_vals <- c(0, cumsum(trap12))
  F13_vals <- c(0, cumsum(trap13))

  # Standalone, vectorized approximations of the integrals
  F12 <- stats::approxfun(t_grid, F12_vals, rule = 2)  # constant extrapolation
  F13 <- stats::approxfun(t_grid, F13_vals, rule = 2)

  ## ------------------------------ P22 as before ------------------------------

  P22 <- function(t, entry_time = 0) {
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

  class(estimators) <- c("idm_estimators", class(estimators))
  estimators
}
#' Create constant hazards
#'
#' Creates constant hazard functions for the illness-death model.
#'
#' @param a12_const Numeric. Constant hazard for 1→2 (healthy to illness) transition.
#' @param a13_const Numeric. Constant hazard for 1→3 (healthy to death) transition.
#' @param a23_const Numeric. Constant hazard for 2→3 (illness to death) transition.
#'
#' @return idm_estimators object.
#'
#' @export

create_constant_hazard <- function(a12, a13, a23) {

  estimators <- list(
    hazard_functions = list(
      a12 = function(t) rep(a12,length(t)),
      a13 = function(t) rep(a13,length(t)),
      a23 = function(t) rep(a23,length(t))
    ),
    cum_hazard_functions = list(
      A12 = function(t) a12 * t,
      A13 = function(t) a13 * t,
      A23 = function(t) a23 * t
    ),
    distribution_functions = list(
      F12 = function(t)
        a12/(a12+a13)*(1 - exp(-(a12 + a13) * t)),
      F13 = function(t)
        a13/(a12+a13)*(1 - exp(-(a12 + a13) * t)),
      P22 = function(t, entry_time = 0) ifelse(t >= entry_time, exp(-a23 * (t - entry_time)), NA_real_)
    )
  )
  class(estimators) <- c("idm_estimators", class(estimators))
  estimators
}



#' Simulate Illness-Death Data with using Joly's approach
#'
#'
#' @param n Integer. Number of subjects to simulate. Default is 1000.
#
#' @return An object of class \code{c("idm_object", "idm_exact_object")} containing:
#'   \describe{
#'     \item{data}{Simulated observed data frame}
#'     \item{model_type}{Character string describing the simulation model}
#'     \item{model_config}{List with simulation parameters and summaries}
#'     \item{estimators}{True hazard, cumulative hazard, and distribution functions}
#'     \item{exact_idm}{Data frame with true (uncensored) event times}
#'   }
#'
#' @export
simulate_idm_joly <- function(n) {
  # Mixture hazard: 0.4*Gamma(37,1.5) + 0.6*Gamma(20,2)
  h_mix <- function(t) {
    f1 <- dgamma(t, shape = 37, rate = 1.5)
    s1 <- pgamma(t, shape = 37, rate = 1.5, lower.tail = FALSE)
    f2 <- dgamma(t, shape = 20, rate = 2)
    s2 <- pgamma(t, shape = 20, rate = 2, lower.tail = FALSE)
    (0.4 * f1 + 0.6 * f2) / (0.4 * s1 + 0.6 * s2)
  }

  h_weibull <- function(t, shape, scale) {
    f <- dweibull(t, shape = shape, scale = scale)
    S <- pweibull(t, shape = shape, scale = scale, lower.tail = FALSE)
    f / S
  }

  shape23 <- 2.5
  scale23 <- 1 / 0.08
  shape13 <- 3
  scale13 <- 1 / 0.04

  a12 <- h_mix
  a13 <- function(t) h_weibull(t, shape13, scale13)
  a23 <- function(t) h_weibull(t, shape23, scale23)

  # Simulate exact data
  exact_idm <- simulate_exact_idm(n = n, a12 = a12, a13 = a13, a23 = a23)

  # Apply Joly censoring and summarize
  result <- simulate_idm(
    exact_idm = exact_idm,
    censoring_fn = add_censoring_joly
  )

  data <- result$obs

  model_config <- list(
    censoring_mechanism = "Joly",
    observation_scheme = result$cens_mechanism,
    censoring_summary = result$summary,
    shape13 = shape13,
    scale13 = scale13,
    shape23 = shape23,
    scale23 = scale23
  )

  A12 = function(t) {
    sapply(t, function(x) {
      stats::integrate(a12, lower = 0, upper = x,
                       stop.on.error = TRUE, abs.tol = 1e-8)$value
    })
  }
  A13 = function(t) (t / scale13)^shape13
  A23 = function(t) (t / scale23)^shape23


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
      F12 = Vectorize(function(t){
        stats::integrate(
          function(s) exp(-A12(s)-A13(s))*a12(s),
          lower = 0, upper = t)$value
      }),
      F13 = Vectorize(function(t){
        stats::integrate(
          function(s) exp(-A12(s)-A13(s))*a13(s),
          lower = 0, upper = t)$value
      }),
      P22 = function(t, entry_time = 0) {
        ifelse(t >= entry_time, exp(-A23(t) + A23(entry_time)), NA_real_)
      }
    )
  )
  class(estimators) <- c("idm_estimators", class(estimators))

  retval <- list(
    data = data,
    model_type = paste0("simulation_joly_n", n),
    model_config = model_config,
    estimators = estimators,
    exact_idm = exact_idm
  )

  class(retval) <- c("idm_object", "idm_exact_object", class(retval))

  retval
}



#' Simulate Illness-Death Data with using Frydman and Szareks's approach
#'
#'
#' @param n Integer. Number of subjects to simulate. Default is 1000.
#' @param scenario Determining censoring scheme default is 1.
#' @return An object of class \code{c("idm_object", "idm_exact_object")} containing:
#'   \describe{
#'     \item{data}{Simulated observed data frame}
#'     \item{model_type}{Character string describing the simulation model}
#'     \item{model_config}{List with simulation parameters and summaries}
#'     \item{estimators}{True hazard, cumulative hazard, and distribution functions}
#'     \item{exact_idm}{Data frame with true (uncensored) event times}
#'   }
#'
#' @export
simulate_idm_frydman <- function(n, scenario = 1L) {
  a12_c <- 0.0008
  a13_c <- 0.0016
  a23_c <- 0.0002


  estimators <- create_constant_hazard(a12_c,a13_c,a23_c)

  # Simulate exact data
  exact_idm <- simulate_exact_idm(
    n = n,
    a12 = estimators$hazard_functions$a12,
    a13 = estimators$hazard_functions$a13,
    a23 = estimators$hazard_functions$a23
  )

  # Apply Frydman censoring and summarize
  result <- simulate_idm(
    exact_idm = exact_idm,
    censoring_fn = add_censoring_frydman,
    scenario = scenario
  )

  data <- result$obs

  model_config <- list(
    censoring_mechanism = "Frydman",
    observation_scheme = result$cens_mechanism,
    censoring_summary = result$summary,
    scenario = scenario,
    a12 = a12_c,
    a13 = a13_c,
    a23 = a23_c
  )



  retval <- list(
    data = data,
    model_type = paste0("simulation_frydman_n", n, "_scenario", scenario),
    model_config = model_config,
    estimators = estimators,
    exact_idm = exact_idm
  )

  class(retval) <- c("idm_object", "idm_exact_object", class(retval))
  retval
}
