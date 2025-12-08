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
#' @seealso \code{\link{add_censoring}}, \code{\link{simulate_idm}}
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
#' Applies an observation schedule with interval censoring and right-censoring
#' to exact illness-death model data. Illness times become interval-censored between
#' visits, and subjects may be censored before death.
#'
#' @param exact_idm Data frame with class \code{"exact_idm"} from
#'   \code{\link{simulate_exact_idm}}, containing true event times.
#' @param average_number_of_visits Numeric. Average number of observation visits
#'   per subject across the follow-up period.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{obs}{Data frame with observed (censored) data containing:
#'       \itemize{
#'         \item \code{id}: Subject ID
#'         \item \code{V_0}: Baseline visit time (always 0)
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
#' @examples
#' set.seed(456)
#' exact_data <- simulate_exact_idm(
#'   n = 50,
#'   a12 = function(t) rep(0.01, length(t)),
#'   a13 = function(t) rep(0.005, length(t)),
#'   a23 = function(t) rep(0.02, length(t))
#' )
#' censored_data <- add_censoring(exact_data, average_number_of_visits = 10)
#' head(censored_data$obs)
#' table(censored_data$obs$case)
#'
#' @seealso \code{\link{simulate_exact_idm}}, \code{\link{simulate_idm}}
#'
#' @export
add_censoring <- function(exact_idm,
                          average_number_of_visits) {

  stopifnot("exact_idm" %in% class(exact_idm))
  obs_time_sd <- 0.1
  # Extract core times
  id <- exact_idm$id
  time_to_illness <- as.numeric(exact_idm$T_ill)
  time_to_death <- as.numeric(exact_idm$T_death)
  n <- length(time_to_illness)

  # ---- Censoring: subject-specific Uniform(0, 3 * death_i) with Inf-safe fallback
  time_to_censor <- pmin(1 + runif(n) * 3 * time_to_death)


  # Observation cutoff per subject
  T_cutoff <- pmin(time_to_death, time_to_censor)

  # ---- Build an observation schedule up to the global max cutoff
  max_follow_up <- max(T_cutoff[is.finite(T_cutoff)], 0)

  obs_interval <- max_follow_up / average_number_of_visits
  grid <- seq(0, max_follow_up, by = obs_interval)

  obs_schedule <- matrix(rep(grid, n), nrow = n, byrow = TRUE)
  average_number_of_visits <- ncol(obs_schedule)

  # Jitter the schedule but preserve monotonicity
  if (average_number_of_visits > 1 && obs_time_sd > 0) {
    noise <- matrix(rnorm(n * average_number_of_visits, mean = 0, sd = obs_time_sd), nrow = n)
    clip <- max(min(obs_interval / 2 - 1e-6, 3 * obs_time_sd), 0) # keep order
    if (clip > 0) {
      noise <- pmin(pmax(noise, -clip), clip)
      obs_schedule <- pmax(obs_schedule + noise, 0)
      # enforce increasing times row-wise
      obs_schedule <- t(apply(obs_schedule, 1, sort))
    }
  }
  # Force first visit at time 0
  obs_schedule[, 1] <- abs(obs_schedule[, 1])
  V_0 <- rep(0, n)

  # ---- Determine last healthy visit before illness (or before cutoff if earlier)
  # We use T_end = min(illness time, cutoff) row-wise
  T_end <- pmin(time_to_illness, T_cutoff)
  # Count visits strictly before T_end, row-wise
  idx_healthy <- rowSums(sweep(obs_schedule, 1, T_end, "<"))

  # V_healthy: if no visit before T_end, use baseline 0
  V_healthy <- V_0
  sel <- idx_healthy > 0
  V_healthy[sel] <- obs_schedule[cbind(which(sel), idx_healthy[sel])]

  # ---- First visit after the last healthy visit (candidate "ill" visit)
  idx_ill <- pmin(idx_healthy + 1L, average_number_of_visits)
  V_ill_candidate <- obs_schedule[cbind(seq_len(n), idx_ill)]

  # Illness occurs before cutoff?
  ill_before_cutoff <- is.finite(time_to_illness) & (time_to_illness <= T_cutoff)

  # Candidate ill-visit is valid only if it happens on/before cutoff and there is a next visit
  has_next_visit_before_cutoff <- ill_before_cutoff &
    (idx_healthy < average_number_of_visits) &
    (V_ill_candidate <= T_cutoff + 1e-12)

  V_ill <- ifelse(has_next_visit_before_cutoff, V_ill_candidate, NA_real_)

  # ---- Status classification (purely time-based; no reliance on 'path' labels)
  # 1: no illness observed (no V_ill), alive at cutoff
  # 2: no illness observed (no V_ill), died at cutoff
  # 3: interval-censored illness (have V_ill), alive at cutoff
  # 4: interval-censored illness (have V_ill), died at cutoff
  died_at_cutoff <- is.finite(time_to_death) & (time_to_death <= time_to_censor)
  has_interval <- !is.na(V_ill)




  status <- character(n)
  # Case A illness and death observed
  # Case B illness observed and censored while ill
  # Case C death observed at state 1 without illness V_healthy == T_obs
  # Case D censored at state 1 without illness V_healthy == T_obs
  # Case E not observed and death observed either at state 1 or 2
  # Case F not observed and censored either at state 1 or 2

  status[has_interval & died_at_cutoff] <- "A"
  status[has_interval & !died_at_cutoff] <- "B"
  status[!has_interval & died_at_cutoff & abs(V_healthy - T_cutoff) < 1e-8] <- "C"
  status[!has_interval & !died_at_cutoff & abs(V_healthy - T_cutoff) < 1e-8] <- "D"
  status[!has_interval & died_at_cutoff & (T_cutoff- V_healthy) >= 1e-8] <- "E"
  status[!has_interval & !died_at_cutoff & (T_cutoff - V_healthy) >= 1e-8] <- "F"

  # For status 1 (alive, never observed ill), the illness time is right-censored at the last healthy visit
  T_obs <- T_cutoff
  #  T_obs[status == 1] <- V_healthy[status == 1]

  # Return
  df_obs_idm <- data.frame(
    id = id,
    V_0 = V_0,
    V_healthy = V_healthy,
    V_ill = V_ill,
    T_obs = T_obs,
    status_dead = as.numeric(died_at_cutoff),
    status_ill = as.numeric(has_interval),
    case = factor(
      status,
      levels = c("A", "B", "C", "D", "E", "F")
      # labels = c(
      #   "illness observed, died@cutoff",
      #   "illness observed, alive@cutoff",
      #   "died@cutoff (healty before death)",
      #   "healthy@cutoff (healty before censoring)",
      #   "died@cutoff (healty or ill before death)",
      #   "alive@cutoff (healty or ill before censoring)"
      # )
    )
  )

  # Store the (possibly wide) schedule as a list-column to avoid unintended column expansion
  df_observation_scheme <- data.frame(
    id = id,
    obs_schedule = obs_schedule,
    time_to_censor = time_to_censor
  )

  list(obs = df_obs_idm, cens_mechanism = df_observation_scheme)
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
    V_0 = 0,
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
    V_0 = 0,
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
#' @seealso \code{\link{simulate_exact_idm}}, \code{\link{add_censoring}}
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
  average_number_of_visits = 10
) {
  # Create constant hazard functions
  a12_const <- function(t) rep(a12, length(t))
  a13_const <- function(t) rep(a13, length(t))
  a23_const <- function(s) rep(a23, length(s))

  # Simulate exact data
  exact_idm <- simulate_exact_idm(
    n = n,
    a12 = a12_const,
    a13 = a13_const,
    a23 = a23_const
  )

  # Apply censoring and summarize
  result <- simulate_idm(
    exact_idm = exact_idm,
    censoring_fn = add_censoring,
    average_number_of_visits = average_number_of_visits
  )

  data <- result$obs

  model_config <- list(
    censoring_mechanism = "constant_hazards",
    observation_scheme = result$cens_mechanism,
    censoring_summary = result$summary,
    a12 = a12,
    a13 = a13,
    a23 = a23,
    average_number_of_visits = average_number_of_visits
  )

  estimators <- list(
    hazard_functions = list(
      a12 = a12_const,
      a13 = a13_const,
      a23 = a23_const
    ),
    cum_hazard_functions = list(
      A12 = function(t) a12 * t,
      A13 = function(t) a13 * t,
      A23 = function(s) a23 * s
    ),
    distribution_functions = list(
      F12 = function(t) 1 - exp(-a12 * t),
      F13 = function(t) 1 - exp(-a13 * t),
      P22 = function(t, entry_time = 0) ifelse(t >= entry_time, exp(-a23 * (t - entry_time)), NA_real_)
    )
  )
  class(estimators) <- c("idm_estimators", class(estimators))

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

#' Simulate Illness-Death Data with Weibull Hazards
#'
#' Generates illness-death model data with Weibull-distributed transition times.
#' Allows flexible specification of time-dependent hazards through shape and scale
#' parameters.
#'
#' @param n Integer. Number of subjects to simulate. Default is 1000.
#' @param shape12 Numeric. Weibull shape parameter for 1→2 transition. Default is 3.
#' @param scale12 Numeric. Weibull scale parameter for 1→2 transition. Default is 1.
#' @param shape13 Numeric. Weibull shape parameter for 1→3 transition. Default is 5.
#' @param scale13 Numeric. Weibull scale parameter for 1→3 transition. Default is 1.
#' @param shape23 Numeric. Weibull shape parameter for 2→3 transition. Default is 2.
#' @param scale23 Numeric. Weibull scale parameter for 2→3 transition. Default is 1.
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
#' The Weibull hazard function is \eqn{h(t) = (k/\lambda)(t/\lambda)^{k-1}}, where
#' \code{k} is the shape and \code{\lambda} is the scale parameter.
#'
#' Shape parameters control the hazard's time dependence:
#' \itemize{
#'   \item k < 1: Decreasing hazard over time
#'   \item k = 1: Constant hazard (exponential distribution)
#'   \item k > 1: Increasing hazard over time
#' }
#'
#' @examples
#' set.seed(101)
#' sim_data <- simulate_idm_weibull(
#'   n = 500,
#'   shape12 = 2, scale12 = 10,
#'   shape13 = 3, scale13 = 15,
#'   shape23 = 1.5, scale23 = 8
#' )
#' # View comprehensive summary with exact vs observed comparison
#' summary(sim_data)
#'
#' @seealso \code{\link{simulate_idm_constant_hazards}}, \code{\link{fit_spline_model}}
#' @export
simulate_idm_weibull <- function(
  n = 1000,
  shape12 = 3, scale12 = 1,
  shape13 = 5, scale13 = 1,
  shape23 = 2, scale23 = 1,
  average_number_of_visits = 10
) {
  # Validate inputs
  params <- c(shape12, scale12, shape13, scale13, shape23, scale23)
  if (any(!is.finite(params)) || any(params <= 0)) {
    stop("All shapes and scales must be positive and finite.")
  }

  # Weibull hazard: (k/λ) * (t/λ)^(k-1), vectorized and safe at t=0
  h_weibull <- function(t, shape, scale) {
    t <- pmax(as.numeric(t), .Machine$double.eps)
    (shape / scale) * (t / scale)^(shape - 1)
  }

  a12 <- function(t) h_weibull(t, shape12, scale12)
  a13 <- function(t) h_weibull(t, shape13, scale13)
  a23 <- function(t) h_weibull(t, shape23, scale23)

  # Simulate exact data
  exact_idm <- simulate_exact_idm(
    n = n,
    a12 = a12,
    a13 = a13,
    a23 = a23
  )

  # Apply censoring and summarize
  result <- simulate_idm(
    exact_idm = exact_idm,
    censoring_fn = add_censoring,
    average_number_of_visits = average_number_of_visits
  )

  data <- result$obs

  model_config <- list(
    censoring_mechanism = "weibull",
    observation_scheme = result$cens_mechanism,
    censoring_summary = result$summary,
    shape12 = shape12, scale12 = scale12,
    shape13 = shape13, scale13 = scale13,
    shape23 = shape23, scale23 = scale23,
    average_number_of_visits = average_number_of_visits
  )

  estimators <- list(
    hazard_functions = list(
      a12 = a12,
      a13 = a13,
      a23 = a23
    ),
    cum_hazard_functions = list(
      A12 = function(t) (t / scale12)^shape12,
      A13 = function(t) (t / scale13)^shape13,
      A23 = function(t) (t / scale23)^shape23
    ),
    distribution_functions = list(
      F12 = function(t) 1 - exp(-(t / scale12)^shape12),
      F13 = function(t) 1 - exp(-(t / scale13)^shape13),
      P22 = function(t, entry_time = 0) ifelse(t >= entry_time, exp(-(t / scale23)^shape23 + (entry_time / scale23)^shape23), NA_real_)
    )
  )
  class(estimators) <- c("idm_estimators", class(estimators))

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


  estimators <- list(
    hazard_functions = list(
      a12 = a12,
      a13 = a13,
      a23 = a23
    ),
    cum_hazard_functions = list(
      A12 = function(t) {
        sapply(t, function(x) {
          stats::integrate(a12, lower = 0, upper = x,
                          stop.on.error = TRUE, abs.tol = 1e-8)$value
        })
      },
      A13 = function(t) (t / scale13)^shape13,
      A23 = function(t) (t / scale23)^shape23
    ),
    distribution_functions = list(
      F12 = function(t) {
        sapply(t, function(x) {
          1 - exp(-stats::integrate(a12, lower = 0, upper = x,
                                    stop.on.error = TRUE, abs.tol = 1e-8)$value)
        })
      },
      F13 = function(t) 1 - exp(-(t / scale13)^shape13),
      P22 = function(t, entry_time = 0) ifelse(t >= entry_time, exp(-(t / scale23)^shape23 + (entry_time / scale23)^shape23), NA_real_)
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
  # Constant hazards as in Frydman paper
  a12 <- function(x) rep(0.0008, length(x))
  a23 <- function(x) rep(0.0016, length(x))
  a13 <- function(x) rep(0.0002, length(x))

  # Simulate exact data
  exact_idm <- simulate_exact_idm(
    n = n,
    a12 = a12,
    a23 = a23,
    a13 = a13
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
    a12 = 0.0008,
    a13 = 0.0002,
    a23 = 0.0016
  )

  estimators <- list(
    hazard_functions = list(
      a12 = a12,
      a13 = a13,
      a23 = a23
    ),
    cum_hazard_functions = list(
      A12 = function(t) 0.0008 * t,
      A13 = function(t) 0.0002 * t,
      A23 = function(s) 0.0016 * s
    ),
    distribution_functions = list(
      F12 = function(t) 1 - exp(-0.0008 * t),
      F13 = function(t) 1 - exp(-0.0002 * t),
      P22 = function(t, entry_time = 0) ifelse(t >= entry_time, exp(-0.0016 * (t - entry_time)), NA_real_)
    )
  )
  class(estimators) <- c("idm_estimators", class(estimators))

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
