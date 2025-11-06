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


add_censoring <- function(exact_idm, target_pct_censoring, average_number_of_visits) {
  stopifnot("exact_idm" %in% class(exact_idm))
  obs_time_sd <- 0.1
  # Extract core times
  id <- exact_idm$id
  time_to_illness <- as.numeric(exact_idm$T_ill)
  time_to_death <- as.numeric(exact_idm$T_death)
  n <- length(time_to_illness)

  # ---- Censoring: subject-specific Uniform(0, 3 * death_i) with Inf-safe fallback
  time_to_censor <- min(1 + runif(n) * 3 * time_to_death)


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
  # V_healthy[sel] <- obs_schedule[cbind(which(sel), idx_healthy[sel])]

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


  status <- integer(n)
  status[!has_interval & !died_at_cutoff] <- 1
  status[!has_interval & died_at_cutoff] <- 2
  status[has_interval & !died_at_cutoff] <- 3
  status[has_interval & died_at_cutoff] <- 4

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
      levels = 1:4,
      labels = c(
        "healthy@cutoff (cens)", "died@cutoff (no illness observed)",
        "interval illness, alive@cutoff", "interval illness, died@cutoff"
      )
    )
  )

  # Store the (possibly wide) schedule as a list-column to avoid unintended column expansion
  df_observation_scheme <- data.frame(
    id = id,
    obs_schedule = obs_schedule,
    time_to_censor = time_to_censor
  )

  list(obs = df_obs_idm, cens_mechanism = df_observation_scheme, exact_idm = exact_idm)
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



  status <- integer(n)
  status[!status_illness & !status_dead] <- 1
  status[!status_illness & status_dead] <- 2
  status[status_illness & !status_dead] <- 3
  status[status_illness & status_dead] <- 4


  # Return
  df_obs_idm <- data.frame(
    id = id,
    V_0 = 0,
    V_healthy = V_healthy,
    V_ill = V_ill,
    T_obs = T_obs,
    status_dead = as.numeric(status_dead),
    status_ill = as.numeric(status_illness),
    case = factor(
      status,
      levels = 1:4,
      labels = c(
        "healthy@cutoff (cens)", "died@cutoff (no illness observed)",
        "interval illness, alive@cutoff", "interval illness, died@cutoff"
      )
    )
  )

  list(
    obs = df_obs_idm,
    exact_idm = exact_idm
  )
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


  status <- integer(n)
  status[!status_illness & !status_dead] <- 1
  status[!status_illness & status_dead] <- 2
  status[status_illness & !status_dead] <- 3
  status[status_illness & status_dead] <- 4


  # Return
  df_obs_idm <- data.frame(
    id = id,
    V_0 = 0,
    V_healthy = V_healthy,
    V_ill = V_ill,
    T_obs = T_obs,
    status_dead = as.numeric(status_dead),
    status_ill = as.numeric(status_illness),
    case = factor(
      status,
      levels = 1:4,
      labels = c(
        "healthy@cutoff (cens)", "died@cutoff (no illness observed)",
        "interval illness, alive@cutoff", "interval illness, died@cutoff"
      )
    )
  )

  list(
    obs = df_obs_idm,
    exact_idm = exact_idm
  )
}


#### Wrappers ####

simulate_idm <- function(n, a12, a13, a23,
                         target_pct_censoring,
                         average_number_of_visits) {
  exact_idm <- simulate_exact_idm(
    n = n,
    a12 = a12,
    a13 = a13,
    a23 = a23,
  )

  simulated_data <- add_censoring(
    exact_idm = exact_idm,
    target_pct_censoring = target_pct_censoring,
    average_number_of_visits = average_number_of_visits
  )

  hazards <- list(
    a12 = a12,
    a13 = a13,
    a23 = a23,
    A12 = function(t) {
      sapply(t, function(x) {
        stats::integrate(a12,
          lower = 0, upper = x,
          stop.on.error = TRUE, abs.tol = 1e-8
        )$value
      })
    },
    A13 = function(t) {
      sapply(t, function(x) {
        stats::integrate(a13,
          lower = 0, upper = x,
          stop.on.error = TRUE, abs.tol = 1e-8
        )$value
      })
    },
    A23 = function(s) {
      sapply(s, function(x) {
        stats::integrate(a23,
          lower = 0, upper = x,
          stop.on.error = TRUE, abs.tol = 1e-8
        )$value
      })
    }
  )

  class(hazards) <- c("idm_hazards", class(hazards))

  res <- list(
    datasets = simulated_data,
    hazards = hazards
  )

  class(res) <- c(class(res), "simulated_idm")
  res
}

simulate_idm_constant_hazards <- function(
  n = 300,
  a12 = 0.0008,
  a13 = 0.0002,
  a23 = 0.0016,
  target_pct_censoring = 0.1,
  average_number_of_visits = 10
) {
  a12_const <- function(t) rep(a12, length(t))
  a13_const <- function(t) rep(a13, length(t))
  a23_const <- function(s) rep(a23, length(s))
  simulate_idm(
    n = n,
    a12 = a12_const,
    a13 = a13_const,
    a23 = a23_const,
    target_pct_censoring = target_pct_censoring,
    average_number_of_visits = average_number_of_visits
  )
}

simulate_idm_weibull <- function(
  n = 1000,
  shape12 = 3, scale12 = 1,
  shape13 = 5, scale13 = 1,
  shape23 = 2, scale23 = 1,
  target_pct_censoring = 0.1,
  average_number_of_visits = 10
) {
  # validate inputs
  params <- c(shape12, scale12, shape13, scale13, shape23, scale23)
  if (any(!is.finite(params)) || any(params <= 0)) {
    stop("All shapes and scales must be positive and finite.")
  }

  # Weibull hazard: (k/λ) * (t/λ)^(k-1), vectorized and safe at t=0
  h_weibull <- function(t, shape, scale) {
    t <- pmax(as.numeric(t), .Machine$double.eps)
    (shape / scale) * (t / scale)^(shape - 1)
  }

  a12 <- function(t) h_weibull(t, shape12, scale12) # 1 -> 2
  a13 <- function(t) h_weibull(t, shape13, scale13) # 1 -> 3
  a23 <- function(t) h_weibull(t, shape23, scale23) # 2 -> 3

  simulate_idm(
    n = n,
    a12 = a12,
    a13 = a13,
    a23 = a23,
    target_pct_censoring = target_pct_censoring,
    average_number_of_visits = average_number_of_visits
  )
}


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
  a13 <- function(t) h_weibull(t, shape13, scale13) # 1 -> 3
  a23 <- function(t) h_weibull(t, shape23, scale23) # 2 -> 3

  res <- simulate_exact_idm(n, a12 = a12, a13 = a13, a23 = a23)

  simulated_data <- add_censoring_joly(exact_idm = res)

  hazards <- list(
    a12 = a12,
    a13 = a13,
    a23 = a23,
    A12 = function(t) {
      sapply(t, function(x) {
        stats::integrate(a12,
          lower = 0, upper = x,
          stop.on.error = TRUE, abs.tol = 1e-8
        )$value
      })
    },
    A13 = function(t) {
      sapply(t, function(x) {
        stats::integrate(a13,
          lower = 0, upper = x,
          stop.on.error = TRUE, abs.tol = 1e-8
        )$value
      })
    },
    A23 = function(s) {
      sapply(s, function(x) {
        stats::integrate(a23,
          lower = 0, upper = x,
          stop.on.error = TRUE, abs.tol = 1e-8
        )$value
      })
    }
  )

  class(hazards) <- c("idm_hazards", class(hazards))

  res <- list(
    datasets = simulated_data,
    hazards = hazards
  )

  class(res) <- c(class(res), "simulated_idm")
  res
}

simulate_idm_frydman <- function(n, scenario) {
  a12 <- function(x) rep(0.0008, length(x))
  a23 <- function(x) rep(0.0016, length(x))
  a13 <- function(x) rep(0.0002, length(x))

  exact_data <- simulate_exact_idm(
    n = n,
    a12 = a12,
    a23 = a23,
    a13 = a13
  )

  simulated_data <- add_censoring_frydman(exact_data, scenario = scenario)

  hazards <- list(
    a12 = a12,
    a13 = a13,
    a23 = a23,
    A12 = function(t) {
      sapply(t, function(x) {
        stats::integrate(a12,
          lower = 0, upper = x,
          stop.on.error = TRUE, abs.tol = 1e-8
        )$value
      })
    },
    A13 = function(t) {
      sapply(t, function(x) {
        stats::integrate(a13,
          lower = 0, upper = x,
          stop.on.error = TRUE, abs.tol = 1e-8
        )$value
      })
    },
    A23 = function(s) {
      sapply(s, function(x) {
        stats::integrate(a23,
          lower = 0, upper = x,
          stop.on.error = TRUE, abs.tol = 1e-8
        )$value
      })
    }
  )

  class(hazards) <- c("idm_hazards", class(hazards))

  res <- list(
    datasets = simulated_data,
    hazards = hazards
  )

  class(res) <- c(class(res), "simulated_idm")
  res
}
