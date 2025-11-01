
simulate_exact_idm <- function(
    n,
    a12, a13, a23,      # all are functions of ABSOLUTE calendar time t
    t0 = 0,
    tmax = Inf,          # if finite: events after tmax are treated as never occurring (Inf)
    init_step = 1,       # initial bracket step for inversion search
    int_abs_tol = 1e-8,  # absolute tolerance for integrate()
    root_tol = 1e-8,     # tolerance for uniroot()
    max_doublings = 60   # safety cap for growing the search bracket
) {

# --- Safeguards ----------
  if (length(formals(a12)) != 1L || length(formals(a13)) != 1L)
    stop("a12 and a13 must be functions of one argument: absolute time t.")
  if (length(formals(a23)) != 1L)
    stop("a23 must be a function of one argument: absolute time t (calendar-time hazard).")

# --- Helper: cumulative hazard from t_start to t_end ---------
  cumhaz <- function(h, t_start, t_end) {
    if (t_end <= t_start) return(0)
    res <- stats::integrate(function(u) h(u), lower = t_start, upper = t_end,
                            stop.on.error = TRUE, abs.tol = int_abs_tol)
    as.numeric(res$value)
  }

  # --- Helper: draw absolute event time by inverting the cumulative hazard ----
  # Given hazard h(t) and starting clock t_start, draw T ≥ t_start with
  #  P(T > t | T ≥ t_start) = exp(-∫_{t_start}^t h(u) du).
  draw_time_from_hazard <- function(h, t_start) {
    target <- stats::rexp(1)  # -log(U) ~ Exp(1)

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
      ch <- ch + cumhaz(h, lb, ub)  # accumulate ∫ on the new segment
      n_doubles <- n_doubles + 1
    }

    # if we couldn’t reach target before tmax, the event never occurs
    if (!is.finite(ch) || (ub >= tmax && ch < target)) return(Inf)

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
    if (val_ub < 0) return(Inf)

    uniroot(F, lower = lb0, upper = ub0, tol = root_tol)$root  # absolute time
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


add_interval_censoring_to_illness <- function(exact_idm, admin_cutoff, n_obs = 10, obs_time_sd = 0.1) {



  stopifnot("exact_idm" %in% class(exact_idm))

  # Extract core times
  id = exact_idm$id
  time_to_illness <- as.numeric(exact_idm$T_ill)
  time_to_death   <- as.numeric(exact_idm$T_death)
  n <- length(time_to_illness)

  # ---- Censoring: subject-specific Uniform(0, 3 * death_i) with Inf-safe fallback
  time_to_censor <- min(1 + runif(n) * 3 * time_to_death, admin_cutoff)


  # Observation cutoff per subject
  T_cutoff <- pmin(time_to_death, time_to_censor)

  # ---- Build an observation schedule up to the global max cutoff
  max_follow_up <- max(T_cutoff[is.finite(T_cutoff)], 0)

  obs_interval <- max_follow_up/n_obs
  grid <- seq(0, max_follow_up, by = obs_interval)

  obs_schedule <- matrix(rep(grid, n), nrow = n, byrow = TRUE)
  n_obs <- ncol(obs_schedule)

  # Jitter the schedule but preserve monotonicity
  if (n_obs > 1 && obs_time_sd > 0) {
    noise <- matrix(rnorm(n * n_obs, mean = 0, sd = obs_time_sd), nrow = n)
    clip <- max(min(obs_interval / 2 - 1e-6, 3 * obs_time_sd), 0)  # keep order
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
  idx_ill <- pmin(idx_healthy + 1L, n_obs)
  V_ill_candidate <- obs_schedule[cbind(seq_len(n), idx_ill)]

  # Illness occurs before cutoff?
  ill_before_cutoff <- is.finite(time_to_illness) & (time_to_illness <= T_cutoff)

  # Candidate ill-visit is valid only if it happens on/before cutoff and there is a next visit
  has_next_visit_before_cutoff <- ill_before_cutoff &
    (idx_healthy < n_obs) &
    (V_ill_candidate <= T_cutoff + 1e-12)

  V_ill <- ifelse(has_next_visit_before_cutoff, V_ill_candidate, NA_real_)

  # ---- Status classification (purely time-based; no reliance on 'path' labels)
  # 1: no illness observed (no V_ill), alive at cutoff
  # 2: no illness observed (no V_ill), died at cutoff
  # 3: interval-censored illness (have V_ill), alive at cutoff
  # 4: interval-censored illness (have V_ill), died at cutoff
  died_at_cutoff <- is.finite(time_to_death) & (time_to_death <= time_to_censor)
  has_interval   <- !is.na(V_ill)



  status <- integer(n)
  status[!has_interval & !died_at_cutoff] <- 1
  status[!has_interval &  died_at_cutoff] <- 2
  status[ has_interval & !died_at_cutoff] <- 3
  status[ has_interval &  died_at_cutoff] <- 4

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
      labels = c("healthy@cutoff (cens)", "died@cutoff (no illness observed)",
                 "interval illness, alive@cutoff", "interval illness, died@cutoff")
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

# add_interval_censoring_to_illness_extension <- function(
#     dt,
#     scenario = 1L,                # 1..4: visit grids & missingness per Frydman–Szarek
#     study_end = 1460,             # days (≈ 4 years)
#     obs_time_sd = 100,            # days; N(0, sd) jitter at each planned visit (except baseline)
#     seed = NULL
# ) {
#   # Required columns
#   stopifnot(all(c("entry2", "death_time") %in% names(dt)))
#   n <- nrow(dt)
#   time_to_illness <- as.numeric(dt$entry2)
#   time_to_death   <- as.numeric(dt$death_time)
#
#   # ---- Scenario-specific visit schedule (in days)
#   # Scen. 1 & 3: every 6 months; Scen. 2 & 4: 0, 2, 4, 6 months; then every 6m through 24m; then yearly
#   six_months <- 182.5
#   base_schedule <- switch(
#     as.character(scenario),
#     "1" = seq(0, study_end, by = six_months),
#     "3" = seq(0, study_end, by = six_months),
#     "2" = {
#       v <- c(0, 60, 120, 180, 360, 540, 720, 1080, 1460)
#       v[v <= study_end]
#     },
#     "4" = {
#       v <- c(0, 60, 120, 180, 360, 540, 720, 1080, 1460)
#       v[v <= study_end]
#     },
#     stop("scenario must be 1, 2, 3, or 4")
#   )
#   m <- length(base_schedule)
#   if (m < 2) stop("Study end too short for any post-baseline visits.")
#
#   # ---- Build per-subject schedule with jitter; keep monotone and within [0, study_end]
#   if (!is.null(seed)) set.seed(seed)
#   obs_schedule <- matrix(rep(base_schedule, each = n), nrow = n, ncol = m)
#   if (obs_time_sd > 0) {
#     noise <- matrix(rnorm(n * (m - 1), mean = 0, sd = obs_time_sd), nrow = n, ncol = (m - 1))
#     obs_schedule[, 2:m] <- obs_schedule[, 2:m] + noise
#     obs_schedule[, 2:m] <- pmax(pmin(obs_schedule[, 2:m], study_end), 0)  # clamp
#     # enforce increasing times row-wise; baseline fixed at 0
#     obs_schedule <- t(apply(obs_schedule, 1, function(x) { x <- sort(x); x[1] <- 0; x }))
#   }
#
#   # ---- Missingness of the intermediate status (illness) by assessment; baseline never missing
#   # p1 = 0.5% in Scen. 1–2; p1 = 3% in Scen. 3–4; doubles each assessment; cap at 0.98.
#   p1 <- if (scenario %in% c(1L, 2L)) 0.005 else 0.03
#   miss_prob <- pmin(0.98, p1 * 2^(0:(m - 2)))     # length m-1 (post-baseline)
#   miss_mat  <- matrix(runif(n * (m - 1)) < rep(miss_prob, each = n), nrow = n, ncol = (m - 1))
#   missing_status <- cbind(rep(FALSE, n), miss_mat)  # include baseline (never missing)
#
#   # ---- Administrative end & death
#   T_end <- pmin(time_to_death, study_end)
#   died_by_end <- is.finite(time_to_death) & (time_to_death <= study_end)
#
#   # ---- Determine interval [V_healthy, V_ill] if illness occurs before end and before death
#   V_0 <- rep(0, n)
#   V_healthy <- rep(0, n)
#   V_ill <- rep(NA_real_, n)
#
#   ill_before_end <- is.finite(time_to_illness) & (time_to_illness <= study_end) & (time_to_illness < time_to_death)
#
#   for (i in seq_len(n)) {
#     t_vis  <- obs_schedule[i, ]
#     known  <- (!missing_status[i, ]) & (t_vis <= (T_end[i] + 1e-12))
#     known[1] <- TRUE  # baseline known
#
#     if (ill_before_end[i]) {
#       idx_before <- which((t_vis < time_to_illness[i]) & known)
#       j_star <- if (length(idx_before)) max(idx_before) else 1L
#       idx_after <- which(seq_len(m) > j_star & known)
#       if (length(idx_after)) {
#         k <- min(idx_after)
#         V_healthy[i] <- t_vis[j_star]
#         V_ill[i]     <- t_vis[k]
#       } else {
#         # No known post-illness visit before T_end
#         V_healthy[i] <- t_vis[j_star]
#         V_ill[i]     <- NA_real_
#       }
#     } else {
#       # No illness before end/death: last known healthy visit before T_end
#       idx_known <- which(known & (t_vis <= (T_end[i] + 1e-12)))
#       V_healthy[i] <- if (length(idx_known)) max(t_vis[idx_known]) else 0
#       V_ill[i]     <- NA_real_
#     }
#   }
#
#   # ---- Status categories (as in your original)
#   status <- integer(n)
#   status[is.na(V_ill) & !died_by_end] <- 1L
#   status[is.na(V_ill) &  died_by_end] <- 2L
#   status[!is.na(V_ill) & !died_by_end] <- 3L
#   status[!is.na(V_ill) &  died_by_end] <- 4L
#
#   # ---- Outputs
#   obs_dt <- data.frame(
#     V_0 = V_0,
#     V_healthy = V_healthy,
#     V_ill = V_ill,
#     T_obs = T_end,
#     status = factor(
#       status, levels = 1:4,
#       labels = c("healthy@admin_end (cens)",
#                  "died@admin_end (no illness observed)",
#                  "interval illness, alive@admin_end",
#                  "interval illness, died@admin_end")
#     )
#   )
#
#   # keep schedules and missingness as list-cols
#   schedule_list <- split(obs_schedule, row(obs_schedule))
#   missing_list  <- split(missing_status, row(missing_status))
#
#   true_dt <- data.frame(
#     scenario = scenario,
#     obs_schedule = schedule_list,
#     missing_status = missing_list,
#     time_to_illness = time_to_illness,
#     time_to_death = time_to_death,
#     study_end = study_end,
#     died_by_end = died_by_end,
#     status = obs_dt$status
#   )
#
#   structure(
#     list(obs = obs_dt, true = true_dt),
#     scenario = scenario,
#     missing_p1 = p1,
#     obs_time_sd = obs_time_sd,
#     class = "fs_interval_obs"
#   )
# }

simulate_idm <- function(n, a12, a13, a23, admin_cutoff = Inf) {
  res <- simulate_exact_idm(
    n,
    a12 = a12,
    a13 = a13,
    a23 = a23,
  )

  simulated_data <- add_interval_censoring_to_illness(res, admin_cutoff)

  hazards <- list(
    a12 = a12,
    a13 = a13,
    a23 = a23,
    A12 = function(t) {
      sapply(t, function(x) {
        stats::integrate(a12, lower = 0, upper = x,
                         stop.on.error = TRUE, abs.tol = 1e-8)$value
      })
    },
    A13 = function(t) {
      sapply(t, function(x) {
        stats::integrate(a13, lower = 0, upper = x,
                         stop.on.error = TRUE, abs.tol = 1e-8)$value
      })
    },
    A23 = function(s) {
      sapply(s, function(x) {
        stats::integrate(a23, lower = 0, upper = x,
                         stop.on.error = TRUE, abs.tol = 1e-8)$value
      })
    }
  )

  class(hazards) <- c("idm_hazards", class(hazards))

  res <- list(datasets = simulated_data,
       hazards = hazards)

  class(res)  <- c(class(res), "simulated_idm")
  res
}

simulate_idm_constant_hazards <- function(
    n = 1000,
    a12 = 0.1,
    a13 = 0.2,
    a23 = 0.3) {

  a12_const <- function(t) rep(a12, length(t))
  a13_const <- function(t) rep(a13, length(t))
  a23_const <- function(s) rep(a23, length(s))
  simulate_idm(n, a12_const, a13_const, a23_const)
}

simulate_idm_weibull <- function(
    n = 1000,
    shape12 = 3, scale12 = 1,
    shape13 = 5, scale13 = 1,
    shape23 = 2, scale23 = 1,
    admin_cutoff = 2.5) {
  # validate inputs
  params <- c(shape12, scale12, shape13, scale13, shape23, scale23)
  if (any(!is.finite(params)) || any(params <= 0))
    stop("All shapes and scales must be positive and finite.")

  # Weibull hazard: (k/λ) * (t/λ)^(k-1), vectorized and safe at t=0
  h_weibull <- function(t, shape, scale) {
    t <- pmax(as.numeric(t), .Machine$double.eps)
    (shape / scale) * (t / scale)^(shape - 1)
  }

  a12 <- function(t) h_weibull(t, shape12, scale12)  # 1 -> 2
  a13 <- function(t) h_weibull(t, shape13, scale13)  # 1 -> 3
  a23 <- function(t) h_weibull(t, shape23, scale23)  # 2 -> 3

  simulate_idm(n, a12, a13, a23, admin_cutoff)
}
