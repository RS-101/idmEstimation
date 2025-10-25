verify_illness_death <- function(
    n,
    a12, a13, a23, # all: function(t_abs); a23 is CALENDAR-TIME hazard
    t0 = 0,
    t_grid = NULL,
    s_grid = NULL,
    ngrid = 600,
    seed = NULL,
    sim_fun = simulate_illness_death) {
  if (!is.null(seed)) set.seed(seed)

  trapz_cum <- function(x, y) {
    n <- length(x)
    out <- numeric(n)
    for (i in 2:n) out[i] <- out[i-1] + 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1])
    out
  }
  trapz_cols <- function(x, Y) {
    # integrate each column of Y over x by trapezoid
    dx <- diff(x)
    Yu <- Y[-1, , drop = FALSE]
    Yd <- Y[-nrow(Y), , drop = FALSE]
    colSums((Yu + Yd) * rep(dx, times = ncol(Y)) / 2)
  }

  # 1) Simulate data (histograms)
  sim <- sim_fun(n = n, a12 = a12, a13 = a13, a23 = a23, t0 = t0)

  first_time <- ifelse(sim$path == "1->2->3", sim$entry2, sim$t13_direct)
  got_ill   <- is.finite(sim$entry2)
  dir_death <- sim$path == "1->3"
  s_23      <- sim$t23_after[!is.na(sim$t23_after)]
  death_123 <- sim$death_time[sim$path == "1->2->3"]

  # 2) Grids
  if (is.null(t_grid)) {
    t_end <- max(first_time[is.finite(first_time)], na.rm = TRUE)
    t_end <- max(t_end, stats::quantile(first_time, 0.995, na.rm = TRUE))
    t_end <- max(t_end, stats::quantile(sim$death_time, 0.995, na.rm = TRUE))
    t_end <- t_end * 1.25 + 1e-8
    t_grid <- seq(t0, t_end, length.out = ngrid)
  }
  if (is.null(s_grid)) {
    if (length(s_23)) {
      s_end <- max(stats::quantile(s_23, 0.995, na.rm = TRUE), max(s_23, na.rm = TRUE))
      s_end <- s_end * 1.25 + 1e-8
    } else {
      s_end <- 1
    }
    s_grid <- seq(0, s_end, length.out = ngrid)
  }

  # 3) First-event theory on t_grid
  h_sum  <- a12(t_grid) + a13(t_grid)
  H_sum  <- trapz_cum(t_grid, h_sum)
  S1     <- exp(-H_sum)
  f12    <- S1 * a12(t_grid)
  f13    <- S1 * a13(t_grid)
  f_first<- f12 + f13
  F12    <- trapz_cum(t_grid, f12)
  F13    <- trapz_cum(t_grid, f13)
  F_first<- 1 - S1
  P12    <- F12[length(F12)]
  P13    <- F13[length(F13)]

  f12_cond <- if (P12 > 0) f12 / P12 else rep(0, length(f12))
  f13_cond <- if (P13 > 0) f13 / P13 else rep(0, length(f13))
  F12_cond <- if (P12 > 0) F12 / P12 else F12
  F13_cond <- if (P13 > 0) F13 / P13 else F13

  # 4) Calendar-time a23: precompute A23(t)
  tA_max   <- max(t_grid) + max(s_grid)
  tA_grid  <- seq(min(t_grid), tA_max, length.out = max(ngrid, length(t_grid)))
  a23_vals <- a23(tA_grid)
  A23      <- trapz_cum(tA_grid, a23_vals)
  A23_fun  <- stats::approxfun(tA_grid, A23, rule = 2)
  a23_fun  <- stats::approxfun(tA_grid, a23_vals, rule = 2)

  te <- t_grid
  s  <- s_grid

  # 4a) Sojourn s = t - te (mixture over entry-time density f12|ill(te))
  te_plus_s   <- outer(te, s, "+")
  A_teps      <- A23_fun(te_plus_s)
  A_te_mat    <- matrix(A23_fun(te), nrow = length(te), ncol = length(s), byrow = FALSE)
  S23_te_s    <- exp(-(A_teps - A_te_mat))                 # S23(s | te)
  a23_teps    <- a23_fun(te_plus_s)                        # a23(te + s)
  dens_te_s   <- matrix(f12_cond, nrow = length(te), ncol = length(s), byrow = FALSE)
  integrand_S <- S23_te_s * dens_te_s
  integrand_f <- S23_te_s * a23_teps * dens_te_s
  S23_mix     <- if (P12 > 0) trapz_cols(te, integrand_S) else rep(1, length(s))
  f23_mix     <- if (P12 > 0) trapz_cols(te, integrand_f) else rep(0, length(s))
  F23_mix     <- 1 - S23_mix

  # 4b) Death-time density conditional on path 1->2->3 (calendar time)
  # f_T(t | 1->2->3) = ∫_{u<=t} f12|ill(u) * S23(t | u) * a23(t) du
  A_t_col     <- matrix(A23_fun(t_grid), nrow = length(te), ncol = length(t_grid), byrow = TRUE)
  A_te_col    <- matrix(A23_fun(te),     nrow = length(te), ncol = length(t_grid), byrow = FALSE)
  S23_tu      <- exp(-(A_t_col - A_te_col))                # S23(t | u)
  mask_ut     <- outer(te, t_grid, function(u, t) as.numeric(u <= t))
  a23_t_col   <- matrix(a23_fun(t_grid), nrow = length(te), ncol = length(t_grid), byrow = TRUE)
  dens_te_t   <- matrix(f12_cond, nrow = length(te), ncol = length(t_grid), byrow = FALSE)
  integrand_ft<- S23_tu * a23_t_col * dens_te_t * mask_ut
  f_death_123 <- if (P12 > 0) trapz_cols(te, integrand_ft) else rep(0, length(t_grid))
  F_death_123 <- trapz_cum(t_grid, f_death_123)

  # 5) Plots

  # Helper for FD binwidth (like breaks = "FD")
  fd_binwidth <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(NULL)
    bw <- 2 * IQR(x) / (length(x)^(1/3))
    if (is.finite(bw) && bw > 0) bw else NULL
  }

  # Theoretical curves as data frames
  df_A  <- data.frame(t = t_grid, f = f_first)
  df_B  <- if (exists("f12_cond")) data.frame(t = t_grid, f = f12_cond) else NULL
  df_C  <- if (exists("f13_cond")) data.frame(t = t_grid, f = f13_cond) else NULL
  df_D  <- data.frame(s = s_grid,  f = f23_mix)
  df_E  <- data.frame(t = t_grid,  f = f_death_123)

  # (A) First event time
  pA <- ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = first_time, y = ggplot2::after_stat(density)),
                   binwidth = fd_binwidth(first_time),
                   color = "grey40", fill = "grey80") +
    ggplot2::geom_line(data = df_A, ggplot2::aes(t, f, color = "theory f_first(t)"), linewidth = 1) +
    ggplot2::labs(title = "First event time (min of illness/death)", x = "t", y = "Density", color = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")

  # (B) Illness time | illness occurs
  if (any(got_ill)) {
    xB <- sim$entry2[got_ill]
    pB <- ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = xB, y = ggplot2::after_stat(density)),
                     binwidth = fd_binwidth(xB),
                     color = "grey40", fill = "grey80") +
      ggplot2::geom_line(data = df_B, ggplot2::aes(t, f, color = "theory f12|ill(t)"), linewidth = 1) +
      ggplot2::labs(title = "Illness time | illness occurs", x = "t", y = "Density", color = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "top")
  } else {
    pB <- ggplot2::ggplot() +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void(base_size = 11) +
      ggplot2::geom_text(ggplot2::aes(0.5, 0.5, label = "No illness events observed"))
  }

  # (C) Direct death time | direct death
  if (any(dir_death)) {
    xC <- sim$death_time[dir_death]
    pC <- ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = xC, y = ggplot2::after_stat(density)),
                     binwidth = fd_binwidth(xC),
                     color = "grey40", fill = "grey80") +
      ggplot2::geom_line(data = df_C, ggplot2::aes(t, f, color = "theory f13|dir(t)"), linewidth = 1) +
      ggplot2::labs(title = "Direct death time | direct death", x = "t", y = "Density", color = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "top")
  } else {
    pC <- ggplot2::ggplot() +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void(base_size = 11) +
      ggplot2::geom_text(ggplot2::aes(0.5, 0.5, label = "No direct deaths observed"))
  }

  # (D) Ill → death elapsed time (sojourn s)
  if (length(s_23)) {
    pD <- ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = s_23, y = ggplot2::after_stat(density)),
                     binwidth = fd_binwidth(s_23),
                     color = "grey40", fill = "grey80") +
      ggplot2::geom_line(data = df_D, ggplot2::aes(s, f, color = "mixture f_S(s)"), linewidth = 1) +
      ggplot2::labs(title = "Ill → death elapsed time (sojourn s)", x = "s", y = "Density", color = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "top")
  } else {
    pD <- ggplot2::ggplot() +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void(base_size = 11) +
      ggplot2::geom_text(ggplot2::aes(0.5, 0.5, label = "No ill→death times observed"))
  }

  # (E) Death time | path 1→2→3 (calendar time)
  if (length(death_123)) {
    pE <- ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = death_123, y = ggplot2::after_stat(density)),
                     binwidth = fd_binwidth(death_123),
                     color = "grey40", fill = "grey80") +
      ggplot2::geom_line(data = df_E, ggplot2::aes(t, f, color = "theory f_T(t | 1→2→3)"), linewidth = 1) +
      ggplot2::labs(title = "Death time | path 1→2→3 (calendar time)", x = "t", y = "Density", color = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "top")
  } else {
    pE <- ggplot2::ggplot() +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void(base_size = 11) +
      ggplot2::geom_text(ggplot2::aes(0.5, 0.5, label = "No 1→2→3 paths observed"))
  }

  # (F) Summary panel as a geom_text plot
  emp_P12 <- mean(got_ill)
  emp_P13 <- mean(dir_death)
  txtF <- sprintf(
    "Estimated path probabilities from simulation (n=%d):\n  P(1→2→3)  empirical = %.4f ; theory = %.4f\n  P(1→3)    empirical = %.4f ; theory = %.4f\n  Check: P12+P13 empirical = %.4f ; theory = %.4f",
    n, emp_P12, P12, emp_P13, P13, emp_P12 + emp_P13, P12 + P13
  )

  pF <- ggplot2::ggplot() +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void(base_size = 11) +
    ggplot2::geom_text(ggplot2::aes(x = 0, y = 1, label = txtF), hjust = 0, vjust = 1) +
    ggplot2::labs(title = "Summary: P12 & P13 shown below")

  # Compose 2×3 layout with patchwork
  my_plot <- (pA | pB | pC) /
    (pD | pE | pF)

  hazards = list(a12 = a12, a13 = a13, a23 = a23)
  class(hazards) <- c(class(hazards), "full_hazard")

  # 6) Return evaluators
  list(
    sim = sim,
    grid = list(t = t_grid, s = s_grid),
    # First-event
    pdf_first = stats::approxfun(t_grid, f_first, rule = 2),
    cdf_first = stats::approxfun(t_grid, F_first, rule = 2),
    # Illness time | illness
    pdf_ill_cond = stats::approxfun(t_grid, f12_cond, rule = 2),
    cdf_ill_cond = stats::approxfun(t_grid, F12_cond, rule = 2),
    # Direct death | direct
    pdf_dirdeath_cond = stats::approxfun(t_grid, f13_cond, rule = 2),
    cdf_dirdeath_cond = stats::approxfun(t_grid, F13_cond, rule = 2),
    # Sojourn s density and CDF
    pdf_ill2death_sojourn = stats::approxfun(s_grid, f23_mix, rule = 2),
    cdf_ill2death_sojourn = stats::approxfun(s_grid, F23_mix, rule = 2),
    # Death time | path 1->2->3 (calendar time)
    pdf_death_given_123 = stats::approxfun(t_grid, f_death_123, rule = 2),
    cdf_death_given_123 = stats::approxfun(t_grid, F_death_123, rule = 2),
    hazards = hazards,
    theory = list(
      t_grid = t_grid, S1 = S1, f12 = f12, f13 = f13, f_first = f_first,
      F12 = F12, F13 = F13, F_first = F_first, P12 = P12, P13 = P13,
      s_grid = s_grid, S23_sojourn = S23_mix
    ),
    plot = my_plot
  )
}
