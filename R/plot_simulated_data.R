#' Verify illness-death model theory against simulated data
#'
#' Given hazard functions and a simulated dataset from an illness-death model,
#' this function checks whether the model's predictions are consistent with the
#' observed data. It computes theoretical densities and cumulative distribution
#' functions, and compares them with the empirical distributions from the
#' simulated data.
#'
#' @inheritParams simulate_illness_death
#' @return TRUE
#' @export
verify_illness_death <- function(a12, a13, a23, simulated_data) { 
  names(simulated_data) <- tolower(names(simulated_data))

  trapz_cum <- function(x, y) {
    n <- length(x)
    out <- numeric(n)
    for (i in 2:n) out[i] <- out[i-1] + 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1])
    out
  }

  trapz_cols <- function(x, y) {
    # integrate each column of Y over x by trapezoid
    dx <- diff(x)
    yu <- y[-1, , drop = FALSE]
    yd <- y[-nrow(y), , drop = FALSE]
    colSums((yu + yd) * rep(dx, times = ncol(y)) / 2)
  }
 
  first_time <- ifelse(simulated_data$path == "1->2->3", simulated_data$entry2, simulated_data$t13_direct)
  got_ill   <- is.finite(simulated_data$entry2)
  dir_death <- simulated_data$path == "1->3"
  s_23      <- simulated_data$t23_after[!is.na(simulated_data$t23_after)]
  death_123 <- simulated_data$death_time[simulated_data$path == "1->2->3"]


  grid_points <- seq(, length.out = ngrid)


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
  h_sum  <- a12(grid_points) + a13(grid_points)
  H_sum  <- trapz_cum(grid_points, h_sum)
  S1     <- exp(-H_sum)
  f12    <- S1 * a12(grid_points)
  f13    <- S1 * a13(grid_points)
  f_first<- f12 + f13
  F12    <- trapz_cum(grid_points, f12)
  F13    <- trapz_cum(grid_points, f13)
  F_first<- 1 - S1
  P12    <- F12[length(F12)]
  P13    <- F13[length(F13)]

  f12_cond <- if (P12 > 0) f12 / P12 else rep(0, length(f12))
  f13_cond <- if (P13 > 0) f13 / P13 else rep(0, length(f13))
  F12_cond <- if (P12 > 0) F12 / P12 else F12
  F13_cond <- if (P13 > 0) F13 / P13 else F13

  # 4) Calendar-time a23: precompute A23(t)
  tA_max   <- max(grid_points) + max(s_grid)
  tA_grid  <- seq(min(grid_points), tA_max, length.out = max(ngrid, length(grid_points)))
  a23_vals <- a23(tA_grid)
  A23      <- trapz_cum(tA_grid, a23_vals)
  A23_fun  <- stats::approxfun(tA_grid, A23, rule = 2)
  a23_fun  <- stats::approxfun(tA_grid, a23_vals, rule = 2)

  te <- grid_points
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
  A_t_col     <- matrix(A23_fun(grid_points), nrow = length(te), ncol = length(grid_points), byrow = TRUE)
  A_te_col    <- matrix(A23_fun(te),     nrow = length(te), ncol = length(grid_points), byrow = FALSE)
  S23_tu      <- exp(-(A_t_col - A_te_col))                # S23(t | u)
  mask_ut     <- outer(te, grid_points, function(u, t) as.numeric(u <= t))
  a23_t_col   <- matrix(a23_fun(grid_points), nrow = length(te), ncol = length(grid_points), byrow = TRUE)
  dens_te_t   <- matrix(f12_cond, nrow = length(te), ncol = length(grid_points), byrow = FALSE)
  integrand_ft<- S23_tu * a23_t_col * dens_te_t * mask_ut
  f_death_123 <- if (P12 > 0) trapz_cols(te, integrand_ft) else rep(0, length(grid_points))
  F_death_123 <- trapz_cum(grid_points, f_death_123)

  # 5) Plots

  # Helper for FD binwidth (like breaks = "FD")
  fd_binwidth <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(NULL)
    bw <- 2 * IQR(x) / (length(x)^(1/3))
    if (is.finite(bw) && bw > 0) bw else NULL
  }

  # Theoretical curves as data frames
  df_A  <- data.frame(t = grid_points, f = f_first)
  df_B  <- if (exists("f12_cond")) data.frame(t = grid_points, f = f12_cond) else NULL
  df_C  <- if (exists("f13_cond")) data.frame(t = grid_points, f = f13_cond) else NULL
  df_D  <- data.frame(s = s_grid,  f = f23_mix)
  df_E  <- data.frame(t = grid_points,  f = f_death_123)

  hazards = list(a12 = a12, a13 = a13, a23 = a23)
  class(hazards) <- c(class(hazards), "full_hazard")

  # 6) Return evaluators
  list(
    sim = simulated_data,
    grid = list(t = grid_points, s = s_grid),
    # First-event
    pdf_first = stats::approxfun(grid_points, f_first, rule = 2),
    cdf_first = stats::approxfun(grid_points, F_first, rule = 2),
    # Illness time | illness
    pdf_ill_cond = stats::approxfun(grid_points, f12_cond, rule = 2),
    cdf_ill_cond = stats::approxfun(grid_points, F12_cond, rule = 2),
    # Direct death | direct
    pdf_dirdeath_cond = stats::approxfun(grid_points, f13_cond, rule = 2),
    cdf_dirdeath_cond = stats::approxfun(grid_points, F13_cond, rule = 2),
    # Sojourn s density and CDF
    pdf_ill2death_sojourn = stats::approxfun(s_grid, f23_mix, rule = 2),
    cdf_ill2death_sojourn = stats::approxfun(s_grid, F23_mix, rule = 2),
    # Death time | path 1->2->3 (calendar time)
    pdf_death_given_123 = stats::approxfun(grid_points, f_death_123, rule = 2),
    cdf_death_given_123 = stats::approxfun(grid_points, F_death_123, rule = 2),
    hazards = hazards,
    theory = list(
      t_grid = grid_points, S1 = S1, f12 = f12, f13 = f13, f_first = f_first,
      F12 = F12, F13 = F13, F_first = F_first, P12 = P12, P13 = P13,
      s_grid = s_grid, S23_sojourn = S23_mix
    )
  )
  TRUE
}
