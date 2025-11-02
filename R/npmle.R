data_to_list_format <- function(data, is_equal_tol = 1e-8) {

  intersect.interval <- function(x, y) {
    if (inherits(y, "interval") & inherits(x, "numeric")) {
      x_temp <- x
      x <- y
      y <- x_temp
    }

    if(inherits(y, "numeric")) {
      L <- matrix(rep(x[,1], length(y)),ncol = nrow(x), byrow = T)
      R <- matrix(rep(x[,2], length(y)),ncol = nrow(x), byrow = T)

      if (inherits(x, "c_c")) {
        sel <- L <= y & y <= R
      } else if (inherits(x, "c_o")) {
        sel <- L <= y & y < R
      } else if (inherits(x, "o_c")) {
        sel <- L < y & y <= R
      } else if (inherits(x, "o_o")) {
        sel <- L < y & y < R
      }
      res = y[1 <= rowSums(sel)]
      return(res)
    } else if (inherits(y, "interval")) {
      stop("not implemented")
      if (inherits(x, "c_c")) {

      } else if (inherits(x, "c_o")) {

      } else if (inherits(x, "o_c")) {

      } else if (inherits(x, "o_o")) {

      }
    }
    stop("y should be numeric or interval")
  }

  get_interval <- function(m, L_open = F, R_open = F) {
    if (ncol(m) != 2) stop("LR mat must be of dim m x 2")
    type = paste(ifelse(L_open, "o", "c"), ifelse(R_open, "o", "c"), sep = "_")
    if(nrow(m) > 1) {

      m <- unique(m)
      m_sorted <- m[order(m[, 1]), ]

      L <- m_sorted[,1][-1]
      R <- m_sorted[,2][-nrow(m)]


      if(L_open & R_open) {
        start_stop <- !(L < R)
      } else {
        start_stop <- !(L <= R)
      }
      interval_start <- c(min(m_sorted[, 1]),L[start_stop])
      interval_end <- c(R[start_stop], max(m_sorted[, 2]))
      res <- matrix(c(interval_start, interval_end), byrow = F, ncol = 2)
    } else {
      res <- m
    }
    class(res) <- c("interval", type, class(res))
    res
  }


  ###### as.interval ######

  as.interval <- function(x, L_open = F, R_open = F) {
    if (ncol(x) != 2) stop("LR mat must be of dim m x 2")
    type = paste(ifelse(L_open, "o", "c"), ifelse(R_open, "o", "c"), sep = "_")

    class(x) <- c("interval", type, class(x))
    x
  }

  ##### From paper specific #####

  make_Q <- function(L_bar, R_bar) {
    L_bar <- sort(L_bar[!is.na(L_bar)])
    R_bar <- sort(R_bar[!is.na(R_bar)])
    Q <- matrix(c(rep(0L, length(L_bar)), rep(1L, length(R_bar)),
                  L_bar, R_bar), ncol = 2)
    Q <- Q[order(Q[,2]), ]
    tag <- which(diff(Q[, 1], 1) == 1)
    Q <- matrix(c(Q[tag, 2], Q[tag + 1, 2]), ncol = 2)

    Q <- as.interval(Q, L_open = F, R_open = F)
    Q
  }






  to_mat <- function(x) if (is.matrix(x)) x else as.matrix(unclass(x))


  # check if data has correct columns
  required_cols <- c("id", "T_obs", "V_healthy", "V_ill", "status_dead", "status_ill")
  missing_cols <- setdiff(required_cols, names(data))
  if(length(missing_cols) > 0) {
    stop(paste("Data is missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Create logical index vectors for each observation type
  is_ill <- data$status_ill == 1
  is_dead <- data$status_dead == 1
  is_T_eq_Vhealthy <- abs(data$T_obs - data$V_healthy) < is_equal_tol
  is_T_gt_Vhealthy <- abs(data$T_obs - data$V_healthy) > is_equal_tol

  # M observations with interval censored times of 1 → 2 transition:
  idx_M <- is_ill
  M <- sum(idx_M)
  L_m <- data$V_healthy[idx_M]
  R_m <- data$V_ill[idx_M]
  t_m <- data$T_obs[idx_M]

  stopifnot(all(L_m < R_m))

  # N_tilde of M also make a transition to state 3: t_m_in_N_tilde ⊆ t_m
  idx_N_tilde <- is_ill & is_dead
  N_tilde <- sum(idx_N_tilde)
  t_m_in_N_tilde <- data$T_obs[idx_N_tilde]

  stopifnot(N_tilde < M)

  # K_tilde observations with direct transition 1 → 3, no missing transitions T_obs = V_healthy:
  idx_K_tilde <- !is_ill & is_dead & is_T_eq_Vhealthy
  K_tilde <- sum(idx_K_tilde)
  e_k <- data$T_obs[idx_K_tilde]

  # J observations censored in state 1, no missing transitions T_obs = V_healthy:
  idx_J <- !is_dead & !is_ill & is_T_eq_Vhealthy
  J <- sum(idx_J)
  s_j <- data$T_obs[idx_J]

  # U observations, last seen in state 1 and then seen in state 3:
  idx_U <- !is_ill & is_dead & is_T_gt_Vhealthy
  U <- sum(idx_U)
  L_u <- data$V_healthy[idx_U]
  t_u <- data$T_obs[idx_U]

  # C observations, last seen in state 1 and then censored:
  idx_C <- !is_ill & !is_dead & is_T_gt_Vhealthy
  C <- sum(idx_C)
  L_c <- data$V_healthy[idx_C]
  t_c <- data$T_obs[idx_C]


  ##### K: E* - Obs and potential 1 -> 3 ####
  E_star <- unique(c(e_k, t_u))
  # c_k should only count exact observations from case 2 (e_k), not t_u
  # sum(c_k) should equal K_tilde, not K_tilde + U
  c_k <- as.numeric(table(factor(e_k, levels = E_star)))
  K <- length(E_star)

  ##### N: T* - Obs and potential entry to state 3 from state 2: 1 -> 2 -> 3 ####
  t_star_n <- unique(c(t_m_in_N_tilde, t_u))
  d_n <- as.numeric(table(factor(c(t_m_in_N_tilde, t_u), levels = t_star_n)))

  N1_obs_of_T_star <- length(unique(t_m_in_N_tilde))
  U_pos_obs_of_T_star <- length(setdiff(unique(t_u), unique(t_m_in_N_tilde)))

  N <- length(t_star_n)

  ##### Total: N* = M + U + C + K_tilde + J ####
  N_star <- M + U + C + K_tilde + J # Total count

  ##### Max 1 -> 2: M' = M + U + C ####
  M_mark <- M + U + C # Max number through 2.

  #### Creation of A sets ####

  ##### M: A_m := [L_m, R_m] ####
  A_m <- as.interval(matrix(c(L_m, R_m), ncol = 2, byrow = F))

  ##### W: M < m <= W := M + U, R_{M+u} = t_{M+u} ####
  A_u <- as.interval(matrix(c(L_u, t_u), ncol = 2, byrow = F))
  W = M + U

  ##### M': W := M + U < m <= M', R_{W+c} = t_{W+c} ####
  A_c <- as.interval(matrix(c(L_c, t_c), ncol = 2, byrow = F))

  ##### full_A_m: A_m ∪ A_u ∪ A_c ####
  full_A_m <- as.interval(rbind(A_m, A_u, A_c))

  ##### A := ⋃_{m=1}^{M'} A_m ####
  A_union <- get_interval(full_A_m)

  #### Data manipulation ####
  ##### I: Q_i = [l_i,r_i] ####

  # s_max = max(s_j, 1 <= j <= J)
  s_max <- ifelse(length(s_j) > 0 ,max(s_j),0)

  # R_max = max(R_m, 1 <= m <= W)
  R_max <- max(A_m[, 2], A_u[, 2])

  # e*_max = max(e*_k, 1 <= k <= K)
  e_star_max <- max(E_star)

  # L_bar ={L_m, 1 <= m <= M'} ∪ {T* ∩ A} ∪ {S_J ∩ A} ∪ {s_max : s_max > R_max ∨ e*_max}
  L_bar <- c(
    full_A_m[, 1],
    intersect.interval(A_union, t_star_n),
    intersect.interval(A_union, s_j),
    na.omit(ifelse(s_max > max(R_max, e_star_max), s_max, NA))
  )

  # R_bar = {R_m, 1 <= m <= W} ∪ {∞}
  R_bar <- c(full_A_m[1:(M + U), 2], Inf)

  # !!!DANGER I AM UNSURE ABOUT THE CREATION OF Q!!!
  Q_i <- make_Q(L_bar, R_bar)

  Q <- get_interval(Q_i)

  I <- nrow(Q_i)


  ##### I'= K + I: Q_i' = e*_i-I ####
  Q_i_mark <- E_star

  ##### Q_full = list ####
  Q_full <- list(Q_i, Q_i_mark)

  ##### C: s_J+c = t_W+c ####
  s_j_c <- t_c

  s_j_full <- c(s_j, s_j_c)


  ##### N: lambda_n and I': z_i ####
  # Comment: I believe we have I' z_i's and N: lambda_n
  I_mark <- I + K

  data_list <- list(
    # ints
    J = J, C = C, K_tilde = K_tilde, U = U, N_tilde = N_tilde, M = M, W = W,
    N1_obs_of_T_star = N1_obs_of_T_star,
    U_pos_obs_of_T_star = U_pos_obs_of_T_star,
    N = N, N_star = N_star, M_mark = M_mark, I = I, K = K, I_mark = I_mark,

    # scalars
    s_max = s_max, R_max = R_max, e_star_max = e_star_max,

    # vectors
    s_j = s_j, L_c = L_c, t_c = t_c, e_k = e_k, L_u = L_u, t_u = t_u,
    t_m_in_N_tilde = t_m_in_N_tilde, L_m = L_m, R_m = R_m, t_m = t_m,
    E_star = E_star, t_star_n = t_star_n, c_k = c_k, d_n = d_n,
    L_bar = L_bar, R_bar = R_bar, s_j_c = s_j_c, s_j_full = s_j_full,
    Q = Q, Q_i_mark = Q_i_mark, A_union = A_union,

    # 2-col “interval” matrices
    A_m = to_mat(A_m), A_u = to_mat(A_u), A_c = to_mat(A_c),
    full_A_m = to_mat(full_A_m), Q_i = to_mat(Q_i)
  )
  
  data_list
}

find_estimator_from_z_and_lambda <- function(grid_points, z_i, lambda_n, Q_i, Q_i_mark, t_star_n, E_star) {

  I <- nrow(Q_i)
  I_mark <- I + length(Q_i_mark)

  # helper: right-continuous step CDF from (times, masses) at s_eval
  step_cdf <- function(grid_points, times, masses, ...) {
    time_order <- order(times)
    times <- times[time_order]
    masses <- masses[time_order]
    cs <- pmax(cumsum(masses),0)
    idx <- findInterval(grid_points, times, ...)
    c(0,cs)[idx+1]
  }

  # use intercept for F12 instead of just right or left endpoint
  F12 <- step_cdf(grid_points, times = Q_i[1:I, 2], masses = z_i[1:I])
  F13 <- step_cdf(grid_points, times = Q_i_mark, masses = z_i[(I+1):I_mark])
  F_total <- F12 + F13

  A23 <- step_cdf(grid_points, t_star_n, masses = lambda_n)

  F12_at_l_i <- step_cdf(Q_i[,1]-1e-6, times = Q_i[1:I, 2], masses = z_i[1:I])
  # A12(s): denominators need F(l_i-) for each i
  denom12 <- 1 - F12_at_l_i
  term12  <- ifelse(denom12 > 0, z_i[1:I] / denom12, 0)
  A12 <- step_cdf(grid_points, times = Q_i[,2], term12)

  F_total_at_e_k <- step_cdf(E_star-1e-6, times = Q_i[1:I, 2], masses = z_i[1:I]) +
    step_cdf(E_star-1e-6, times = Q_i_mark, masses = z_i[(I+1):I_mark])
  # A12(s): denominators need F(l_i-) for each i
  denom13 <- 1 - F_total_at_e_k
  term13  <- ifelse(denom12 > 0, z_i[(I+1):I_mark] / denom13, 0)
  A13 <- step_cdf(grid_points, times = E_star, term13)


  # Turn vectors into step functions over grid_points
  stepify <- function(x, y, side = c("left", "right"), extend = TRUE) {
    side <- match.arg(side)
    f <- if (side == "left") 0 else 1
    rule <- if (extend) 2 else 1  # 2 = hold ends constant, 1 = NA outside
    approxfun(x, y, method = "constant", f = f, rule = rule)
  }

  hazards_as_functions <- list(
    A12  = stepify(grid_points, A12, side = "left"),
    A13  = stepify(grid_points, A13, side = "left"),
    A23  = stepify(grid_points, A23, side = "left")
  )
  class(hazards_as_functions) <- c("idm_hazards", class(hazards_as_functions))


  list(
    distributions_as_functions = list(
      F12 = stepify(grid_points, F12, side = "left"),
      F13 = stepify(grid_points, F13, side = "left"),
      F = stepify(grid_points, F_total, side = "left")),
    hazards_as_functions = hazards_as_functions
  )
}

fit_npmle <- function(data,
                      max_iter = 200,
                      tol = 1e-4,
                      verbose = FALSE) {


  data_list <- data_to_list_format(data)
  mdl_ptr <- make_model_data(data_list)

  z_init <- runif(data_list$I_mark)

  z_init <- z_init/sum(z_init)

  lambda_init <- runif(data_list$N,min = 0.00001, max = 0.05)

  fit <- em_fit(mdl_ptr,
                z_init = z_init,
                lambda_init = lambda_init,
                max_iter = max_iter,
                tol = tol,
                verbose = verbose)

  estimators <- find_estimator_from_z_and_lambda(
    grid_points = seq(0, max(data$T_obs), length.out = 512),
    z_i = fit$z_i,
    lambda = fit$lambda_n,
    data_list$Q,
    data_list$Q_i_mark,
    data_list$t_star_n,
    data_list$E_star
  )

  list(
    hazards = estimators$hazards_as_functions,
    distributions_functions = estimators$distributions_as_functions,
    raw_estimators = list(
      z_i = fit$z_i,
      lambda = fit$lambda_n
    ),
    settings = list(
      data = data,
      max_iter = max_iter,
      tol = tol
    ))
}
