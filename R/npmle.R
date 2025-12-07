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

  # N_AB observations with interval censored times of 1 → 2 transition:
  idx_M <- is_ill
  N_AB <- sum(idx_M)
  L_AB <- data$V_healthy[idx_M]
  R_AB <- data$V_ill[idx_M]
  t_AB <- data$T_obs[idx_M]

  stopifnot(all(L_AB < R_AB))

  # N_A of N_AB also make a transition to state 3: t_A ⊆ t_AB
  idx_N_tilde <- is_ill & is_dead
  N_A <- sum(idx_N_tilde)
  t_A <- data$T_obs[idx_N_tilde]

  N_B <- N_AB - N_A
  stopifnot(N_A <= N_AB)

  # N_C observations with direct transition 1 → 3, no missing transitions T_obs = V_healthy:
  idx_K_tilde <- !is_ill & is_dead & is_T_eq_Vhealthy
  N_C <- sum(idx_K_tilde)
  t_C <- data$T_obs[idx_K_tilde]

  # N_D observations censored in state 1, no missing transitions T_obs = V_healthy:
  idx_J <- !is_dead & !is_ill & is_T_eq_Vhealthy
  N_D <- sum(idx_J)
  t_D <- data$T_obs[idx_J]

  # N_E observations, last seen in state 1 and then seen in state 3:
  idx_U <- !is_ill & is_dead & is_T_gt_Vhealthy
  N_E <- sum(idx_U)
  L_E <- data$V_healthy[idx_U]
  t_E <- data$T_obs[idx_U]

  # N_F observations, last seen in state 1 and then censored:
  idx_C <- !is_ill & !is_dead & is_T_gt_Vhealthy
  N_F <- sum(idx_C)
  L_F <- data$V_healthy[idx_C]
  t_F <- data$T_obs[idx_C]

  t_DF <- c(t_D,t_F)

  ##### N_CE_star: E* - Obs and potential 1 -> 3 ####
  t_CE_star <- unique(c(t_C, t_E))
  # r_C should only count exact observations from case 2 (t_C), not t_E
  # sum(r_C) should equal N_C, not N_C + N_E
  r_C <- as.numeric(table(factor(t_C, levels = t_CE_star)))
  N_CE_star <- length(t_CE_star)

  ##### N_AE_star: T* - Obs and potential entry to state 3 from state 2: 1 -> 2 -> 3 ####
  t_AE_star <- unique(c(t_A, t_E))
  r_A <- as.numeric(table(factor(c(t_A), levels = t_AE_star)))

  N_A_star <- length(unique(t_A))

  N_AE_star <- length(t_AE_star)

  ##### Total: N_ABEFCD = N_AB + N_E + N_F + N_C + N_D ####
  N_ABEFCD <- N_AB + N_E + N_F + N_C + N_D # Total count

  ##### Max 1 -> 2: N_ABEF = N_AB + N_E + N_F ####
  N_ABEF <- N_AB + N_E + N_F # Max number through 2.

  #### Creation of A sets ####

  ##### N_AB: LR_AB := [L_AB, R_AB] ####
  LR_AB <- as.interval(matrix(c(L_AB, R_AB), ncol = 2, byrow = F))

  ##### N_ABE: N_AB < m <= N_ABE := N_AB + N_E, R_{N_AB+u} = t_{N_AB+u} ####
  LR_E <- as.interval(matrix(c(L_E, t_E), ncol = 2, byrow = F))
  N_ABE = N_AB + N_E

  ##### N_ABEF: N_ABE := N_AB + N_E < m <= N_ABEF, R_{N_ABE+c} = t_{N_ABE+c} ####
  LR_F <- as.interval(matrix(c(L_F, t_F), ncol = 2, byrow = F))

  ##### LR_ABEF: LR_AB ∪ LR_E ∪ LR_F ####
  LR_ABEF <- as.interval(rbind(LR_AB, LR_E, LR_F))
  I_union <- get_interval(LR_ABEF)

  #### Data manipulation ####
  ##### I: Q_j = [l_i,r_i] ####

  # s_max = max(t_D, 1 <= j <= N_D)
  s_max <- ifelse(length(t_D) > 0 ,max(t_D),0)

  # R_max = max(R_AB, 1 <= m <= N_ABE)
  R_max <- max(LR_AB[, 2], LR_E[, 2])

  # e*_max = max(e*_k, 1 <= k <= N_CE_star)
  if(length(t_CE_star) > 0) {
    e_star_max <- max(t_CE_star)
  } else {
    e_star_max = 0
  }

  # L_bar ={L_m, 1 <= m <= M'} ∪ {T* ∩ A} ∪ {S_J ∩ A} ∪ {s_max : s_max > R_max ∨ e*_max}
  L_bar <- c(
    LR_ABEF[, 1],
    intersect.interval(I_union, t_AE_star),
    intersect.interval(I_union, t_D),
    na.omit(ifelse(s_max > max(R_max, e_star_max), s_max, NA))
  )

  # R_bar = {R_m, 1 <= m <= N_ABE} ∪ {∞}
  R_bar <- c(LR_AB[, 2], t_E, Inf)

  Q_j <- make_Q(L_bar, R_bar)

  I <- nrow(Q_j)

  ##### N_AE_star: lambda_u and I': z_j ####
  # Comment: I believe we have I' z_j's and N_AE_star: lambda_u
  I_mark <- I + N_CE_star
  data_list <- list(
    # ints
    N_D = N_D, N_F = N_F, N_C = N_C, N_E = N_E, N_A = N_A, N_B = N_B,
    N_A_star = N_A_star, N_AB = N_AB, N_ABE = N_ABE, N_AE_star = N_AE_star,
    N_ABEF = N_ABEF, I = I,
    N_CE_star = N_CE_star, I_mark = I_mark,
    N = N_ABEFCD,     # only if you also use it in C++

    # scalars
    s_max = s_max, R_max = R_max, e_star_max = e_star_max,

    # vectors
    t_D = t_D, L_F = L_F, t_F = t_F, t_C = t_C,
    L_E = L_E, t_E = t_E,
    L_AB = L_AB, R_AB = R_AB,     # then in C++ read x["R_AB"]
    t_AB = t_AB,
    t_CE_star = t_CE_star, t_AE_star = t_AE_star,
    t_DF = t_DF,                  # if you really need t_DF

    r_C = r_C, r_A = r_A,
    t_CE_star = t_CE_star,

    # 2-col matrices
    LR_AB = to_mat(LR_AB), LR_E = to_mat(LR_E), LR_F = to_mat(LR_F),
    LR_ABEF = to_mat(LR_ABEF), Q_j = to_mat(Q_j)
  )

  mod_data_list <- list(
    counts = list(
      N_D = N_D, N_F = N_F, N_C = N_C, N_E = N_E, N_A = N_A,
      N_B = N_B, N_A_star = N_A_star, N_AB = N_AB, N_ABE = N_ABE,
      N_AE_star = N_AE_star, N_ABEF = N_ABEF, I = I,
      N_CE_star = N_CE_star, I_mark = I_mark,
      N = N_ABEFCD
    ),

    maxs = list(
      s_max = s_max, R_max = R_max, e_star_max = e_star_max
    ),

    # vectors
    data_vectors = list(
      t_D = t_D, L_F = L_F, t_F = t_F, t_C = t_C,
      L_E = L_E, t_E = t_E,
      L_AB = L_AB, R_AB = R_AB,     # then in C++ read x["R_AB"]
      t_AB = t_AB,
      t_DF = t_DF,r_C = r_C, r_A = r_A
    ),

    support_point = list(
      Q_j = Q_j, t_CE_star = t_CE_star, t_AE_star = t_AE_star,
      Q = Q
    )
  )

  list(data_list_to_cpp = data_list,
       data_list_to_read = mod_data_list)
}

find_estimator_from_z_and_lambda <- function(grid_points, z_j, lambda_u, Q_j, t_CE_star, t_AE_star) {

  I <- nrow(Q_j)
  I_mark <- I + length(t_CE_star)

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
  F12 <- step_cdf(grid_points, times = Q_j[1:I, 2], masses = z_j[1:I])
  F13 <- step_cdf(grid_points, times = t_CE_star, masses = z_j[(I+1):I_mark])
  F_total <- F12 + F13

  A23 <- step_cdf(grid_points, t_AE_star, masses = lambda_u)

  F12_at_l_i <- step_cdf(Q_j[,1]-1e-6, times = Q_j[1:I, 2], masses = z_j[1:I])
  # A12(s): denominators need F(l_i-) for each i
  denom12 <- 1 - F12_at_l_i
  term12  <- ifelse(denom12 > 0, z_j[1:I] / denom12, 0)
  A12 <- step_cdf(grid_points, times = Q_j[,2], term12)

  F_total_at_e_k <- step_cdf(t_CE_star-1e-6, times = Q_j[1:I, 2], masses = z_j[1:I]) +
    step_cdf(t_CE_star-1e-6, times = t_CE_star, masses = z_j[(I+1):I_mark])
  # A12(s): denominators need F(l_i-) for each i
  denom13 <- 1 - F_total_at_e_k
  term13  <- ifelse(denom13 > 0, z_j[(I+1):I_mark] / denom13, 0)
  A13 <- step_cdf(grid_points, times = t_CE_star, term13)


  # Turn vectors into step functions over grid_points
  stepify <- function(x, y, side = c("left", "right"), extend = TRUE) {
    side <- match.arg(side)
    f <- if (side == "left") 0 else 1
    rule <- if (extend) 2 else 1  # 2 = hold ends constant, 1 = NA outside
    approxfun(x, y, method = "constant", f = f, rule = rule)
  }

  F23 <- function(s,t) {
    stopifnot(length(s) == 1)

    t_star_n_order <- order(t_AE_star)
    t_AE_star <- t_AE_star[t_star_n_order]
    lambda_u <- lambda_u[t_star_n_order]

    eligeble <- s < t_AE_star & t_AE_star <= max(t)

    t_AE_star <- t_AE_star[eligeble]
    lambda_u  <- lambda_u[eligeble]

    cp <- cumprod(1-lambda_u)
    idx <- findInterval(t, t_AE_star)
    c(0,cp)[idx+1]
  }



  estimators <- list(
    A12  = stepify(grid_points, A12, side = "left"),
    A13  = stepify(grid_points, A13, side = "left"),
    A23  = stepify(grid_points, A23, side = "left"),
    F12 = stepify(grid_points, F12, side = "left"),
    F13 = stepify(grid_points, F13, side = "left"),
    F = stepify(grid_points, F_total, side = "left"),
    F23 = F23
  )
  class(estimators) <- c("idm_estimators", class(estimators))
  estimators
}

fit_npmle <- function(data,
                      max_iter = 200,
                      tol = 1e-4,
                      verbose = FALSE,
                      eval_likelihood = FALSE,
                      use_frydman = FALSE) {


  lists <- data_to_list_format(data)

  data_list_to_cpp <- lists$data_list_to_cpp

  mdl_ptr <- make_model_data(data_list_to_cpp)

  z_init <- runif(data_list_to_cpp$I_mark)

  z_init <- z_init/sum(z_init)

  lambda_init <- runif(data_list_to_cpp$N_AE_star,min = 0.00001, max = 0.1)

  fit <- em_fit(mdl_ptr,
                z_init = z_init,
                lambda_init = lambda_init,
                max_iter = max_iter,
                tol = tol,
                verbose = verbose,
                eval_likelihood = eval_likelihood,
                use_frydman)

  estimators <- find_estimator_from_z_and_lambda(
    grid_points = seq(0, max(data$T_obs), length.out = 512),
    z_j = fit$z_j,
    lambda = fit$lambda_u,
    data_list_to_cpp$Q,
    data_list_to_cpp$t_CE_star,
    data_list_to_cpp$t_AE_star,
    data_list_to_cpp$t_CE_star
  )

  list(
    estimators = estimators,
    raw_estimators = list(
      z_j = fit$z_j,
      lambda = fit$lambda_u,
      Q = data_list_to_cpp$Q,
      t_CE_star = data_list_to_cpp$t_CE_star,
      t_AE_star = data_list_to_cpp$t_AE_star,
    ),
    settings = list(
      data = data,
      data_list = lists$data_list_to_read,
      max_iter = max_iter,
      tol = tol,
      alpha_ji = fit$alpha_ji,
      gamma_ji = fit$gamma_ji,
      likelihoods = fit$likelihoods,
      z_history = fit$z_history,
      lambda_history = fit$lambda_history
    ),
    converged = fit$converged)
}



