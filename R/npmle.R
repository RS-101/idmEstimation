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



  case_data_list <- create_case_data(data)


  N_A <- length(case_data_list$case_A$T_obs)
  N_B <- length(case_data_list$case_B$T_obs)
  N_C <- length(case_data_list$case_C$T_obs)
  N_D <- length(case_data_list$case_D$T_obs)
  N_E <- length(case_data_list$case_E$T_obs)
  N_F <- length(case_data_list$case_F$T_obs)

  # N_AB observations with interval censored times of 1 → 2 transition:
  N_AB <- N_A + N_B
  if("case_A" %in% names(case_data_list)) {
    N_A <- length(case_data_list$case_A$T_obs)
    L_A <- case_data_list$case_A$V_healthy
    R_A <- case_data_list$case_A$V_ill
    t_A <- case_data_list$case_A$T_obs
  } else {
    N_A <- 0
    L_A <- c()
    R_A <- c()
    t_A <- c()
  }

  if("case_B" %in% names(case_data_list)) {
    N_B <- length(case_data_list$case_B$T_obs)
    L_B <- case_data_list$case_B$V_healthy
    R_B <- case_data_list$case_B$V_ill
    t_B <- case_data_list$case_B$T_obs
  } else {
    N_B <- 0
    L_B <- c()
    R_B <- c()
    t_B <- c()
  }

  if("case_C" %in% names(case_data_list)) {
    N_C <- length(case_data_list$case_C$T_obs)
    t_C <- case_data_list$case_C$T_obs
  } else {
    N_C <- 0
    t_C <- c()
  }

  if("case_D" %in% names(case_data_list)) {
    N_D <- length(case_data_list$case_D$T_obs)
    t_D <- case_data_list$case_D$T_obs
  } else {
    N_D <- 0
    t_D <- c()
  }

  if("case_E" %in% names(case_data_list)) {
    N_E <- length(case_data_list$case_E$T_obs)
    L_E <- case_data_list$case_E$V_healthy
    t_E <- case_data_list$case_E$T_obs
  } else {
    N_E <- 0
    L_E <- c()
    t_E <- c()
  }

  if("case_F" %in% names(case_data_list)) {
    N_F <- length(case_data_list$case_F$T_obs)
    L_F <- case_data_list$case_F$V_healthy
    t_F <- case_data_list$case_F$T_obs
  } else {
    N_F <- 0
    L_F <- c()
    t_F <- c()
  }


  N_AB <- N_A + N_B
  L_AB <- c(L_A, L_B)
  R_AB <- c(R_A, R_B)
  t_AB <- c(t_A, t_B)

  stopifnot(all(L_AB < R_AB))
  stopifnot(N_A <= N_AB)

  t_DF <- c(t_D, t_F)

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
  C_D_max <- ifelse(length(t_D) > 0 ,max(t_D),0)

  # R_max = max(R_AB, 1 <= m <= N_ABE)
  R_max <- max(LR_AB[, 2], LR_E[, 2])

  # e*_max = max(e*_k, 1 <= k <= N_CE_star)
  if(length(t_CE_star) > 0) {
    T_CE_star_max <- max(t_CE_star)
  } else {
    T_CE_star_max = 0
  }

  L_F_max <- max(L_F)

  # L_bar ={L_m, 1 <= m <= M'} ∪ {T* ∩ A} ∪ {S_J ∩ A} ∪ {s_max : s_max > R_max ∨ e*_max}
  L_bar <- c(
    L_AB, L_E,
    intersect.interval(I_union, t_AE_star),
    intersect.interval(I_union, t_D),
    na.omit(ifelse(L_F <= max(R_max, T_CE_star_max), L_F, NA )),
    na.omit(ifelse(C_D_max >= max(R_max, T_CE_star_max, L_F_max), C_D_max, NA)),
    na.omit(ifelse(L_F_max >= max(R_max, T_CE_star_max, C_D_max), L_F_max, NA))
  )

  # R_bar = {R_m, 1 <= m <= N_ABE} ∪ {∞}
  R_bar <- c(R_AB, t_E, Inf)

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
    s_max = C_D_max, R_max = R_max, T_CE_star_max = T_CE_star_max,

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
      s_max = C_D_max, R_max = R_max, T_CE_star_max = T_CE_star_max
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
      Q_j = Q_j, t_CE_star = t_CE_star, t_AE_star = t_AE_star
    )
  )

  list(data_list_to_cpp = data_list,
       data_list_to_read = mod_data_list)
}

create_npmle_estimators <- function(z_j, lambda_u, Q_j, t_CE_star, t_AE_star) {
  I <- nrow(Q_j)
  I_mark <- I + length(t_CE_star)

  y_k <- z_j[(I + 1):I_mark]
  z_j <- z_j[1:I]

  # we make sure that y_k and lambda_u are in the order of t_CE_star and t_AE_star
  t_CE_star_order <- order(t_CE_star)
  t_CE_star <- t_CE_star[t_CE_star_order]
  y_k <- y_k[t_CE_star_order]
  t_AE_star_order <- order(t_AE_star)
  t_AE_star <- t_AE_star[t_AE_star_order]
  lambda_u <- lambda_u[t_AE_star_order]

  q_m <- Q_j[I, 1]
  p_m_is_inf <- is.infinite(Q_j[I, 2])

  F13 <- function(t) {
    # We calculate sum 1(t_CE_star <= t) * y_k
    # create the matrix of t_CE_star < t
    eligeble <- outer(t_CE_star, t, FUN = "<=")
    # multiply each row with y_k
    weighted <- eligeble * matrix(y_k, nrow = length(t_CE_star), ncol = length(t), byrow = FALSE)
    # sum over rows
    retval <- colSums(weighted)
    if(p_m_is_inf)  retval[t > q_m] <- NA

    retval
  }

  t_AE_star_max <- max(t_AE_star)

  A23 <- function(t) {
    # if t > t_AE_star_max return NA

    # We calculuate sum 1(t_AE_star <= t) * lambda_u
    # create the matrix of t_AE_star < t
    eligeble <- outer(t_AE_star, t, FUN = "<=")
    # multiply each row with lambda_u
    weighted <- eligeble * matrix(lambda_u, nrow = length(t_AE_star), ncol = length(t), byrow = FALSE)
    # sum over rows
    retval <- colSums(weighted)

    retval[t > t_AE_star_max] <- NA
    retval
  }

  F12 <- function(t) {
    # cumulative sum part: sum z_j for all j with p_j <= t
    eligeble <- outer(Q_j[, 2], t, FUN = "<=")
    weighted <- eligeble * matrix(z_j, nrow = I, ncol = length(t), byrow = FALSE)
    retval   <- colSums(weighted)

    # outside global range
    retval[t < Q_j[1, 1]] <- 0
    retval[t > Q_j[I, 2]] <- NA

    # set NA inside each interval [q_j, p_j]
    for (j in seq_len(I)) {
      inside <- t > Q_j[j, 1] & t <= Q_j[j, 2]
      retval[inside] <- NA
    }

    retval
  }

  F12_at_q_j <- F12(Q_j[,1])

  A12 <- function(t) {
    # cumulative sum part: sum z_j for all j with p_j <= t
    eligeble <- outer(Q_j[, 2], t, FUN = "<=")
    weighted <- eligeble * matrix(z_j / (1 - F12_at_q_j), nrow = I, ncol = length(t), byrow = FALSE)
    retval   <- colSums(weighted)

    # outside global range
    retval[t < Q_j[1, 1]] <- 0
    retval[t > Q_j[I, 2]] <- NA

    retval
  }

  F13_at_t_CE_star <- F13(t_CE_star - 1e-8)

  t_CE_star_max <- max(t_CE_star)

  A13 <- function(t) {
    # We calculuate sum 1(t_CE_star <= t) * lambda_u / (1 - F13(t_CE_star))
    eligeble <- outer(t_CE_star, t, FUN = "<=")
    weighted <- eligeble * matrix(y_k / (1 - F13_at_t_CE_star), nrow = length(t_CE_star), ncol = length(t), byrow = FALSE)
    retval <- colSums(weighted)

    retval[t > t_CE_star_max] <- NA
    retval
  }


  P22 <- function(t, entry_time = 1e-4) {
    stopifnot(length(entry_time) == 1)

    eligeble <- entry_time < t_AE_star & t_AE_star <= max(t)

    t_AE_star_filtered <- t_AE_star[eligeble]
    lambda_u_filtered  <- lambda_u[eligeble]

    # Handle edge case: no eligible time points
    if (length(t_AE_star_filtered) == 0) {
      return(rep(NA, length(t)))
    }

    # Compute survival function: S(t) = prod(1 - lambda_u) for all u <= t
    survival <- cumprod(1 - lambda_u_filtered)

    # P22(t) = S(t) = survival probability in state 2
    idx <- findInterval(t, t_AE_star_filtered, rightmost.closed = TRUE)

    retval <- numeric(length(t))
    # Before first event after entry: P22 = 1 (still in state 2)
    retval[idx == 0] <- 1
    # After first event: P22 = S(t)
    valid_idx <- idx > 0
    retval[valid_idx] <- survival[idx[valid_idx]]
    # After last event: set to NA if needed
    retval[t > max(t_AE_star)] <- NA

    retval
  }



  estimators <- list(
    distribution_functions = list(
      F12 = F12,
      F13 = F13,
      P22 = P22
    ),
    cum_hazard_functions = list(
      A12 = A12,
      A13 = A13,
      A23 = A23
    )
  )
  class(estimators) <- c("idm_estimators", class(estimators))
  estimators
}

fit_npmle <- function(data,
                      max_iter = 200,
                      tol = 1e-4,
                      verbose = FALSE,
                      eval_likelihood = FALSE,
                      use_true_EM = FALSE) {


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
                use_true_EM)

  estimators <- create_npmle_estimators(
    z_j = fit$z_j,
    lambda_u = fit$lambda_u,
    Q_j = data_list_to_cpp$Q_j,
    t_CE_star = data_list_to_cpp$t_CE_star,
    t_AE_star = data_list_to_cpp$t_AE_star
  )



  idm_fit <- list(
    estimators = estimators,
    data = data,
    model_type = "non_parametric_mle",
    model_config = list(
      data_list = lists$data_list_to_read,
      max_iter = max_iter,
      tol = tol
    ),
    raw_estimators = list(
      z_j = fit$z_j[1:data_list_to_cpp$I],
      z_k = fit$z_j[(data_list_to_cpp$I + 1):data_list_to_cpp$I_mark],
      lambda_u = fit$lambda_u,
      Q_j = data_list_to_cpp$Q_j,
      t_CE_star = data_list_to_cpp$t_CE_star,
      t_AE_star = data_list_to_cpp$t_AE_star
    ),
    model_specific = list(
      alpha_ji = fit$alpha_ji,
      gamma_ji = fit$gamma_ji,
      likelihoods = fit$likelihoods,
      z_history = fit$z_history,
      lambda_history = fit$lambda_history
    ),
    converged = fit$converged
  )
  class(idm_fit) <- c("idm_object", class(idm_fit))
  return(idm_fit)
}



