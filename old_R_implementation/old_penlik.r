library(splines2)
#### Create spline hazard ####
make_spline_funs <- function(param, knots, degree) {

  spline_hazard <- function(x, theta, knots, degree = 3, cumulative = F) {
    n_knots <- length(knots)
    if (length(theta) != degree + n_knots - 1) stop("theta wrong length")

    if (cumulative) {
      I <- iSpline(x,
                   degree = degree,
                   knots = knots[-c(1,n_knots)],
                   Boundary.knots = knots[c(1,n_knots)],
                   intercept = T,
                   warn.outside = F)
      return(as.vector(I %*% theta))
    } else {
      M <- mSpline(x,
                   degree = degree,
                   knots = knots[-c(1,n_knots)],
                   Boundary.knots = knots[c(1,n_knots)],
                   intercept = T,
                   warn.outside = F)
      return(as.vector(M %*% theta))
    }
  }



  return_obj <- list(
    hazard = function(x) spline_hazard(x, param, knots, degree, cumulative = FALSE),
    cumhaz = function(x) spline_hazard(x, param, knots, degree, cumulative = TRUE)
  )


  class(return_obj) <- c(class(return_obj), "single_spline_hazard")
  return_obj
}


# we assume order a12, a13 and a23
make_all_3_spline_funs <- function(theta, knots, degree) {
  # define all transitions
  a12_funs <- make_spline_funs(theta$a12, knots$a12, degree)
  a13_funs <- make_spline_funs(theta$a13, knots$a13, degree)
  a23_funs <- make_spline_funs(theta$a23, knots$a23, degree)

  # hazards
  a12 <- a12_funs$hazard
  a13 <- a13_funs$hazard
  a23 <- a23_funs$hazard

  # cumulative hazards
  A12 <- a12_funs$cumhaz
  A13 <- a13_funs$cumhaz
  A23 <- a23_funs$cumhaz

  return_obj <- list(a12 = a12,
                     a13 = a13,
                     a23 = a23,
                     A12 = A12,
                     A13 = A13,
                     A23 = A23)

  class(return_obj) <- c(class(return_obj), "full_spline_hazard")
  return_obj
}

#### Likelihood functions for each case ####
case_1_likelihood <- function(V_0, V_healthy, V_ill, T_obs, a12, a13, a23, A12, A13, A23) {
  if (all(V_healthy == T_obs)) {
    exp(-A12(V_healthy) - A13(V_healthy)) / exp(-A12(V_0) - A13(V_0))
  } else {
    p1 <- 1 / exp(-A12(V_0) - A13(V_0))
    p2 <- exp(-A12(T_obs) - A13(T_obs))

    factored_out <- exp(-A23(T_obs))
#    grid <- sort(unique(c(V_0, V_healthy, T_obs)))
    grid <- seq(min(V_healthy), max(T_obs), length.out = 250)


    integrand <- exp(-A12(grid) - A13(grid)) * a12(grid) * exp(A23(grid))

    dx <- diff(grid)
    mid <- (integrand[-1] + integrand[-length(grid)]) / 2
    integral_at_grid <- c(0, cumsum(dx * mid))

    integral_fun <- splinefun(grid, integral_at_grid, method = "natural")

    p3 <- factored_out * (integral_fun(T_obs) - integral_fun(V_healthy))

    p1 * (p2 + p3)
  }
}
case_2_likelihood <- function(V_0, V_healthy, V_ill, T_obs, a12, a13, a23, A12, A13, A23) {
  p1 <- 1 / exp(-A12(V_0) - A13(V_0))
  p2 <- exp(-A12(T_obs) - A13(T_obs)) * a13(T_obs)

  factored_out <- exp(-A23(T_obs)) * a23(T_obs)
  #grid <- sort(unique(c(V_0, V_healthy, T_obs)))
  grid <- seq(min(V_healthy), max(T_obs), length.out = 250)
  integrand <- exp(-A12(grid) - A13(grid)) * a12(grid) * exp(A23(grid))

  dx <- diff(grid)
  mid <- (integrand[-1] + integrand[-length(grid)]) / 2
  integral_at_grid <- c(0, cumsum(dx * mid))

  integral_fun <- splinefun(grid, integral_at_grid, method = "natural")

  p3 <- factored_out * (integral_fun(T_obs) - integral_fun(V_healthy))

  p1 * (p2 + p3)
}
case_3_likelihood <- function(V_0, V_healthy, V_ill, T_obs, a12, a13, a23, A12, A13, A23) {
  p1 <- 1 / exp(-A12(V_0) - A13(V_0))

  factored_out <- exp(-A23(T_obs))
  #grid <- sort(unique(c(V_0, V_healthy, V_ill, T_obs)))
  grid <- seq(min(V_healthy), max(V_ill), length.out = 250)
  integrand <- exp(-A12(grid) - A13(grid)) * a12(grid) * exp(A23(grid))

  dx <- diff(grid)
  mid <- (integrand[-1] + integrand[-length(grid)]) / 2
  integral_at_grid <- c(0, cumsum(dx * mid))

  integral_fun <- splinefun(grid, integral_at_grid, method = "natural")



  p2 <- factored_out * (integral_fun(V_ill) - integral_fun(V_healthy))

  p1 * p2
}
case_4_likelihood <- function(V_0, V_healthy, V_ill, T_obs, a12, a13, a23, A12, A13, A23) {
  p1 <- 1 / exp(-A12(V_0) - A13(V_0))

  ####

  factored_out <- exp(-A23(T_obs)) * a23(T_obs)
  grid <- seq(min(V_healthy), max(V_ill), length.out = 250)
  integrand <- exp(-A12(grid) - A13(grid)) * a12(grid) * exp(A23(grid))

  dx <- diff(grid)
  mid <- (integrand[-1] + integrand[-length(grid)]) / 2
  integral_at_grid <- c(0, cumsum(dx * mid))

  integral_fun <- splinefun(grid, integral_at_grid, method = "natural")



  p2 <- factored_out * (integral_fun(V_ill) - integral_fun(V_healthy))

  p1 * p2
}


#### Penalizing ####

pen_mat_m_splines <- function(input_knots, degree = 3) {
  n_knots <- length(input_knots)
  d <- diff(input_knots)
  g_ab <- mSpline(x = input_knots,
                  knots = input_knots[-c(1,n_knots)],
                  Boundary.knots = range(input_knots),
                  degree = degree,
                  derivs = 2,
                  intercept = TRUE
  )
  knots_mid <- input_knots[-length(input_knots)] + d / 2
  g_ab_mid <- mSpline(x = knots_mid,
                      knots = input_knots[-c(1,n_knots)],
                      Boundary.knots = range(input_knots),
                      degree = degree,
                      derivs = 2,
                      intercept = TRUE
  )
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) +
      4 * crossprod(d * g_ab_mid, g_ab_mid) +
      crossprod(d * g_b, g_b)) / 6
}

#### Full log likelihood ####

# V_healthy = V_k or V_m
# V_ill = V_k+1


cal_log_likehood <- function(data, theta, degree, knots) {

  stopifnot(all(c("V_0", "V_healthy", "V_ill", "T_obs", "status_dead", "status_ill") %in% names(data)))

  V_0       = data$V_0
  V_healthy = data$V_healthy
  V_ill     = data$V_ill
  T_obs     = data$T_obs
  status_dead = data$status_dead
  status_ill  = data$status_ill

  # define all transitions
  a12_funs <- make_spline_funs(theta$a12, knots$a12, degree)
  a13_funs <- make_spline_funs(theta$a13, knots$a13, degree)
  a23_funs <- make_spline_funs(theta$a23, knots$a23, degree)

  # hazards
  a12 <- a12_funs$hazard
  a13 <- a13_funs$hazard
  a23 <- a23_funs$hazard

  # cumulative hazards
  A12 <- a12_funs$cumhaz
  A13 <- a13_funs$cumhaz
  A23 <- a23_funs$cumhaz


  status <- ifelse(status_ill == 0 & status_dead == 0, 1,
                    ifelse(status_ill == 0 & status_dead == 1, 2,
                        ifelse(status_ill == 1 & status_dead == 0, 3, 4)))

  ll_value <- 0
  if (any(as.integer(status) == 1)) ll_value <- ll_value + sum(log(case_1_likelihood(V_0[as.integer(status) == 1],
                                                                         V_healthy[as.integer(status) == 1],
                                                                         V_ill[as.integer(status) == 1],
                                                                         T_obs[as.integer(status) == 1],
                                                                         a12, a13, a23,
                                                                         A12, A13, A23)))

  if (any(as.integer(status) == 2)) ll_value <- ll_value + sum(log(case_2_likelihood(V_0[as.integer(status) == 2],
                                                                                   V_healthy[as.integer(status) == 2],
                                                                                   V_ill[as.integer(status) == 2],
                                                                                     T_obs[as.integer(status) == 2],
                                                                                     a12, a13, a23,
                                                                                     A12, A13, A23)))

  if (any(as.integer(status) == 3)) ll_value <- ll_value + sum(log(case_3_likelihood(V_0[as.integer(status) == 3],
                                                                                     V_healthy[as.integer(status) == 3],
                                                                                     V_ill[as.integer(status) == 3],
                                                                                     T_obs[as.integer(status) == 3],
                                                                                     a12, a13, a23,
                                                                                     A12, A13, A23)))

  if (any(as.integer(status) == 4)) ll_value <- ll_value + sum(log(case_4_likelihood(V_0[as.integer(status) == 4],
                                                                                     V_healthy[as.integer(status) == 4],
                                                                                     V_ill[as.integer(status) == 4],
                                                                                     T_obs[as.integer(status) == 4],
                                                                                     a12, a13, a23,
                                                                                     A12, A13, A23)))

  hazards = list(a12 = a12,
                 a13 = a13,
                 a23 = a23,
                 A12 = A12,
                 A13 = A13,
                 A23 = A23)

  class(hazards) <- c(class(hazards), "full_spline_hazard")

  return(list(ll_value = ll_value,
              hazards = hazards,
              theta = theta,
              knots = knots,
              degree = degree))
}


# Provide P_mat, list of P12, P13, P23, for speed
cal_pen_log_likehood <- function(data,
                                 theta,
                                 degree,
                                 knots,
                                 kappa_term,
                                 P_mat = NULL) {

  stopifnot(all(c("V_0", "V_healthy", "V_ill", "T_obs", "status_dead", "status_ill") %in% names(data)))

  stopifnot(length(kappa_term) == 3)

  ll <- cal_log_likehood(data, theta, degree, knots)

  ll_value <- ll$ll_value

  if(is.null(P_mat)) {
    P12 <- pen_mat_m_splines(knots$a12)
    P13 <- pen_mat_m_splines(knots$a13)
    P23 <- pen_mat_m_splines(knots$a23)
  } else {
    P12 <- P_mat$P12
    P13 <- P_mat$P13
    P23 <- P_mat$P23
  }

  raw_penalty <- c(drop(t(theta$a12) %*% P12 %*% theta$a12),
                   drop(t(theta$a13) %*% P13 %*% theta$a13),
                   drop(t(theta$a23) %*% P23 %*% theta$a23))

  penalty <- as.numeric(crossprod(raw_penalty, kappa_term))

  return(list(pl_value = ll_value - penalty,
              ll_value = ll_value,
              hazards = ll$hazards,
              theta = theta,
              knots = knots,
              penalty = penalty,
              degree = degree))
}




#### Fit with penalty using base::optim ####
do_likelihood_optim <- function(data,
                                kappa_term,
                                n_knots,
                                degree = 3,
                                knots = NULL,
                                long_theta_0 = NULL) {

  stopifnot(all(c("V_0", "V_healthy", "V_ill", "T_obs", "status_dead", "status_ill") %in% names(data)))


  n_par <- n_knots + degree - 1
  stopifnot(length(kappa_term) == 3)

  if (is.null(knots)) {
    knots_a12 <- seq(min(data$V_0),
                     max(data$T_obs),
                     length.out = n_knots)
    knots_a13 <- knots_a12
    knots_a23 <- knots_a12
    knots <- list(a12 = knots_a12, a13 = knots_a13, a23 = knots_a23)
  } else {
    stopifnot(length(knots$a12) == n_knots,
              length(knots$a13) == n_knots,
              length(knots$a23) == n_knots)
  }

  P_mat <- list(
    P12 = pen_mat_m_splines(knots$a12),
    P13 = pen_mat_m_splines(knots$a13),
    P23 = pen_mat_m_splines(knots$a23)
  )


  obj_fun <- function(long_theta) {
    theta <- list(
      a12 = long_theta[1:n_par],
      a13 = long_theta[(n_par+1):(2*n_par)],
      a23 = long_theta[(2*n_par+1):(3*n_par)]
    )

    pl <- cal_pen_log_likehood(data = data,
                               theta = theta,
                               degree = degree,
                               knots = knots,
                               kappa_term = kappa_term,
                               P_mat = P_mat)


    if (!is.finite(pl$pl_value)) return(1e10)
    -pl$pl_value
  }

  if(is.null(long_theta_0)) {
    long_theta_0 <- rep(0.5, 3*n_par)
  }

  stopifnot(length(long_theta_0) == 3*n_par)

  res <- optim(
    par   = long_theta_0,
    fn    = obj_fun,
    method= "L-BFGS-B",
    lower = 0
  )

  theta_hat <- list(
    a12 = res$par[1:n_par],
    a13 = res$par[(n_par+1):(2*n_par)],
    a23 = res$par[(2*n_par+1):(3*n_par)]
  )

  pl_at_theta_hat <- cal_pen_log_likehood(data = data,
                                          theta = theta_hat,
                                          degree = degree,
                                          knots = knots,
                                          kappa_term = kappa_term,
                                          P_mat = P_mat)


  res <- list(
    theta_hat = theta_hat,
    knots = knots,
    degree = degree,
    kappa_term = kappa_term,
    hazards = pl_at_theta_hat$hazards,
    ll_value = pl_at_theta_hat$ll_value,
    pl_value = pl_at_theta_hat$pl_value,
    data = data,
    P_mat = P_mat,
    n_par = n_par)

  class(res) <- c("pl_optim", class(res))
  res
}


do_likelihood_optim_unconstrained <- function(data,
                                kappa_term,
                                n_knots,
                                degree = 3,
                                knots = NULL,
                                long_theta_0 = NULL) {

  stopifnot(all(c("V_0", "V_healthy", "V_ill", "T_obs", "status_dead", "status_ill") %in% names(data)))


  n_par <- n_knots + degree - 1
  stopifnot(length(kappa_term) == 3)

  if (is.null(knots)) {
    knots_a12 <- seq(min(data$V_0),
                     max(data$T_obs),
                     length.out = n_knots)
    knots_a13 <- knots_a12
    knots_a23 <- knots_a12
    knots <- list(a12 = knots_a12, a13 = knots_a13, a23 = knots_a23)
  } else {
    stopifnot(length(knots$a12) == n_knots,
              length(knots$a13) == n_knots,
              length(knots$a23) == n_knots)
  }

  P_mat <- list(
    P12 = pen_mat_m_splines(knots$a12),
    P13 = pen_mat_m_splines(knots$a13),
    P23 = pen_mat_m_splines(knots$a23)
  )


  obj_fun <- function(long_theta) {

    long_theta <- exp(long_theta)

    theta <- list(
      a12 = long_theta[1:n_par],
      a13 = long_theta[(n_par+1):(2*n_par)],
      a23 = long_theta[(2*n_par+1):(3*n_par)]
    )

    pl <- cal_pen_log_likehood(data = data,
                               theta = theta,
                               degree = degree,
                               knots = knots,
                               kappa_term = kappa_term,
                               P_mat = P_mat)


    if (!is.finite(pl$pl_value)) return(1e10)
    -pl$pl_value
  }

  if(is.null(long_theta_0)) {
    long_theta_0 <- rep(-1, 3*n_par)
  }

  stopifnot(length(long_theta_0) == 3*n_par)

  res <- optim(
    par   = long_theta_0,
    fn    = obj_fun,
    method = "BFGS"
  )

  res$par <- exp(res$par)

  theta_hat <- list(
    a12 = res$par[1:n_par],
    a13 = res$par[(n_par+1):(2*n_par)],
    a23 = res$par[(2*n_par+1):(3*n_par)]
  )

  pl_at_theta_hat <- cal_pen_log_likehood(data = data,
                                          theta = theta_hat,
                                          degree = degree,
                                          knots = knots,
                                          kappa_term = kappa_term,
                                          P_mat = P_mat)


  res <- list(
    theta_hat = theta_hat,
    knots = knots,
    degree = degree,
    kappa_term = kappa_term,
    hazards = pl_at_theta_hat$hazards,
    ll_value = pl_at_theta_hat$ll_value,
    pl_value = pl_at_theta_hat$pl_value,
    data = data,
    P_mat = P_mat,
    n_par = n_par)

  class(res) <- c("pl_optim", class(res))
  res
}


safe_solve <- function(H_pl, H_ll, ridge_start = 0, tol = 1e-8, max_iter = 8) {
  # Ensure symmetry (Hessians should be symmetric, but numerics can drift)
  H <- (H_pl + t(H_pl)) / 2

  # Scale for the ridge so lambda is magnitude-aware
  sc <- mean(abs(diag(H)))
  if (!is.finite(sc) || sc == 0) sc <- 1

  # Try Cholesky with increasing ridge
  lambda <- ridge_start
  for (i in 0:max_iter) {
    H_reg <- H + (lambda * sc + .Machine$double.eps) * diag(nrow(H))
    R <- try(chol(H_reg), silent = TRUE)  # H_reg = t(R) %*% R
    if (!inherits(R, "try-error")) {
      # Solve H_reg X = H_ll via two triangular solves (no explicit inverse)
      Y <- forwardsolve(t(R), H_ll)   # t(R) Y = H_ll
      X <- backsolve(R, Y)            # R X = Y
      return(list(X = X, method = "chol", lambda = lambda * sc, iters = i))
    }
    lambda <- if (lambda == 0) tol else lambda * 10
  }

  # Fallback: SVD pseudo-inverse with thresholding
  sv <- svd(H)
  thr <- max(tol, max(sv$d) * 1e-12)
  d_inv <- ifelse(sv$d > thr, 1 / sv$d, 0)
  X <- sv$v %*% (d_inv * t(sv$u) %*% H_ll)
  list(X = X, method = "svd", lambda = lambda * sc, rank = sum(sv$d > thr))
}




approx_cv <- function(pl_optim) {

  stopifnot(inherits(pl_optim, "pl_optim"))

  theta_hat <- pl_optim$theta_hat
  ll_value <- pl_optim$ll_value
  pl_value <- pl_optim$pl_value
  data <- pl_optim$data
  degree <- pl_optim$degree
  knots <- pl_optim$knots
  kappa_term <- pl_optim$kappa_term
  P_mat <- pl_optim$P_mat
  data <- pl_optim$data
  n_par = pl_optim$n_par

  cal_ll_in_long_theta <- function(long_theta){

    theta <- list(
      a12 = long_theta[1:n_par],
      a13 = long_theta[(n_par+1):(2*n_par)],
      a23 = long_theta[(2*n_par+1):(3*n_par)]
    )

    cal_log_likehood(data = data, theta = theta, degree = degree, knots = knots)$ll_value
  }


  cal_pl_in_long_theta <- function(long_theta){

    theta <- list(
      a12 = long_theta[1:n_par],
      a13 = long_theta[(n_par+1):(2*n_par)],
      a23 = long_theta[(2*n_par+1):(3*n_par)]
    )

    cal_pen_log_likehood(data = data, theta = theta, degree = degree, knots = knots, kappa_term = kappa_term, P_mat = P_mat)$pl_value
  }
  long_theta_hat <- unlist(theta_hat)


  H_pl <- numDeriv::hessian(cal_pl_in_long_theta, long_theta_hat)
  H_ll  <- numDeriv::hessian(cal_ll_in_long_theta, long_theta_hat)

  res <- safe_solve(H_pl, H_ll) # solves H_pl X = H_ll
  tr_val <- sum(diag(res$X))


  approx_cv <- ll_value - tr_val

  pl_optim$approx_cv_value <- approx_cv
  class(pl_optim) <- c(class(pl_optim), "pl_optim_w_cv")
  pl_optim
}

fit_idm <- function(data,
                    n_knots = 7,
                    degree  = 3,
                    knots   = NULL,
                    kappa_values = 10^(-2:2),
                    kappa_values_full = NULL,
                    verbose = TRUE) {
  if (is.null(kappa_values_full)) {
    kappa_values <- expand.grid(kappa_values, kappa_values, kappa_values, KEEP.OUT.ATTRS = FALSE)
  } else {
    kappa_values <- as.data.frame(kappa_values_full, stringsAsFactors = FALSE)
  }

  n <- nrow(kappa_values)
  cvs  <- numeric(n)
  fits <- vector("list", n)

  for (i in seq_len(n)) {
    kappa_term <- as.numeric(kappa_values[i, ])

    fit <- do_likelihood_optim(
      data,
      n_knots    = n_knots,
      kappa_term = kappa_term,
      knots      = knots,
      degree     = degree
    )

    cv_i <- approx_cv(fit)$approx_cv_value

    fits[[i]] <- fit
    cvs[i]    <- cv_i

    if (verbose) {
      kappa_str <- paste(formatC(kappa_term, format = "g"), collapse = ", ")
      message(sprintf("[%d/%d] kappa=(%s)  cv=%g", i, n, kappa_str, cv_i))
    }
  }

  best <- which.min(cvs)
  list(
    fit   = fits[[best]],
    kappa = as.numeric(kappa_values[best, ]),
    cv    = cvs[best],
    full_data = list(cvs = cvs,
                     fits = fits)
  )
}



fit_idm_greedy <- function(data,
                           n_knots = 7,
                           degree  = 3,
                           knots   = NULL,
                           kappa_values = NULL,
                           fracs = c(0.25, 0.5, 1.0),
                           drop_rate = 0.5,
                           verbose = TRUE,
                           seed = 1) {
  if (is.null(kappa_values)) {
    kv <- 10^(-2:2)
    kappa_values <- expand.grid(kv, kv, kv, KEEP.OUT.ATTRS = FALSE)
  } else {
    kappa_values <- expand.grid(kappa_values, kappa_values, kappa_values, KEEP.OUT.ATTRS = FALSE)
  }

  set.seed(seed)
  n_cand <- nrow(kappa_values)
  alive  <- seq_len(n_cand)
  cvs    <- rep(Inf, n_cand)
  fits   <- vector("list", n_cand)

  N <- nrow(data)

  for (r in seq_along(fracs)) {
    m <- max(1, floor(fracs[r] * N))
    idx <- if (m < N) sample.int(N, m) else seq_len(N)
    obs_r <- data[idx, , drop = FALSE]

    for (i in alive) {
      kappa_term <- as.numeric(kappa_values[i, ])
      fit <- do_likelihood_optim(
        obs_r,
        n_knots    = n_knots,
        kappa_term = kappa_term,
        knots      = knots,
        degree     = degree
      )
      cv_i <- approx_cv(fit)$approx_cv_value
      fits[[i]] <- fit
      cvs[i]    <- cv_i

      if (verbose) {
        kappa_str <- paste(formatC(kappa_term, format = "g"), collapse = ", ")
        message(sprintf("[round %d] i=%d/%d  m=%d  kappa=(%s)  cv=%g",
                        r, i, n_cand, m, kappa_str, cv_i))
      }
    }

    if (length(alive) > 1 && r < length(fracs)) {
      keep_n <- max(1, min(length(alive) - 1L,
                           floor((1 - drop_rate) * length(alive))))
      alive  <- alive[order(cvs[alive])][seq_len(keep_n)]
    }
  }

  best <- which.min(cvs)
  list(
    fit   = fits[[best]],
    kappa = as.numeric(kappa_values[best, ]),
    cv    = cvs[best],
    full_data = list(cvs = cvs, fits = fits)
  )
}




