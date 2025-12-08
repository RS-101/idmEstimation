calculate_penalty_matrix <- function(knots, degree = 3) {
  n_knots <- length(knots)
  d <- diff(knots)
  g_ab <- splines2::mSpline(x = knots,
                  knots = knots[-c(1, n_knots)],
                  Boundary.knots = range(knots),
                  degree = degree,
                  derivs = 2,
                  intercept = TRUE
  )
  knots_mid <- knots[-length(knots)] + d / 2
  g_ab_mid <- splines2::mSpline(x = knots_mid,
                      knots = knots[-c(1, n_knots)],
                      Boundary.knots = range(knots),
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



max_pen_likelihood <- function(
  cpp_pointer,
  cpp_pointer_penalty,
  model_config,
  kappa_12,
  kappa_13,
  kappa_23,
  long_theta_0 = NULL,
  ...) {
  n_theta_12 <- model_config$n_theta_12
  n_theta_13 <- model_config$n_theta_13
  n_theta_23 <- model_config$n_theta_23

  n_theta <- n_theta_12 + n_theta_13 + n_theta_23

  obj_fun <- function(long_theta) {

    res <- calc_penalized_log_likelihood(
      cpp_pointer,
      cpp_pointer_penalty,
      long_theta[1:n_theta_12],
      long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
      long_theta[(n_theta_12 + n_theta_13 + 1):n_theta],
      kappa_12,
      kappa_13,
      kappa_23
    )$penalized_log_likelihood

    if (!is.finite(res)) return(1e10)
    -res
  }


  if(is.null(long_theta_0)) {
    long_theta_0 <- rep(0.5, n_theta)
  }

  res <- optim(
    par = long_theta_0,
    fn = obj_fun,
    method = "L-BFGS-B",
    lower = 0,
    ...)

  theta_hat_list <- list(
    theta_12 = res$par[1:n_theta_12],
    theta_13 = res$par[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
    theta_23 = res$par[(n_theta_12 + n_theta_13 + 1):n_theta]
  )

  pl_at_theta_hat <- calc_penalized_log_likelihood(
    cpp_pointer,
    cpp_pointer_penalty,
    theta_hat_list$theta_12,
    theta_hat_list$theta_13,
    theta_hat_list$theta_23,
    kappa_12,
    kappa_13,
    kappa_23
  )

  list(
    theta_hat_list = theta_hat_list,
    kappa_12 = kappa_12,
    kappa_13 = kappa_13,
    kappa_23 = kappa_23,
    log_likelihood = max(pl_at_theta_hat$log_likelihood, -1e10),
    penalized_log_likelihood = max(pl_at_theta_hat$penalized_log_likelihood, -1e10),
    penalty = pl_at_theta_hat$penalty,
    optim_res = res
  )
}


approx_cv <- function(model_config, fit) {

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

  theta_hat_list <- fit$theta_hat_list
  ll_value <- fit$log_likelihood
  pl_value <- fit$penalized_log_likelihood

  if(abs(ll_value) > 1e9 | abs(pl_value) > 1e9){
    return(-1e9)
  }

  model_pointer <- model_config$model_pointer
  penalty_pointer <- model_config$penalty_pointer
  n_theta_12 <- model_config$n_theta_12
  n_theta_13 <- model_config$n_theta_13
  n_theta_23 <- model_config$n_theta_23
  n_theta <- n_theta_12 + n_theta_13 + n_theta_23


  ll_in_long_theta <- function(long_theta){
    calc_log_likelihood(
      model_pointer,
      theta_12 = long_theta[1:n_theta_12],
      theta_13 = long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
      theta_23 = long_theta[(n_theta_12 + n_theta_13 + 1):n_theta])
  }

  pl_in_long_theta <- function(long_theta){
    calc_penalized_log_likelihood(
      model_pointer,
      penalty_pointer,
      theta_12 = long_theta[1:n_theta_12],
      theta_13 = long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
      theta_23 = long_theta[(n_theta_12 + n_theta_13 + 1):n_theta],
      kappa_12 = fit$kappa_12,
      kappa_13 = fit$kappa_13,
      kappa_23 = fit$kappa_23)$penalized_log_likelihood
  }


  long_theta_hat <- unlist(fit$theta_hat_list)

  H_pl <- numDeriv::hessian(pl_in_long_theta, long_theta_hat)
  H_ll  <- numDeriv::hessian(ll_in_long_theta, long_theta_hat)

  res <- safe_solve(H_pl, H_ll) # solves H_pl X = H_ll
  tr_val <- sum(diag(res$X))

  approx_cv <- ll_value - tr_val
  approx_cv
}

fit_spline_model <- function(data,
                    knots_12 = NULL,
                    knots_13 = NULL,
                    knots_23 = NULL,
                    degree = 3,
                    n_knots = 7,
                    kappa_12 = NULL,
                    kappa_13 = NULL,
                    kappa_23 = NULL,
                    verbose = TRUE,
                    run_in_parallel = FALSE,
                    ...) {

  data_object <- create_case_data(data)
  summary_data_object <- summarise_data_object(data_object)

  # Set default knots if not provided
  if (is.null(knots_12)) {
    knots_12 <- seq(
      0,
      summary_data_object$last_possible_12,
      length.out = n_knots)
  }
  if (is.null(knots_13)) {
    knots_13 <- seq(
      0,
      summary_data_object$last_possible_13,
      length.out = n_knots)
  }
  if (is.null(knots_23)) {
    knots_23 <- seq(
      summary_data_object$first_possible_23,
      summary_data_object$last_possible_23, length.out = n_knots)
  }


  model_config <- list(
    knots_12 = knots_12,
    knots_13 = knots_13,
    knots_23 = knots_23,
    n_theta_12 = length(knots_12) + degree - 1,
    n_theta_13 = length(knots_13) + degree - 1,
    n_theta_23 = length(knots_23) + degree - 1,
    degree = degree,
    n_knots = n_knots
  )

  penalty_config <- list(
    penalty_matrix_12 = calculate_penalty_matrix(knots_12, degree),
    penalty_matrix_13 = calculate_penalty_matrix(knots_13, degree),
    penalty_matrix_23 = calculate_penalty_matrix(knots_23, degree),
    kappa_12 = if (!is.null(kappa_12)) kappa_12 else 10^seq(-2, 20, length.out = 6),
    kappa_13 = if (!is.null(kappa_13)) kappa_13 else 10^seq(-2, 20, length.out = 6),
    kappa_23 = if (!is.null(kappa_23)) kappa_23 else 10^seq(-2, 20, length.out = 6)
  )



  # Setup model configuration
  cpp_pointer <- setup_cpp_model(data_object, model_config)

  cpp_pointer_penalty <- create_penalty_config(penalty_config)

  model_config$model_pointer <- cpp_pointer
  model_config$penalty_pointer <- cpp_pointer_penalty

  # make kappa grid
  kappa_grid <- expand.grid(
    kappa_12 = penalty_config$kappa_12,
    kappa_13 = penalty_config$kappa_13,
    kappa_23 = penalty_config$kappa_23
  )

  kappa_iterations <- nrow(kappa_grid)

  cvs <- numeric(kappa_iterations)
  fits <- vector("list", kappa_iterations)

  if(run_in_parallel){
    res <- parallel::mclapply(seq_len(kappa_iterations), function(i) {
      kappa_12 <- kappa_grid$kappa_12[i]
      kappa_13 <- kappa_grid$kappa_13[i]
      kappa_23 <- kappa_grid$kappa_23[i]

      fit <- max_pen_likelihood(cpp_pointer, cpp_pointer_penalty, model_config, kappa_12, kappa_13, kappa_23,...)
      approx_cv_value <- approx_cv(model_config, fit)

      list(fit = fit, cv = approx_cv_value)
    }, mc.cores = max(1, parallel::detectCores() - 1))

    fits <- lapply(res, `[[`, "fit")
    cvs  <- vapply(res, function(x) x$cv, numeric(1))


  } else {
    for (i in seq_len(kappa_iterations)) {
      kappa_12 <- kappa_grid$kappa_12[i]
      kappa_13 <- kappa_grid$kappa_13[i]
      kappa_23 <- kappa_grid$kappa_23[i]

      fit <- max_pen_likelihood(cpp_pointer, cpp_pointer_penalty, model_config, kappa_12, kappa_13, kappa_23,...)
      approx_cv_value <- approx_cv(model_config, fit)

      fits[[i]] <- fit
      cvs[i] <- approx_cv_value

      if (verbose) {
        message(sprintf(
          "[%d/%d] kappa=(%.2g, %.2g, %.2g)  cv=%.4f",
          i, kappa_iterations, kappa_12, kappa_13, kappa_23, cvs[i]
        ))
      }
    }
  }

  best <- which.max(cvs)
  final_kappas <- c(
    kappa_12 = kappa_grid$kappa_12[best],
    kappa_13 = kappa_grid$kappa_13[best],
    kappa_23 = kappa_grid$kappa_23[best]
  )
  final_fit <- fits[[best]]

  estimators <- create_estimators(model_config, final_fit$theta_hat_list)


    # idm fit object
  idm_fit <- list(
    estimators = estimators,
    data = data,
    model_type = "penalized_spline",
    raw_estimates = final_fit$theta_hat_list,
    model_config = model_config,
    converged = final_fit$optim_res$convergence == 0,
    model_specific = list(
      final_kappas = final_kappas,
      optim_res = final_fit$optim_res,
      penalty_config = penalty_config,
      full_cv_results = list(
            cv_values = cvs,
            fits = fits
          )
    )
  )
  class(idm_fit) <- c("idm_object", class(idm_fit))

  return(idm_fit)
}
