max_pc_likelihood <- function(cpp_pointer, model_config, long_theta_0 = NULL) {
  n_theta_12 <- model_config$n_theta_12
  n_theta_13 <- model_config$n_theta_13
  n_theta_23 <- model_config$n_theta_23

  n_theta <- n_theta_12 + n_theta_13 + n_theta_23

  obj_fun <- function(long_theta) {

    res <- calc_log_likelihood(
      cpp_pointer,
      long_theta[1:n_theta_12],
      long_theta[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
      long_theta[(n_theta_12 + n_theta_13 + 1):n_theta]
    )

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
    lower = 0)

  theta_hat_list <- list(
    theta_12 = res$par[1:n_theta_12],
    theta_13 = res$par[(n_theta_12 + 1):(n_theta_12 + n_theta_13)],
    theta_23 = res$par[(n_theta_12 + n_theta_13 + 1):n_theta]
  )

  list(
    theta_hat_list = theta_hat_list,
    optim_res = res
  )
}

fit_pc_model <- function(data,
                    knots_12 = NULL,
                    knots_13 = NULL,
                    knots_23 = NULL,
                    n_knots = 7,
                    use_bSpline = FALSE) {


  data_object <- create_case_data(data)
  summary_data_object <- summarise_data_object(data_object)

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
      0,
      summary_data_object$last_possible_23,
      length.out = n_knots)
  }


  model_config <- list(
    knots_12 = knots_12,
    knots_13 = knots_13,
    knots_23 = knots_23,
    n_theta_12 = length(knots_12) - 1,
    n_theta_13 = length(knots_13) - 1,
    n_theta_23 = length(knots_23) - 1,
    n_knots = n_knots,
    degree = 0,
    use_bSpline = use_bSpline
  )

  # Setup model configuration
  cpp_pointer <- setup_cpp_model(data_object, model_config)

  fit <- max_pc_likelihood(cpp_pointer, model_config)

  estimators <- create_estimators(model_config, fit$theta_hat_list)



  # idm fit object
  idm_fit <- list(
    estimators = estimators,
    data = data,
    model_type = "piecewise_constant",
    raw_estimates = fit$theta_hat_list,
    model_config = model_config,
    converged = fit$optim_res$convergence == 0,
    model_specific = list(
      optim_res = fit$optim_res
    )
  )
  class(idm_fit) <- c("idm_fit", class(idm_fit))

  return(idm_fit)
}
