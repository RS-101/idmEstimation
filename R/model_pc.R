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

#' Fit Piecewise-Constant Hazard Model for Illness-Death Data
#'
#' Estimates transition hazards using a piecewise-constant (histogram) model.
#' Hazards are constant within intervals defined by knots and optimized via
#' maximum likelihood.
#'
#' @param data Data frame with observed illness-death data (see \code{\link{fit_npmle}}
#'   for required columns).
#' @param knots_12 Numeric vector of knot positions for the 1→2 transition.
#'   If \code{NULL}, \code{n_knots} equally-spaced knots are created. Default is \code{NULL}.
#' @param knots_13 Numeric vector of knot positions for the 1→3 transition.
#'   Default is \code{NULL}.
#' @param knots_23 Numeric vector of knot positions for the 2→3 transition.
#'   Default is \code{NULL}.
#' @param n_knots Integer. Number of knots for each transition if knots are not
#'   specified. Default is 7.
#' @param use_bSpline Logical. If \code{TRUE}, uses B-spline basis (degree 0 for
#'   piecewise constant). If \code{FALSE}, uses I-spline basis. Default is \code{FALSE}.
#'
#' @return An object of class \code{"idm_object"} containing:
#'   \describe{
#'     \item{estimators}{List with hazard, cumulative hazard, and distribution functions}
#'     \item{data}{Original input data}
#'     \item{model_type}{Character: \code{"piecewise_constant"}}
#'     \item{raw_estimates}{List of theta parameters for each transition}
#'     \item{model_config}{Knot positions and model configuration}
#'     \item{converged}{Logical from \code{\link[stats]{optim}}}
#'     \item{model_specific}{Optimization details}
#'   }
#'
#' @details
#' The piecewise-constant model assumes hazards are constant within intervals
#' \eqn{[t_{j}, t_{j+1})} defined by the knots. This is a flexible semi-parametric
#' approach that can approximate complex hazard shapes.
#' 
#' The model uses \code{\link[stats]{optim}} with L-BFGS-B to maximize the likelihood
#' subject to non-negativity constraints on the hazard parameters.
#' 
#' Knot placement affects model fit: too few knots may underfit, while too many
#' may overfit. Use \code{\link{fit_spline_model}} with cross-validation for
#' automatic complexity selection.
#'
#' @examples
#' # Simulate data
#' set.seed(321)
#' sim_data <- simulate_idm_weibull(n = 300, shape12 = 2, shape13 = 3, shape23 = 1.5)
#' 
#' # Fit piecewise-constant model
#' fit_pc <- fit_pc_model(sim_data$data, n_knots = 6)
#' 
#' # Evaluate estimated hazard at specific times
#' time_points <- c(1, 5, 10, 15)
#' fit_pc$estimators$hazard_functions$a12(time_points)
#' 
#' # Plot results
#' plot(fit_pc)
#'
#' @seealso \code{\link{fit_spline_model}}, \code{\link{fit_npmle}}
#' @export
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
  class(idm_fit) <- c("idm_object", class(idm_fit))

  return(idm_fit)
}
