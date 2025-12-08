#' Plot Illness-Death Model Estimates
#'
#' Creates visualizations of estimated hazard functions, cumulative hazards, and
#' distribution functions from fitted illness-death models. Can overlay multiple
#' model fits for comparison.
#'
#' @param object An object of class \code{"idm_object"} from \code{\link{fit_npmle}},
#'   \code{\link{fit_pc_model}}, or \code{\link{fit_spline_model}}.
#' @param ... Additional \code{idm_object}s to overlay on the same plot, or other
#'   arguments.
#' @param time_points Integer. Number of time points for evaluating functions.
#'   Default is 100.
#'
#' @return A \code{ggplot} object or a list of \code{ggplot} objects (one for each
#'   function type: hazards, cumulative hazards, distributions). Use with
#'   \code{\link[patchwork]{plot_layout}} to arrange multiple plots.
#'
#' @details
#' The function evaluates estimated functions on a grid spanning the observed data
#' range and creates line plots using \code{\link[ggplot2]{ggplot}}. Different
#' models are distinguished by line type, and different transitions by color.
#'
#' For NPMLE estimates, functions may return \code{NA} in indeterminate regions
#' (Turnbull intervals), which appear as gaps in the plots.
#'
#' @examples
#' \donttest{
#' # Simulate data and fit multiple models
#' set.seed(333)
#' sim_data <- simulate_idm_constant_hazards(n = 250)
#'
#' fit_pc <- fit_pc_model(sim_data$data, n_knots = 5)
#' fit_npm <- fit_npmle(sim_data$data, max_iter = 50)
#'
#' # Plot both fits together
#' plot(fit_pc, fit_npm)
#'
#' # Plot only one model
#' plot(fit_pc)
#' }
#'
#' @seealso \code{\link{fit_npmle}}, \code{\link{fit_pc_model}}, \code{\link{fit_spline_model}}
#' @export
plot_fit <- function(idm_fit, ..., time_points = 100) {
  library(ggplot2)

  # Separate idm_fit objects from function arguments
  dots <- list(...)
  is_idm_fit <- sapply(dots, function(x) inherits(x, "idm_object"))

  # Collect all idm_fit objects
  all_fits <- c(list(idm_fit), dots[is_idm_fit])

  # Capture additional arguments for functions (non-idm_fit objects)
  func_args <- dots[!is_idm_fit]

  # Determine global max_time across all fits
  max_time <- max(sapply(all_fits, function(fit) max(fit$data$T_obs)))
  time_grid <- seq(0, max_time, length.out = time_points)

  # Initialize data frames for each plot type
  hazard_data <- data.frame()
  cum_hazard_data <- data.frame()
  dist_data <- data.frame()

  # Extract data from each fit
  for (i in seq_along(all_fits)) {
    fit <- all_fits[[i]]
    model_name <- fit$model_type

    # Process hazard functions if they exist
    if (!is.null(fit$estimators$hazard_functions)) {
      for (func_name in names(fit$estimators$hazard_functions)) {
        func <- fit$estimators$hazard_functions[[func_name]]
        values <- sapply(time_grid, function(t) {
          tryCatch(
            do.call(func, c(list(t), func_args)),
            error = function(e) func(t)
          )
        })
        hazard_data <- rbind(hazard_data, data.frame(
          time = time_grid,
          value = values,
          transition = func_name,
          estimator = model_name
        ))
      }
    }

    # Process cumulative hazard functions
    if (!is.null(fit$estimators$cum_hazard_functions)) {
      for (func_name in names(fit$estimators$cum_hazard_functions)) {
        func <- fit$estimators$cum_hazard_functions[[func_name]]
        values <- sapply(time_grid, function(t) {
          tryCatch(
            do.call(func, c(list(t), func_args)),
            error = function(e) func(t)
          )
        })
        cum_hazard_data <- rbind(cum_hazard_data, data.frame(
          time = time_grid,
          value = values,
          transition = func_name,
          estimator = model_name
        ))
      }
    }

    # Process distribution functions
    if (!is.null(fit$estimators$distribution_functions)) {
      for (func_name in names(fit$estimators$distribution_functions)) {
        func <- fit$estimators$distribution_functions[[func_name]]
        values <- sapply(time_grid, function(t) {
          tryCatch(
            do.call(func, c(list(t), func_args)),
            error = function(e) func(t)
          )
        })
        dist_data <- rbind(dist_data, data.frame(
          time = time_grid,
          value = values,
          transition = func_name,
          estimator = model_name
        ))
      }
    }
  }

  # Create plots
  plots <- list()

  if (nrow(hazard_data) > 0) {
    plots$hazards <- ggplot(hazard_data, aes(x = time, y = value, color = transition, linetype = estimator)) +
      geom_line(linewidth = 1) +
      labs(title = "Hazard Functions", x = "Time", y = "Hazard", color = "Transition", linetype = "Estimator") +
      theme_minimal()
  }

  if (nrow(cum_hazard_data) > 0) {
    plots$cumulative_hazards <- ggplot(cum_hazard_data, aes(x = time, y = value, color = transition, linetype = estimator)) +
      geom_line(linewidth = 1) +
      labs(title = "Cumulative Hazard Functions", x = "Time", y = "Cumulative Hazard", color = "Transition", linetype = "Estimator") +
      theme_minimal()
  }

  if (nrow(dist_data) > 0) {
    plots$distributions <- ggplot(dist_data, aes(x = time, y = value, color = transition, linetype = estimator)) +
      geom_line(linewidth = 1) +
      labs(title = "Distribution Functions", x = "Time", y = "Probability", color = "Transition", linetype = "Estimator") +
      theme_minimal()
  }

  # Return single plot or list of plots
  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    return(plots)
  }
}
