# create a plotting method for the idm_hazards class using ggplot2
#' @export
plot.idm_hazards <- function(x, max_time = 1, cumulative = FALSE, add = NULL, label = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for plotting.")
  }

  # create a data frame for plotting
  time_seq <- seq(0, max_time, length.out = 256)
  if (cumulative) {
    hazard_df <- data.frame(
      time = time_seq,
      A12 = x$A12(time_seq),
      A13 = x$A13(time_seq),
      A23 = x$A23(time_seq)
    )
    hazard_df_long <- tidyr::pivot_longer(hazard_df, cols = c("A12", "A13", "A23"),
                                            names_to = "hazard", values_to = "value")
  } else  {
    if(!"a12" %in% names(x) || !"a13" %in% names(x) || !"a23" %in% names(x)) {
      stop("Hazard functions a12, a13, and a23 must be available in the idm_hazards object to plot non-cumulative hazards.")
    }
    hazard_df <- data.frame(
      time = time_seq,
      a12 = x$a12(time_seq),
      a13 = x$a13(time_seq),
      a23 = x$a23(time_seq)
    )
    hazard_df_long <- tidyr::pivot_longer(hazard_df, cols = c("a12", "a13", "a23"),
                                            names_to = "hazard", values_to = "value")
  }

  # Add estimator label if provided
  if (!is.null(label)) {
    hazard_df_long$estimator <- label
  } else {
    hazard_df_long$estimator <- "Estimator 1"
  }

  # If adding to existing plot
  if (!is.null(add) && inherits(add, "ggplot")) {
    # Extract existing data and combine
    p <- add +
      ggplot2::geom_line(data = hazard_df_long,
                        ggplot2::aes(x = time, y = value, color = hazard, linetype = estimator))
    return(p)
  }

  # Create new plot and store parameters as attributes
  p <- ggplot2::ggplot(hazard_df_long, ggplot2::aes(x = time, y = value, color = hazard, linetype = estimator)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = if(cumulative) "Cumulative Hazard Functions" else "Hazard Functions",
                  x = "Time",
                  y = if(cumulative) "Cumulative Hazard" else "Hazard Rate") +
    ggplot2::theme_minimal()

  # Store plot parameters as attributes for later use
  attr(p, "idm_max_time") <- max_time
  attr(p, "idm_cumulative") <- cumulative

  return(p)
}


# We want to be able to add plots together using ggplot2's + operator
#' @export
`+.idm_hazards` <- function(e1, e2) {
  if (inherits(e1, "idm_hazards") && inherits(e2, "ggplot")) {
    # idm_hazards + ggplot: create plot from idm_hazards then add ggplot layer
    return(plot(e1, label = deparse(substitute(e1))) + e2)
  } else if (inherits(e1, "idm_hazards") && inherits(e2, "idm_hazards")) {
    # idm_hazards + idm_hazards: create base plot and add second
    p <- plot(e1, label = deparse(substitute(e1)))
    p <- plot(e2, add = p, label = deparse(substitute(e2)))
    return(p)
  } else {
    stop("Can only combine idm_hazards objects with ggplot objects or other idm_hazards objects.")
  }
}

# Register a method for ggplot + idm_hazards
#' @export
ggplot_add.idm_hazards <- function(object, plot, object_name) {
  # When ggplot + idm_hazards, add the idm_hazards data to the plot
  # Extract parameters from the existing plot
  max_time <- attr(plot, "idm_max_time")
  cumulative <- attr(plot, "idm_cumulative")

  # Use default values if attributes not found
  if (is.null(max_time)) max_time <- 100
  if (is.null(cumulative)) cumulative <- FALSE

  # Add the new idm_hazards to the existing plot with same parameters
  result <- plot(object, max_time = max_time, cumulative = cumulative,
                 add = plot, label = object_name)

  # Preserve attributes
  attr(result, "idm_max_time") <- max_time
  attr(result, "idm_cumulative") <- cumulative

  return(result)
}
