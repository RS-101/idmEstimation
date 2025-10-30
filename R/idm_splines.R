library(splines2)
#### Create spline hazard ####
make_spline_mat <- function(x, knots, degree = 3) {
  n_knots <- length(knots)
  i_spline_mat <- splines2::iSpline(x,
    degree = degree,
    knots = knots[-c(1, n_knots)],
    Boundary.knots = knots[c(1, n_knots)],
    intercept = TRUE,
    warn.outside = FALSE
  )

  m_spline_mat <- splines2::mSpline(x,
    degree = degree,
    knots = knots[-c(1, n_knots)],
    Boundary.knots = knots[c(1, n_knots)],
    intercept = TRUE,
    warn.outside = FALSE
  )

  spline_mat_list <- list(
    i_spline = i_spline_mat,
    m_spline = m_spline_mat
  )

  class(spline_mat_list) <- c("spline_mat_list", class(spline_mat_list))
  spline_mat_list
}


