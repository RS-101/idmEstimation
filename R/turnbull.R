#' Fit Turnbull NPMLE for Interval-Censored Data
#'
#' Computes the Turnbull non-parametric maximum likelihood estimator for
#' interval-censored observations using an EM-like algorithm.
#'
#' @param obs Data frame with interval observations containing \code{L}, \code{R},
#'   \code{L_tilde}, \code{R_tilde}.
#' @param use_EM Logical. If \code{TRUE}, uses EM updates; if \code{FALSE}, uses
#'   an alternative update scheme. Default is \code{TRUE}.
#' @param s Numeric vector of initial mass values. If \code{NULL}, initialized randomly.
#'
#' @return List with:
#'   \describe{
#'     \item{s}{Vector of probability masses on equivalence classes}
#'     \item{Q}{Matrix of equivalence class intervals (support points)}
#'   }
#'
#' @references
#' Turnbull, B. W. (1976). The empirical distribution function with arbitrarily
#' grouped, censored and truncated data. \emph{Journal of the Royal Statistical
#' Society: Series B}, 38(3), 290-295.
fit_turnbull <- function(obs, use_EM = T, s = NULL) {
  make_Q <- function() {
    LL <- sort(unique(union(obs$L, obs$R_tilde)))
    LL <- LL[is.finite(LL)]
    RR <- sort(unique(union(obs$R, obs$L_tilde)))

    Q <- matrix(c(rep(0L, length(LL)), rep(1L, length(RR)),
                  LL, RR), ncol = 2)

    Q <- Q[order(Q[,2]), ]
    tag <- which(diff(Q[, 1], 1) == 1)
    Q <- matrix(c(Q[tag, 2], Q[tag + 1, 2]), ncol = 2)

    Q
  }
  Q <- make_Q()

  c_ij <- outer(obs$L,Q[,1], "<=") & outer(obs$R, Q[,2], ">=")
  d_ij <- outer(obs$L_tilde,Q[,1], "<=") & outer(obs$R_tilde, Q[,2], ">=")

  # print(d_ij)


  ll <- function(s) sum(log(colSums(s*t(c_ij)))) -  sum(log(colSums(s*t(d_ij))))

  d_ll <- Vectorize(function(s, index) {
    C <- colSums(s*t(c_ij))
    D <- colSums(s*t(d_ij))

    sum(c_ij[,index]/C - d_ij[,index]/D)
  }, vectorize.args = "index")

  M <- function(s) sum(1/colSums(s*t(d_ij)))

  m <- nrow(Q)
  if(is.null(s)) s <- runif(m)
  s <- s/sum(s)
  not_conv <- TRUE

  print(M(s))
  while(not_conv){
    if(M(s) > 0) {
      if (use_EM) {
        s_new <- s*(d_ll(s,1:m)+M(s))/(sum(s*(d_ll(s,1:m)+M(s))))
      } else{
        s_new <- s*(1+d_ll(s,1:m)/M(s))
      }
    } else {
      s_new = rep(0,m)
    }

    not_conv = sum(abs(s_new-s)) > 1e-05
    s = s_new
    print(s)
  }

  print(d_ll(s,1:m))
  if( any(!(s > 0 & abs(d_ll(s,1:m)) < 1e-4) | (s = 0 & d_ll(s,1:m) < 1e-4)) )
    warning("KKT NOT STATISFIED")

  list(s = s, Q = Q)
}

plot_turnbull <- function(s, Q,
                          obs_intervals = NULL,        # matrix/data.frame with columns L,R
                          show_extrema = TRUE,
                          box_border = "black",
                          box_fill = rgb(0, 0, 0, maxColorValue = 255, alpha = 40),
                          step_col  = "black",
                          lower_col = "gray40",
                          upper_col = "gray40",
                          obs_col   = "steelblue",
                          obs_lwd   = 2,
                          obs_y     = -0.06,          # vertical level for observation line
                          obs_point_as_tick = TRUE,   # if L==R, draw a small tick
                          lwd = 2,
                          xlab = "x", ylab = expression(hat(F)(x)),
                          main = "Turnbull estimator with uncertainty boxes",
                          xlim = NULL, ylim = c(0, 1)) {

  if (!is.numeric(s) || !is.matrix(Q) || ncol(Q) != 2)
    stop("Provide numeric s and a two-column matrix Q = [p, q].")
  if (length(s) != nrow(Q))
    stop("Length(s) must equal nrow(Q).")
  if (any(Q[,1] > Q[,2]))
    stop("Each row of Q must satisfy p <= q.")
  if (abs(sum(s) - 1) > 1e-8)
    warning("sum(s) != 1; plotting proceeds with given masses.")

  # Order classes by increasing q, then p (common convention)
  o <- order(Q[,2], Q[,1])
  s <- s[o]
  Q <- Q[o, , drop = FALSE]
  m <- length(s)
  p <- Q[,1]; q <- Q[,2]
  cs <- cumsum(s)
  cs_prev <- c(0, head(cs, -1))

  if (is.null(xlim)) {
    xlim <- range(Q[,1])
    xlim[2] <- xlim[2] + 1
    if (!is.null(obs_intervals)) {
      O <- as.matrix(obs_intervals)
      xlim <- range(c(xlim, O), finite = TRUE)
    }
  }

  # Extend ylim to make room for the observation line (if provided)
  if (!is.null(obs_intervals)) {
    ylim <- c(min(ylim[1], obs_y - 0.02), ylim[2])
  }

  # Frame
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, axes = FALSE)
  axis(1); axis(2); box()

  # Constant parts between equivalence classes (determinate regions)
  segments(x0 = min(xlim), y0 = 0, x1 = q[1], y1 = 0, lwd = lwd, col = step_col)
  if (m >= 2) {
    for (j in 1:(m-1)) {
      if (q[j] < p[j+1]) {
        segments(x0 = q[j], y0 = cs[j], x1 = p[j+1], y1 = cs[j], lwd = lwd, col = step_col)
      }
    }
  }
  segments(x0 = p[m], y0 = 1, x1 = max(xlim), y1 = 1, lwd = lwd, col = step_col)

  # Uncertainty boxes over each class [p_j, q_j]
  for (j in 1:m) {
    rect(xleft = p[j], ybottom = cs_prev[j], xright = q[j], ytop = cs[j],
         border = box_border, col = box_fill)
  }

  # Optional extremal versions consistent with the data
  if (show_extrema) {
    # Lower envelope: all jumps at q_j
    for (j in 1:m) {
      segments(x0 = if (j==1) min(xlim) else q[j-1], y0 = cs_prev[j],
               x1 = q[j], y1 = cs_prev[j], lwd = 1, col = lower_col, lty = 2)
      segments(x0 = q[j], y0 = cs_prev[j], x1 = q[j], y1 = cs[j], lwd = 1, col = lower_col, lty = 2)
    }
    segments(x0 = q[m], y0 = cs[m], x1 = max(xlim), y1 = cs[m], lwd = 1, col = lower_col, lty = 2)

    # Upper envelope: all jumps at p_j
    for (j in 1:m) {
      segments(x0 = if (j==1) min(xlim) else p[j-1], y0 = cs_prev[j],
               x1 = p[j], y1 = cs_prev[j], lwd = 1, col = upper_col, lty = 3)
      segments(x0 = p[j], y0 = cs_prev[j], x1 = p[j], y1 = cs[j], lwd = 1, col = upper_col, lty = 3)
    }
    segments(x0 = p[m], y0 = cs[m], x1 = max(xlim), y1 = cs[m], lwd = 1, col = upper_col, lty = 3)
  }

  # Observation intervals along a line below
  if (!is.null(obs_intervals)) {
    O <- as.matrix(obs_intervals)
    if (ncol(O) != 2) stop("obs_intervals must have two columns: L, R.")
    L <- O[,1]; R <- O[,2]
    # Baseline for readability
    segments(x0 = min(xlim), y0 = obs_y, x1 = max(xlim), y1 = obs_y, col = "gray70", lwd = 1, lty = 1)
    # Draw each observation interval
    for (i in seq_along(L)) {
      if (is.finite(L[i]) && is.finite(R[i])) {
        if (R[i] > L[i]) {
          segments(L[i], obs_y, R[i], obs_y, col = obs_col, lwd = obs_lwd)
          points(c(L[i], R[i]), c(obs_y, obs_y), pch = 16, cex = 0.6, col = obs_col)
        } else { # exact observation (L==R)
          if (obs_point_as_tick) {
            segments(L[i], obs_y - 0.01, L[i], obs_y + 0.01, col = obs_col, lwd = obs_lwd)
          } else {
            points(L[i], obs_y, pch = 16, col = obs_col)
          }
        }
      }
    }
    # Label
    mtext("observations (interval-censored)", side = 4, line = 0.5, at = obs_y, las = 1, cex = 0.8)
  }

  # Legend
  leg_items <- c("uncertainty box")
  leg_lty   <- c(NA)
  leg_lwd   <- c(NA)
  leg_col   <- c(box_fill)
  leg_pch   <- c(15)
  if (show_extrema) {
    leg_items <- c(leg_items, "lower (jumps at q_j)", "upper (jumps at p_j)")
    leg_lty   <- c(leg_lty, 2, 3)
    leg_lwd   <- c(leg_lwd, 1, 1)
    leg_col   <- c(leg_col, lower_col, upper_col)
    leg_pch   <- c(leg_pch, NA, NA)
  }
  if (!is.null(obs_intervals)) {
    leg_items <- c(leg_items, "observation intervals")
    leg_lty   <- c(leg_lty, 1)
    leg_lwd   <- c(leg_lwd, obs_lwd)
    leg_col   <- c(leg_col, obs_col)
    leg_pch   <- c(leg_pch, NA)
  }
  legend("bottomright", legend = leg_items, lty = leg_lty, lwd = leg_lwd,
         pch = leg_pch, pt.cex = 2, col = leg_col, bty = "n")
}
