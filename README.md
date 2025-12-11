# idmEstimation

**idmEstimation** provides comprehensive tools for analyzing illness-death models with interval-censored illness times and right-censored death times. The package implements multiple estimation approaches including non-parametric maximum likelihood estimation (NPMLE), piecewise-constant hazards, and penalized spline-based methods.

## Illness-Death Model

The illness-death model is three state model commonly used in medical research and survival analysis:

- **State 1 (Healthy)** → **State 2 (Illness)**
- **State 1 (Healthy)** → **State 3 (Death)**
- **State 2 (Illness)** → **State 3 (Death)**

### Data Structure

The package handles a common censoring patterns:

- **Illness times** are interval-censored (observed between two visits)
- **Death times** are right-censored
- **Missing transitions** when illness status is unknown at observation end

## Installation

You can install the development version of idmEstimation from GitHub:

```r
# install.packages("devtools")
devtools::install_github("RS-101/idmEstimation")
```

## Quick Start

### Simulate Data

``` r
library(idmEstimation)

# Simulate illness-death data with constant hazards
set.seed(123)
sim_data <- simulate_idm_constant_hazards(
  n = 300,
  a12 = 0.0008,
  a13 = 0.0002,
  a23 = 0.0016,
  mean_time_between_visits = 200,
  sd_time_between_visits = 20,
  max_right_censoring_time = 10000,
  prob_censoring_at_last_visit = 0.2
)

# Print summary
summary(sim_data)   # % missing transition etc.

# Save summary as list
summary_stats <- summary(sim_data)
summary_stats$tables$true_vs_observed

# Observed data
head(sim_data$data)
```

### Fit Models

``` r
# Non-parametric maximum likelihood estimation
fit_npmle <- fit_npmle(sim_data$data, max_iter = 500, verbose = FALSE)

# Piecewise-constant hazard model
fit_pc <- fit_pc_model(sim_data$data, n_knots = 6)

# Penalized spline model (remove kappa values to run CV)
fit_spline <- fit_spline_model(
  sim_data$data, 
  n_knots = 4,
  kappa_12 = 4e15,
  kappa_13 = 1e11,
  kappa_23 = 4e15,
  verbose = FALSE,
  run_in_parallel = FALSE # Only tested on debian system.
)
```

### Visualize Results

``` r
# Plot one model alone
p <- plot_fit(fit_pc, list(entry_time = 100))
p$hazards


# Combine with data generating 
p1 <- plot_fit(fit_pc, sim_data, list(entry_time = 100))
p1$cumulative_hazards

# Combine all
p2 <- plot_fit(fit_npmle, fit_spline, fit_pc, sim_data, list(entry_time = 100))
p2$distributions

```

### Extract Estimates

``` r
# Evaluate hazard functions at specific times
time_points <- c(5, 10, 20, 30, 40)

# True hazards (from simulation)
true_a12 <- sim_data$estimators$hazard_functions$a12(time_points)

# Estimated hazards (from spline model)
est_a12 <- fit_spline$estimators$hazard_functions$a12(time_points)

# Compare
data.frame(
  time = time_points,
  true = true_a12,
  estimated = est_a12
)

# Cumulative hazards
fit_spline$estimators$cum_hazard_functions$A12(time_points)

# Probability of illness by time t
fit_spline$estimators$distribution_functions$F12(time_points)

# Survival in state 2 conditional on entry at t=10
fit_spline$estimators$distribution_functions$P22(time_points, entry_time = 10)
```

## Advanced Usage

### Custom Hazard Functions

``` r

estimators <- create_weibull_hazard()


exact_data <- simulate_exact_idm(
  n = 300,
  a12 = estimators$hazard_functions$a12,
  a13 = estimators$hazard_functions$a13,
  a23 = estimators$hazard_functions$a23,
  )

censored_data <- add_censoring_type_1(
  exact_idm = exact_data,
  mean_time_between_visits = 10,
  sd_time_between_visits = 2,
  max_right_censoring_time = 300,
  prob_censoring_at_last_visit = 0.2
  )

fit <- fit_pc_model(censored_data$obs)
p <- plot_fit(fit_pc, list(entry_time = 100))
p$distributions
```

## Key Features

- **Flexible simulation**: Generate data with arbitrary time-dependent hazards
- **Multiple censoring mechanisms**: Including Frydman and Joly scenarios from the literature
- **Efficient C++ backend**: Fast likelihood computation via Rcpp and RcppArmadillo
- **Automatic smoothing selection**: Cross-validation for penalized splines
- **Visualization**: Compare multiple model fits easily

## Getting Help

- Package documentation: `help(package = "idmEstimation")`
- Function help: `?fit_spline_model`
- Report issues: [GitHub Issues](https://github.com/RS-101/idmEstimation/issues)

## References

Frydman and Szarek 2009 Nonparametric estimation in a Markov "illness-death" process from interval censored observations with missing intermediate transition status

Joly et al. 2002 - A penalized likelihood approach for an illness–death model with interval‐censored data: application to age‐specific incidence of dementia

## License

MIT + file LICENSE
