# idmEstimation

**idmEstimation** provides comprehensive tools for analyzing illness-death models with interval-censored illness times and right-censored death times. The package implements multiple estimation approaches including non-parametric maximum likelihood estimation (NPMLE), piecewise-constant hazards, and penalized spline-based methods.

## Illness-Death Model

The illness-death model is a three-state progressive model commonly used in medical research and survival analysis:

- **State 1 (Healthy)** → **State 2 (Illness)**: transition hazard α₁₂(t)
- **State 1 (Healthy)** → **State 3 (Death)**: transition hazard α₁₃(t)  
- **State 2 (Illness)** → **State 3 (Death)**: transition hazard α₂₃(t)

State 3 (Death) is absorbing. Subjects can die directly from State 1 without experiencing illness.

### Data Structure

The package handles complex censoring patterns:

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
  a12 = 0.001,  # Healthy → Illness
  a13 = 0.0005, # Healthy → Death
  a23 = 0.002,  # Illness → Death
  average_number_of_visits = 10
)

# Examine data structure
head(sim_data$data)

# Summarize observation patterns
summary_stats <- summarise_obs_data(sim_data$data)
```

### Fit Models

``` r
# Non-parametric maximum likelihood estimation
fit_npmle <- fit_npmle(sim_data$data, max_iter = 100, verbose = FALSE)

# Piecewise-constant hazard model
fit_pc <- fit_pc_model(sim_data$data, n_knots = 6)

# Penalized spline model (automatic smoothing selection)
fit_spline <- fit_spline_model(
  sim_data$data, 
  degree = 3,
  n_knots = 7,
  verbose = FALSE
)
```

### Visualize Results

``` r
# Compare all fitted models with true hazards
plot(fit_npmle, fit_pc, fit_spline, sim_data)

# Plot individual model
plot(fit_spline)
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
# Simulate with Weibull hazards
sim_weibull <- simulate_idm_weibull(
  n = 500,
  shape12 = 2, scale12 = 10,
  shape13 = 3, scale13 = 20,
  shape23 = 1.5, scale23 = 8
)

# Fit and compare
fit_weibull_pc <- fit_pc_model(sim_weibull$data, n_knots = 8)
fit_weibull_spline <- fit_spline_model(
  sim_weibull$data, 
  n_knots = 10,
  verbose = FALSE
)
```

## Key Features

- **Flexible simulation**: Generate data with arbitrary time-dependent hazards
- **Multiple censoring mechanisms**: Including Frydman and Joly scenarios from the literature
- **Efficient C++ backend**: Fast likelihood computation via Rcpp and RcppArmadillo
- **Automatic smoothing selection**: Cross-validation for penalized splines
- **Comprehensive visualization**: Compare multiple model fits easily
- **Handles complex patterns**: Missing transitions, interval censoring, right censoring

## Getting Help

- Package documentation: `help(package = "idmEstimation")`
- Function help: `?fit_spline_model`
- Report issues: [GitHub Issues](https://github.com/RS-101/idmEstimation/issues)

## References
Frydman and Szarek 2009 Nonparametric estimation in a Markov "illness-death" process from interval censored observations with missing intermediate transition status
Joly et al. 2002 - A penalized likelihood approach for an illness–death model with interval‐censored data: application to age‐specific incidence of dementia
## License

MIT + file LICENSE
