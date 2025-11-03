# Check if CPP implementation of likelihood calculations are correct.

library(idmEstimation)
source("../idmEstimation/old_R_implementation/old_penlik.r")
dat <- idmEstimation::simulate_idm_weibull(1000,
                            shape12 = 1.1, scale12 = 1/0.0008,
                            shape13 = 1.8, scale13 = 1/0.0002,
                            shape23 = 1.3, scale23 = 1/0.0016)

names(dat$datasets$obs)
summary(dat$datasets$obs$case)

observed_data <- dat$datasets$obs

calc_difference_case_in_ll <- function(case = 1, observed_data, theta_val, ll_function) {
  case_data <- observed_data[as.numeric(observed_data$case) == case, ]
  summary(case_data)
  case_max_T_obs <- max(case_data$T_obs)
  case_knots <- seq(0, case_max_T_obs, length.out = 7)

  model_config <- setup_cpp_model(case_data$V_0,
                      case_data$V_healthy,
                      case_data$V_ill,
                      case_data$T_obs,
                      case_data$status_dead,
                      case_data$status_ill,
                      n_knots = 7,
                      knots_12 = case_knots,
                      knots_13 = case_knots,
                      knots_23 = case_knots,
                      degree = 3)

  theta_init <- rep(theta_val, model_config$n_theta_12)

  new_case <- idmEstimation::calc_log_likelihood(model_config$model_pointer,
                                    theta_12 = theta_init,
                                    theta_13 = theta_init,
                                    theta_23 = theta_init)


  theta_init_list <- list(a12 = theta_init,
                    a13 = theta_init,
                    a23 = theta_init)

  knots_list <- list(a12 = case_knots,
                    a13 = case_knots,
                    a23 = case_knots)


  case_hazards <- make_all_3_spline_funs(theta = theta_init_list, knots = knots_list, degree = 3)

  old_case <- sum(log(ll_function(case_data$V_0,
                    case_data$V_healthy,
                    case_data$V_ill,
                    case_data$T_obs,
                    a12 = case_hazards$a12,
                    a13 = case_hazards$a13,
                    a23 = case_hazards$a23,
                    A12 = case_hazards$A12,
                    A13 = case_hazards$A13,
                    A23 = case_hazards$A23)))

  print(old_case)
  print(new_case)

  (new_case - old_case)/old_case
}

calc_difference_in_full_ll <- function(observed_data, theta_val) {
  case_max_T_obs <- max(observed_data$T_obs)
  case_knots <- seq(0, case_max_T_obs, length.out = 7)

  model_config <- setup_cpp_model(observed_data$V_0,
                                  observed_data$V_healthy,
                                  observed_data$V_ill,
                                  observed_data$T_obs,
                                  observed_data$status_dead,
                                  observed_data$status_ill,
                                  n_knots = 7,
                                  knots_12 = case_knots,
                                  knots_13 = case_knots,
                                  knots_23 = case_knots,
                                  degree = 3)

  theta_init <- rep(theta_val, model_config$n_theta_12)

  new_case <- idmEstimation::calc_log_likelihood(model_config$model_pointer,
                                                 theta_12 = theta_init,
                                                 theta_13 = theta_init,
                                                 theta_23 = theta_init)


  theta_init_list <- list(a12 = theta_init,
                          a13 = theta_init,
                          a23 = theta_init)

  knots_list <- list(a12 = case_knots,
                     a13 = case_knots,
                     a23 = case_knots)

  old_case <- cal_log_likehood(observed_data, theta = theta_init_list, knots = knots_list, degree = 3)$ll_value

  print(new_case)
  print(old_case)
  (new_case - old_case)/old_case
}

calc_difference_in_pen_ll <- function(observed_data, theta_val, kappa_term = c(1e8,1e8,1e8)) {
  case_max_T_obs <- max(observed_data$T_obs)
  case_knots <- seq(0, case_max_T_obs, length.out = 7)

  model_config <- setup_cpp_model(observed_data$V_0,
                                  observed_data$V_healthy,
                                  observed_data$V_ill,
                                  observed_data$T_obs,
                                  observed_data$status_dead,
                                  observed_data$status_ill,
                                  n_knots = 7,
                                  knots_12 = case_knots,
                                  knots_13 = case_knots,
                                  knots_23 = case_knots,
                                  degree = 3)

  theta_init <- rep(theta_val, model_config$n_theta_12)

  new_case <- idmEstimation::calc_full_penalty(model_config$model_pointer,
                                               theta_12 = theta_init,
                                               theta_13 = theta_init,
                                               theta_23 = theta_init,
                                               kappa_12= kappa_term[1],
                                               kappa_13 = kappa_term[2],
                                               kappa_23 = kappa_term[3])


  theta_init_list <- list(a12 = theta_init,
                          a13 = theta_init,
                          a23 = theta_init)

  knots_list <- list(a12 = case_knots,
                     a13 = case_knots,
                     a23 = case_knots)

  old_case <- cal_pen_log_likehood(observed_data, theta = theta_init_list, knots = knots_list, degree = 3, kappa_term)

  print(new_case)
  print(old_case$penalty)
  old_case$penalty - new_case
}

#### Test loop ####
for(i in 1:100){
  theta_val <- runif(1)

  val <- calc_difference_case_in_ll(case = 1, observed_data, theta_val, case_1_likelihood)
  val <- val + calc_difference_case_in_ll(case = 2, observed_data, theta_val, case_2_likelihood)
  val <- val + calc_difference_case_in_ll(case = 3, observed_data, theta_val, case_3_likelihood)
  val <- val + calc_difference_case_in_ll(case = 4, observed_data, theta_val, case_4_likelihood)
  val <- val + calc_difference_in_full_ll(observed_data, theta_val = 0.5)
  val <- val + calc_difference_in_pen_ll(observed_data, theta_val = 0.5)

  if(val > 5e-4) stop("stopped")
}


#### Diff in optim ####
calc_difference_in_optim <- function(observed_data, theta_val, kappa_term = c(1e8,1e8,1e8)) {
  case_max_T_obs <- max(observed_data$T_obs)
  case_knots <- seq(0, case_max_T_obs, length.out = 7)

  model_config <- setup_cpp_model(observed_data$V_0,
                                  observed_data$V_healthy,
                                  observed_data$V_ill,
                                  observed_data$T_obs,
                                  observed_data$status_dead,
                                  observed_data$status_ill,
                                  n_knots = 7,
                                  knots_12 = case_knots,
                                  knots_13 = case_knots,
                                  knots_23 = case_knots,
                                  degree = 3)

  theta_init <- rep(theta_val, model_config$n_theta_12)

  new_case <- max_pen_likelihood(model_config,
                                               kappa_12= kappa_term[1],
                                               kappa_13 = kappa_term[2],
                                               kappa_23 = kappa_term[3])

  new_case_un <- max_pen_likelihood_unconstrained(model_config,
                                 kappa_12= kappa_term[1],
                                 kappa_13 = kappa_term[2],
                                 kappa_23 = kappa_term[3])

  knots_list <- list(a12 = case_knots,
                     a13 = case_knots,
                     a23 = case_knots)

  old_case <- do_likelihood_optim(observed_data,kappa_term, knots = knots_list, n_knots = 7, degree = 3)


  old_case_un <- do_likelihood_optim_unconstrained(observed_data,kappa_term, knots = knots_list, n_knots = 7, degree = 3)


  bm <- bench::mark(
    new_un =max_pen_likelihood_unconstrained(model_config,
                                             kappa_12= kappa_term[1],
                                             kappa_13 = kappa_term[2],
                                             kappa_23 = kappa_term[3]),
    new =max_pen_likelihood(model_config,
                            kappa_12= kappa_term[1],
                            kappa_13 = kappa_term[2],
                            kappa_23 = kappa_term[3]),
    old =  do_likelihood_optim(observed_data,kappa_term, knots = knots_list, n_knots = 7, degree = 3),
    check = F
  )

  print(bm)

  list(old = old_case,
       new = new_case,
       old_un = old_case_un,
       new_un = new_case_un,
       bm = bm)








}

res <- calc_difference_in_optim(observed_data,theta_val = 0.5)
res$new$theta_hat$theta_12
res$old$theta_hat$a12
res$new$penalized_log_likelihood
res$old$pl_value

res$new_un$theta_hat$theta_12
res$old_un$theta_hat$a12
res$new_un$penalized_log_likelihood
res$old_un$pl_value



res$new$penalized_log_likelihood
res$old$pl_value

res$new$log_likelihood
res$old$ll_value




