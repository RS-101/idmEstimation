set.seed(17)
library(idmEstimation)

n <- 500
i <- 100
sim_idm_list <- idmEstimation::simulate_idm_constant_hazards(
  n,a12 = 0.0001, a13 = 0.0003, a23 = 0.0005)

exact_data <- sim_idm_list$datasets$exact_idm



my_data <- data.frame(
  id = 1:n,
  V_0 = rep(0, n),
  V_healthy = ifelse(!is.infinite(exact_data$T_ill), exact_data$T_ill -  runif(n,0,i), exact_data$T_death),
  V_ill = ifelse(is.finite(exact_data$T_ill), exact_data$T_ill + runif(n,0,i), NA),
  T_obs = exact_data$T_death,
  status_dead = 1,
  status_ill = as.numeric(is.finite(exact_data$T_ill))
)

fit_np <- fit_npmle(data = my_data, tol = 1e-4, verbose = T, max_iter = 500)


sim_idm_list$hazards$A12(500)
fit_np$hazards$A12(500)



sim_idm_list$hazards$A23(500)
fit_np$hazards$A23(500)
