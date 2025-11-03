set.seed(17)
sim_idm_list <- idmEstimation::simulate_idm_weibull(
  300,
  shape12 = 1.5, shape13 = 1.3, shape23 = 1.7,
  scale12 = 1000, scale13 = 1500, scale23 = 1240)

max_time <- max(sim_idm_list$datasets$obs$T_obs)



my_data <- sim_idm_list$datasets$obs
summary(my_data)


fit_pl <- fit_idm(data = my_data, run_in_parallel = T)

fit_np <- fit_npmle(data = my_data, tol = 1e-4, verbose = T, max_iter = 500)


p_true <- plot.idm_hazards(sim_idm_list$hazards,
                           estimator_name = "true",
                           max_time = max_time, cumulative = T)

p_true_pl <- plot.idm_hazards(fit_pl$hazards,
                 estimator_name = "pl",
                 max_time = max_time,
                 add =  p_true,
                 cumulative = T)

plot.idm_hazards(fit_np$hazards,
                 estimator_name = "NPMLE",
                 max_time = max_time,
                 add = p_true_pl,
                 cumulative = T)
