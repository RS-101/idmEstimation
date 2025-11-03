set.seed(17)
sim_idm_list <- idmEstimation::simulate_idm_constant_hazards(
  300,a12 = 0.0001, a13 = 0.0003, a23 = 0.0005)

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
