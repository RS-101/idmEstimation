
load("testing/plots/testing_starting_vals_npmle.rdata")
plot(p)


sim_idm_list <- idmEstimation::simulate_idm_weibull(
  300,
  shape12 = 1.5, shape13 = 1.3, shape23 = 1.7,
  scale12 = 1000, scale13 = 1500, scale23 = 1240)

max_time <- max(sim_idm_list$datasets$obs$T_obs)



my_data <- sim_idm_list$datasets$obs
summary(my_data)

set.seed(17)
fits <- list(
  fit_np_1 = fit_npmle(data = my_data, tol = 1e-4, verbose = T, max_iter = 400),
  fit_np_2 = fit_npmle(data = my_data, tol = 1e-4, verbose = T, max_iter = 400),
  fit_np_3 = fit_npmle(data = my_data, tol = 1e-4, verbose = T, max_iter = 400),
  fit_np_4 = fit_npmle(data = my_data, tol = 1e-4, verbose = T, max_iter = 400),
  fit_np_5 = fit_npmle(data = my_data, tol = 1e-4, verbose = T, max_iter = 400)
)


p <- plot.idm_hazards(sim_idm_list$hazards,
                           estimator_name = "true",
                           max_time = max_time, cumulative = T)
j = 1
for(fit in fits) {
  j = j + 1
  p <- plot.idm_hazards(fit$hazards,
                   estimator_name = paste0("NPMLE_", j),
                   max_time = max_time,
                   add = p,
                   cumulative = T)
}

p


save(p, file = "testing/plots/testing_starting_vals_npmle.rdata")

