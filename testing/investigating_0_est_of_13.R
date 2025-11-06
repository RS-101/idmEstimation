sim_list <- idmEstimation::simulate_idm_joly(300)
obs_data <- sim_list$datasets$obs
max_T_obs <- max(obs_data$T_obs)

sum_obs_data <- summarise_obs_data(observed_data = obs_data)

npmle_fit <- fit_npmle(obs_data, max_iter = 300, verbose = TRUE)
#### Comparing to naive estimator ####
naive_obs_data <- obs_data
naive_obs_data$V_healthy = ifelse(naive_obs_data$status_ill == 0, naive_obs_data$T_obs, naive_obs_data$V_healthy)

sum_naive_obs_data <- summarise_obs_data(observed_data = naive_obs_data)
sum_naive_obs_data$tables$table_death

npmle_fit_naive <- fit_npmle(naive_obs_data, max_iter = 300, verbose = TRUE)

npmle_fit_naive$hazards$A13(100)
npmle_fit_naive$hazards$A13(1000)

npmle_fit_naive$distributions_functions$F(33)
npmle_fit_naive$distributions_functions$F12(33)
npmle_fit_naive$distributions_functions$F13(33)

