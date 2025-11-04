library(idmEstimation)

set.seed(100)

n <- 1000
missing_distance = 10
sim_idm_list <- idmEstimation::simulate_idm_weibull(
  n,
  shape12 = 1.5, shape13 = 1.3, shape23 = 1.7,
  scale12 = 700, scale13 = 1500, scale23 = 1000)


exact_data <- sim_idm_list$datasets$exact_idm

summary(exact_data[exact_data$path == "1->3",])
summary(exact_data[exact_data$path == "1->2->3",])



head(exact_data)
head(sim_idm_list$datasets$obs)

p <- plot.idm_hazards(sim_idm_list$hazards,
                      estimator_name = "true",
                      max_time = rep(max(exact_data$T_death),3), cumulative = T)

my_data <- data.frame(
  id = 1:n,
  V_0 = rep(0, n),
  V_healthy = ifelse(!is.infinite(exact_data$T_ill), exact_data$T_ill -  runif(n,0,missing_distance), exact_data$T_death),
  V_ill = ifelse(is.finite(exact_data$T_ill), exact_data$T_ill + runif(n,0,missing_distance), NA),
  T_obs = exact_data$T_death,
  status_dead = 1,
  status_ill = as.numeric(is.finite(exact_data$T_ill))
)

max_time <- c(quantile(my_data$V_ill,probs = 0.8, na.rm = T),
              max(my_data$T_obs[my_data$status_ill == 0]),
              max(my_data$T_obs))

summary(my_data)




fit_np = fit_npmle(data = my_data, tol = 1e-4, verbose = T, max_iter = 1000)

p <- plot.idm_hazards(fit_np$hazards,
                      estimator_name = paste("NPMLE_", missing_distance),
                      max_time = max_time,
                      add = p,
                      cumulative = T)



p + ggplot2::xlim(c(0,1000)) + ggplot2::ylim(c(0,2))

