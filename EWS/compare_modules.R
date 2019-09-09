
## Data  ####
N <- 5
times <- 1:100
# length_scale <- rep(4, N)
length_scale <- rgamma(N, 10, 4)
lambda <- 1/(length_scale)
kappa <- rgamma(N, 5, 10)
# kappa <- rep(.25, N)
stat_var <- kappa
error <- rnorm_trunc(N, .25, 0.1)
# error <- rep(.25, N)
sigma <- sqrt(2*kappa*lambda)
mu <- rep(0, N)

set.seed(1)
sim_data_model <- stan_model("stan_models/generate_gp.stan")
gp_set <- generate_gp_set(N_series = N, times, covariance = "se",
                          length_scale = length_scale, stat_var = stat_var, error = error, modules = FALSE)



ggplot(gp_set) +
  geom_line(aes(x = x, y = f), color = "blue") +
  geom_line(data = gp_set %>% filter((x%%1 == 0)), aes(x = x, y = y), color = "red") +
  facet_wrap(~index)



gp_set_thin <- gp_set %>% thin

# gp_set_thin <- gp_set

# Thin matrix
gp_set_thin_matrix <- lapply(1:N, function(i) {
  
  gp_set_thin %>%
    filter(index == i) %>% 
    pull(y)
  
}) %>% do.call(cbind, .) %>% as.data.frame()




# empirical_length_scale <- lapply(i:N, function(i) {
#   
#   x <- gp_set_thin_matrix[, i]
#   ac <- acf(x)
#   dev.off()
#   ((ac$acf[2])/var(x) %>% log %>% -.)^(-1)
# }) %>% unlist 
# empirical_length_scale <- empirical_length_scale[empirical_length_scale>=0 & empirical_length_scale < 10] %>%
#   mean_and_var()
# empirical_stat_var <-   lapply(i:N, function(i) {
#   
#   x <- gp_set_thin_matrix[, i]
#   
#   x %>% var
#   
# }) %>% unlist



## Stan models ####
# gp_model <- stan_model("stan_models/gp_fit.stan")
# gp_model_with_priors <- stan_model("stan_models/gp_fit_with_priors.stan")
# oup_model <- stan_model("stan_models/oup_single_transition.stan")
# oup_model <- stan_model("stan_models/student_oup_single_transition.stan")
# oup_model <- stan_model("stan_models/student_oup_single_covariance.stan")
oup_model <- stan_model("stan_models/student_oup_single_covariance.stan")
# oup_model <- stan_model("stan_models/student_se_single_covariance.stan")
# oup_model_stat_var_prior <- stan_model("stan_models/oup_single_transition_stat_var_prior.stan")
# oup_hierarchical <- stan_model("stan_models/oup_hierarchical_transition.stan")
# oup_hierarchical <- stan_model("stan_models/student_oup_hierarchical_transition.stan")
# oup_hierarchical <- stan_model("stan_models/student_oup_hierarchical_covariance.stan")
oup_hierarchical <- stan_model("stan_models/student_se_hierarchical_covariance.stan")
# oup_hierarchical_stat_var_prior <- stan_model("stan_models/oup_hierarchical_transition_stat_var_prior.stan")



##  with and without prior information: single ******
oup_sep_res <- lapply(1:N, function(i) {
  
  
  x <- sampling(oup_model,
                list(T = length(times),
                     time = times,
                     Y = gp_set_thin %>% filter(index == i) %>% pull(y),
                     kernel = 1, 
                     kappa_prior = c(kappa[i], 0.1)),
                iter = 1000, chains = 1, control = list(adapt_delta = 0.95))
  
  # x_with_prior <- sampling(oup_model_stat_var_prior,
  #                          list(T = length(times),
  #                               time = times,
  #                               Y = gp_set_thin %>% filter(index == i) %>% pull(y),
  #                               kernel = 0, 
  #                               stat_var_approx = empirical_stat_var[i], 
  #                               length_scale_approx = empirical_length_scale),
  #                          iter = 1000, chains = 1, control = list(adapt_delta = 0.95))
  
  
  res_df <- lapply(c("length_scale", "stat_var", "epsilon", "sigma"), function(par) {
    res <- get_stan_results(x, paste0("^", par), regex = T) %>% select(lower2.5, mode, upper97.5)
    # res_with_prior <- get_stan_results(x_with_prior, paste0("^", par), regex = T) %>%
    #   select(lower2.5, mode, upper97.5)
    
    
    
    return(rbind(res) %>% mutate(parameter = par))
    
  }) %>% do.call(rbind, .)
  
  return(res_df)
}
) %>% do.call(rbind, .) %>% mutate(model = "separate_gp")

# oup_sep_res <- oup_sep_res %>% mutate(index = rep(1:N, each = 8), prior = rep(c("no", "yes"), 4*N), model = "separate")

oup_sep_res <- oup_sep_res %>% mutate(index = rep(1:N, each = 4), model = "separate")

## Hierarchical ********
oup_hierarchical_samples <- sampling(oup_hierarchical,
                                     list(N_times = length(times),
                                          N_series = N,
                                          time = times,
                                          Y = gp_set_thin_matrix %>% t, 
                                          kappa_prior = data.frame(kappa, rep(0.1, N))),
                                     iter = 1000, chains = 1, control = list(adapt_delta = 0.95))
# Hierarchical results
oup_hier_res <- lapply(c("sigma", "stat_var", "length_scale"), function(par) {
  
  cbind(rep("hierarchical", N),
        get_stan_results(oup_hierarchical_samples, paste0("^", par, "\\["), regex = T) %>% select(lower2.5, mode, upper97.5)) %>%
    set_colnames(c("model", "lower2.5", "mode", "upper97.5")) %>%
    mutate(index = 1:N, parameter = par)
  
  
}
) %>% do.call(rbind, .) %>% mutate(prior = "no")




# oup_hierarchical_with_prior_samples <- sampling(oup_hierarchical_stat_var_prior,
#                                                 list(N_times = length(times),
#                                                      N_series = N,
#                                                      time = times,
#                                                      Y = gp_set_thin_matrix %>% t, 
#                                                      empirical_stat_var = empirical_stat_var %>% mean_and_var,
#                                                      empirical_length_scale = empirical_length_scale),
#                                                 iter = 1000, chains = 1, control = list(adapt_delta = 0.95))
# # Hierarchical results
# oup_hier_with_prior_res <- lapply(c("sigma", "stat_var", "length_scale"), function(par) {
#   
#   cbind(rep("hierarchical", N),
#         get_stan_results(oup_hierarchical_with_prior_samples, paste0("^", par, "\\["), regex = T) %>% select(lower2.5, mode, upper97.5)) %>%
#     set_colnames(c("model", "lower2.5", "mode", "upper97.5")) %>%
#     mutate(index = 1:N, parameter = par)
#   
#   
# }
# ) %>% do.call(rbind, .) %>% mutate(prior = "yes")




hier_res <- oup_hier_res %>% mutate(model = "hierarchical")


true_pars <- data.frame(epsilon = rep(error, N),
                        length_scale = length_scale,
                        stat_var = stat_var,
                        index = 1:N, 
                        sigma = sigma)

ggplot() + 
  geom_errorbar(data = oup_sep_res, aes(x = index, ymin = lower2.5, ymax = upper97.5, color = model), position ="dodge") +
  geom_line(data = true_pars %>% melt(id.vars = "index") %>% mutate(parameter = variable), aes(x = index, y = value)) +
  geom_errorbar(data = hier_res, aes(x = index, ymin = lower2.5, ymax = upper97.5, color = model), position ="dodge") +
  facet_grid(~parameter, scales = "free", labeller = labeller(.rows =label_both)) +
  labs()








