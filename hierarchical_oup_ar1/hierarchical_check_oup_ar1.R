# Check hierarchical AR1 and OUP models' performance on data generated from same models


## OUP ********************************* ####

## Data ******************

generate_oup_set <- function(N_series, times, mu, lambda, sigma, epsilon = 0, seed = 1) {
  
  set.seed(seed)
  
  res_df <- lapply(1:N_series, function(i) {
    
    kappa <- (sigma[i]^2)/(2*lambda[i])
    
    x <- rep(NA, length(times))
    x[1] <- rnorm(1, mu[i], sqrt(kappa))
    
    for(t in 2:length(times)) {
      dt <- times[t] - times[t-1]
      x[t] <- rnorm(1, 
                    mu[i] - (mu[i] - x[t-1])*exp(-lambda[i]*dt),
                    sqrt(kappa))
      
    }
    
    # Add measurement error
    if(epsilon !=0) {
      x <- x + rnorm(length(times), 0, epsilon)
    }
    
    return(x)
    
  }) %>% do.call(cbind, .) %>% 
    as.data.frame()
  
  # %>% 
  #   cbind(., times)
  # 
  
  return(res_df)
}
generate_gp_set <- function(N_series, times, covariance, length_scale, stat_var, error = 0, seed = 1) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Model
  # sim_data_model <- stan_model("stan_models/generate_gp.stan")
  
  # Set covariance
  if(covariance == "oup") {
    kernel <- 0
  } else if(covariance == "se") {
    kernel <- 1
  }
  
  
  time_grid <- seq(from = times[1], to = times[length(times)], by = 0.1)
  
  # Generate
  res <- lapply(1:N_series, function(i) {
    
    dat_list <- list(T = length(time_grid),
                     time = time_grid,
                     stat_var = stat_var[i],
                     length_scale = length_scale[i],
                     error = error,
                     kernel = kernel)
    
    
    draw <- sampling(sim_data_model,
                     iter=1,
                     algorithm='Fixed_param',
                     chains = 1,
                     data = dat_list)
    
    
    samps <- rstan::extract(draw)
    plt_df = with(samps, data.frame(y = y[1,], f = f[1,])) %>% 
      mutate(index = i)
    
    return(plt_df)
    
  }) %>%
    do.call(rbind, .) %>% 
    cbind(x = time_grid)
  
  return(res)
}

N <- 10
times <- 1:50

length_scale <- rgamma(N, 10, 3)
lambda <- 1/(length_scale)
kappa <- rgamma(N, 5, 20)
stat_var <- kappa
error <- .25
sigma <- sqrt(2*kappa*lambda)
mu <- rep(0, N)

set.seed(1)
gp_set <- generate_gp_set(N_series = N, times, covariance = "se",
                           length_scale = length_scale, stat_var = stat_var, error)



ggplot(gp_set) + 
  geom_line(aes(x = x, y = f), color = "blue") + 
  geom_line(data = gp_set %>% filter((x%%1 == 0)), aes(x = x, y = y), color = "red") + 
  facet_wrap(~index)



gp_set_thin <- gp_set %>% thin



# Stan*******************
oup_hierarchical <- stan_model("stan_models/oup_hierarchical_transition.stan")
# oup_hierarchical_covariance <- stan_model("stan_models/oup_hierarchical_covariance.stan")
oup_model <- stan_model("stan_models/oup_single_transition.stan")
# oup_model_covariance <- stan_model("stan_models/oup_single_covariance.stan")
# oup_model_covariance <- stan_model("stan_models/oup_single_covariance2.stan")
# oup_model_scaled <- stan_model("stan_models/oup_single_transition_transformed.stan")
gp_model <- stan_model("stan_models/gp_fit.stan")


# oup_hierarchical_samples <- sampling(oup_hierarchical,
#                                  list(N_times = length(times),
#                                       N_series = N,
#                                       time = times,
#                                       Y = oup_set %>% t),
#                                  iter = 1000, chains = 1)
# Hierarchical results
# oup_hier_res <- lapply(c("lambda", "sigma",  "kappa", "inv_lambda"), function(par) {
# 
#   cbind(rep("hierarchical", N),
#         get_stan_results(oup_hierarchical_samples, paste0("^", par, "\\["), regex = T) %>% select(lower2.5, mode, upper97.5)) %>%
#     set_colnames(c("model", "lower2.5", "mode", "upper97.5")) %>%
#     mutate(index = 1:N, parameter = par)
# 
# 
# }
# ) %>% do.call(rbind, .)





# oup_hierarchical_cov_samples <- sampling(oup_hierarchical_covariance,
#                                  list(N_times = length(times),
#                                       N_series = N,
#                                       time = times,
#                                       Y = oup_set %>% t),
#                                  iter = 1000, chains = 1)
# # Hierarchical results
# oup_hier_cov_res <- lapply(c("lambda", "sigma", "mu", "kappa"), function(par) {
#   
#   cbind(rep("hierarchical_covariance", N),
#         get_stan_results(oup_hierarchical_cov_samples, paste0("^", par, "\\["), regex = T) %>% select(lower2.5, mode, upper97.5)) %>%
#     set_colnames(c("model", "lower2.5", "mode", "upper97.5")) %>%
#     mutate(index = 1:N, parameter = par)
#   
#   
# }
# ) %>% do.call(rbind, .)



# # Separate samples and results

oup_sep_res <- lapply(1:N, function(i) {
  
  
  x <- sampling(gp_model,
                list(T = length(times),
                     time = times,
                     Y = gp_set_thin %>% filter(index == 1) %>% pull(y),
                     kernel = 1),
                iter = 1000, chains = 1, control = list(adapt_delta = 0.95))
  
  
  res_df <- lapply(c("length_scale", "stat_var", "error"), function(par) {
    res <- get_stan_results(x, paste0("^", par), regex = T) %>% select(lower2.5, mode, upper97.5)
    
    return(res %>% mutate(parameter = par))
    
  }) %>% do.call(rbind, .)
  
  return(res_df)
}
) %>% do.call(rbind, .) %>% mutate(model = "separate_gp", index = rep(1:N, each = 3))



# oup_sep_res <- lapply(1:N, function(i) {
# 
#   x <- sampling(oup_model,
#                 list(T = length(times),
#                      time = times,
#                      Y = oup_set[, i]),
#                 iter = 2000, chains = 1)
# 
# 
#   res_df <- lapply(c("lambda", "sigma", "kappa", "inv_lambda"), function(par) {
#     res <- get_stan_results(x, paste0("^", par), regex = T) %>% select(lower2.5, mode, upper97.5)
# 
#     return(res %>% mutate(parameter = par))
# 
#   }) %>% do.call(rbind, .)
# 
#   return(res_df)
# }
# ) %>% do.call(rbind, .) %>% mutate(model = "separate_transition", index = rep(1:N, each = 4))

# With covariance parameterixation
# oup_cov_res <- lapply(1:N, function(i) {
# 
#   x <- sampling(oup_model_covariance,
#                 list(T = length(times),
#                      time = times,
#                      Y = oup_set[, i]),
#                 iter = 2000, chains = 1, 
#                 control = list(adapt_delta = 0.95))
# 
# 
#   res_df <- lapply(c("lambda", "sigma", "eta", "kappa"), function(par) {
#     res <- get_stan_results(x, paste0("^", par), regex = T) %>% select(lower2.5, mode, upper97.5)
# 
#     return(res %>% mutate(parameter = par))
# 
#   }) %>% do.call(rbind, .)
# 
#   return(res_df)
# }
# ) %>% do.call(rbind, .) %>% mutate(model = "separate_covariance")

# Separate samples and results: SCALED
# oup_sep_res_scaled <- lapply(1:N, function(i) {
# 
#   x <- sampling(oup_model_scaled,
#                 list(T = length(times),
#                      time = times,
#                      Y = oup_set[, i]),
#                 iter = 2000, chains = 1)
# 
# 
#   res_df <- lapply(c("lambda", "sigma", "mu", "kappa"), function(par) {
#     res <- get_stan_results(x, paste0("^", par), regex = T) %>% select(lower2.5, mode, upper97.5)
# 
#     return(res %>% mutate(parameter = par))
# 
#   }) %>% do.call(rbind, .)
# 
#   return(res_df)
# }
# ) %>% do.call(rbind, .) %>% mutate(model = "separate_scaled", index = rep(1:N, each = 4))



  

# true_pars <- data.frame(value = c(lambda, sigma, mu, kappa),
#                         parameter = rep(c("lambda", "sigma", "mu", "kappa"), each = N),
#                         index = rep(1:N, 4),
#                         dummy = rep(rank(sigma), 4)


true_pars <- data.frame(value = c(length_scale, stat_var, rep(error, N)),
                        parameter = rep(c("length_scale", "stat_var", "error"), each = N),
                        index = rep(1:N, 3),
                        dummy = rep(rank(sigma), 3))


#   
# 
# 
# oup_hier_res <- oup_hier_res %>% 
#   mutate(dummy = rep(rank(sigma), 4))
# 
# oup_sep_res <- oup_sep_res %>% 
#   arrange(parameter) %>% 
#   mutate(dummy = rep(rank(sigma), 4))

df <- rbind(oup_sep_res)


ggplot() +
  geom_line(data = true_pars, aes(x = index,  y= value)) +
  geom_errorbar(data = df, aes(x = index, ymin = lower2.5, ymax = upper97.5, color = model), position = "dodge") +
  facet_wrap(~parameter, scales = "free")



ggplot() +
  geom_line(data = true_pars, aes(x = dummy,  y= value)) +
  geom_errorbar(data = df, aes(x = dummy, ymin = lower2.5, ymax = upper97.5, color = model), position = "dodge") +
  facet_wrap(~parameter, scales = "free")



## Results ***************

## AR1 ********************************* ####

## Data ************
generate_ar1_set <- function(N_series, times, mu, lambda, sigma, seed = 1) {
  
  set.seed(seed)
  
  res_df <- lapply(1:N_series, function(i) {
    
    x <- rep(NA, N_series)
    x[1] <- rnorm(1, (lambda[i]*mu[i])/(1+lambda[i]), sigma/sqrt(1 - lambda^2))
    
    for(t in 2:length(times)) {
      dt <- times[t] - times[t-1]
      x[t] <- rnorm(1, 
                    lambda[i]*mu[i] - lambda[i]*x[t-1],
                    sigma[i]^2)
      
    }
    
    return(x)
    
  }) %>% do.call(cbind, .) %>% 
    as.data.frame()
  
  # %>% 
  #   cbind(., times)
  # 
  
  return(res_df)
}
lambda <- runif(10, 0, 1)
sigma <- runif(10, 0, 1)
ar1_set <- generate_ar1_set(10, 1:100, mu, lambda, sigma)


ar1_hierarchical <- stan_model("stan_models/ar1_hierarchical.stan")


ar1_hierarchical_res <- sampling(ar1_hierarchical, 
                                 list(N_times = 100, 
                                      N_series = 10,
                                      time = 1:100,
                                      Y = ar1_set %>% t), 
                                 iter = 500, chains = 1)


plot(sigma, get_stan_results(oup_hierarchical_res, parameter = "sigma\\[", regex = TRUE) %>% pull(mode))
