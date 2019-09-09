set.seed(1037584)

library(tidyverse)
library(magrittr)
library(reshape2)
library(rstan)

## Data ************************ ####
N_series <- 10
x <- seq(from = 0, to = 15, by = 0.1)
x_int <- x[x %% 1 == 0]

rho <- rgamma(N_series, 5, 2)     # length-scale
alpha <- rgamma(N_series, 2, 2)   # standard deviation

sigma <- rep(.5, N_series)   # observation error


true_pars <- data.frame(sigma = sigma,
                        rho = rho,
                        alpha = alpha,
                        series = 1:N_series)


# Data simulator
sim_data_model <- stan_model("EWS/stan_forum_gp_generator.stan")


data_set <- lapply(1:N_series, function(i) {
  
  dat_list <- list(T = length(x),
                   x = x,
                   alpha = alpha[i],
                   rho = rho[i],
                   error = sigma[i])
  
  
  draw <- sampling(sim_data_model,
                   iter=1,
                   algorithm='Fixed_param',
                   chains = 1,
                   data = dat_list)
  
  
  samps <- rstan::extract(draw)
  df <-  with(samps, data.frame(y = y[1,], f = f[1,])) %>% 
    mutate(index = i, x = x)
  
  return(df)
  
}) %>%
  do.call(rbind, .)



## Plot data
series_plot <- ggplot(data_set) +
  geom_line(aes(x = x, y = f), color = "blue") +
  geom_line(data = data_set %>% filter((x %in% x_int)), aes(x = x, y = y), color = "red") +
  facet_wrap(~index)

series_plot



## Take integer observation times
data_set_thin <- data_set %>% 
  filter(x %in% x_int)

## Stan models ####
single_model <- stan_model("stan_models/se_single_covariance.stan")
hierarchical_model <- stan_model("stan_models/se_hierarchical_covariance.stan")


## Fit separately
separate_results <- lapply(1:N_series, function(i) {
  
  series_df <- data_set_thin %>% filter(index == i)
  
  x <- sampling(single_model,
                list(N = length(series_df$y),
                     x = series_df$x,
                     y = series_df$y),
                     iter = 1000,
                     chains = 1,
                     control = list(adapt_delta = 0.95)
                )

  
  res_df <- lapply(c("rho", "alpha", "sigma"), function(par) {
    
    res <- get_stan_results(x, paste0("^", par), regex = T) %>% select(lower2.5, mode, upper97.5)
    
    return(rbind(res) %>% mutate(parameter = par))
    
  }) %>% do.call(rbind, .) %>% mutate(series = i)
  
  return(res_df)
}
) %>%
  do.call(rbind, .) %>%
  mutate(model = "separate_models")



## Fit hierarchical model
hierarchical_samples <- sampling(hierarchical_model,
                                     list(N = length(x_int),
                                          N_series = N_series,
                                          x = x_int,
                                          y = matrix(data_set_thin[, "y"], ncol = N_series)),
                                     iter = 1000, chains = 1, control = list(adapt_delta = 0.95))
# Hierarchical results
hierarchical_results <- lapply(c("sigma", "alpha", "rho"), function(par) {
  
  cbind(rep("hierarchical", N_series),
        get_stan_results(hierarchical_samples, paste0("^", par, "\\["), regex = T) %>% select(lower2.5, mode, upper97.5)) %>%
    set_colnames(c("model", "lower2.5", "mode", "upper97.5")) %>%
    mutate(series = 1:N_series, parameter = par)
  
  
}
) %>% do.call(rbind, .) %>% mutate(model = "hierarchical")



# Plot results
estimate_plot <- ggplot() + 
  geom_errorbar(data = separate_results, aes(x = series, ymin = lower2.5, ymax = upper97.5, color = model), position ="dodge") +
  geom_line(data = true_pars %>% melt(id.vars = "series") %>% mutate(parameter = variable), aes(x = series, y = value)) +
  geom_errorbar(data = hierarchical_results, aes(x = series, ymin = lower2.5, ymax = upper97.5, color = model), position ="dodge") +
  facet_grid(~parameter, scales = "free", labeller = labeller(.rows =label_both)) +
  labs()

estimate_plot


lapply(c("alpha", "rho", "sigma"), function(par) {
  
  df <- rbind(separate_results %>% filter(parameter == par), 
  hierarchical_results %>% filter(parameter == par)) %>%
    mutate(true = rep(true_pars[, par], 2)) %>% 
    arrange(true) %>% 
    mutate(dummy = rep(1:N_series, each = 2))
    
  
  
  ggplot(df) +
    geom_line(aes(x = dummy, y = true)) + 
    geom_errorbar(aes(x = dummy, ymin = lower2.5, ymax =upper97.5, color =model), position = "dodge")
    
  
})




ggplot() +
  geom_point(data = true_pars, aes(x = rho, y = alpha)) +
  geom_point(data = separate_results %>% select(-c(lower2.5, upper97.5)) %>% spread(parameter, mode), 
             aes(x = rho, y = alpha), color = "blue") +
  geom_point(data = hierarchical_results %>% select(-c(lower2.5, upper97.5)) %>% spread(parameter, mode), 
             aes(x = rho, y = alpha), color = "red") 
