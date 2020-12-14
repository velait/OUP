## Cusp model: one series, no observation model


## Data ********************************************** ####
N <- 9
puolet <- 5

cusp_parameters <- c("alpha", "beta", "lambda", "epsilon", "r" )

r <- runif(N, .5, 2)
# alpha <- c(0, 0, 0, 0, 0,
#            -0.5, -0.5, -1, 1)
# beta <- c(2, -2, 1.5, -1.5, -1,
#           1, -1,  0, 0)
alpha <- rnorm(N, 0, .5)
beta <- rnorm(N, 0, 2)
lambda <- rnorm(N, 0, .25)
log_epsilon <- rnorm(N, 0, 0.5)
epsilon <- exp(log_epsilon)

n_points <- 100

grid <- seq(from = 0, to = n_points, by = 0.1)
seed = sample(1:n_points, 1)


inits <- rnorm(N, lambda[1], 0.1)

df_shoji <- lapply(1:N, function(i) shoji_generator(y0 = inits[i], times = grid,
                                                    r = r, alpha = alpha[i], beta = beta[i],
                                                    lambda = lambda[i], epsilon = epsilon[i],
                                                    seed = sample(1:1000, 1))) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  cbind(time = grid) %>%
  set_colnames(c(paste0("y", 1:N), "x"))


df_shoji <- thin(df_shoji)

p <- df_shoji %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value)) +
  geom_line() +
  facet_wrap(~variable) +
  labs(title = "Cusp", subtitle = "100 observations, dt = 1")


q <- df_shoji %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x =  value)) +
  geom_density() +
  facet_wrap(~variable) +
  labs(title = "Cusp", subtitle = "100 observations, dt = 1")



## Inference ***************************************** ####

shoji_cusp_with_r <- stan_model("cusp/shoji_cusp_with_r.stan")
shoji_cusp <- stan_model("cusp/shoji_cusp.stan")


shoji_cusp_with_r_samples <- lapply(1:N, function(i) {
  
  stan_data <- list(N_obs = nrow(df_shoji),
                    N_OTUs = 1,
                    x = df_shoji$x,
                    y = df_shoji[, i] %>% as.matrix() %>% t)
  
  
  # samples <- sampling(shoji_cusp_with_r,
  #                                           stan_data,
  #                                           iter = 1000,
  #                                           chains = 1,
  #                                           control = list(adapt_delta = 0.8),
  #                     init = list(list(alpha = alpha[i] %>% as.array(),
  #                                      beta = beta[i] %>% as.array(),
  #                                      lambda = lambda[i] %>% as.array(),
  #                                      epsilon = epsilon[i] %>% as.array(),
  #                                      r = r[i] %>% as.array)))
  
  samples <- sampling(shoji_cusp_with_r,
                      stan_data,
                      iter = 1000,
                      chains = 1,
                      control = list(adapt_delta = 0.8))

  
  samples
})

shoji_cusp_samples <- lapply(1:N, function(i) {
  
  stan_data <- list(N_obs = nrow(df_shoji),
                    N_OTUs = 1,
                    x = df_shoji$x,
                    y = df_shoji[, i] %>% as.matrix() %>% t)
  
  
  # samples <- sampling(shoji_cusp,
  #                                           stan_data,
  #                                           iter = 1000, 
  #                                           chains = 1, 
  #                                           control = list(adapt_delta = 0.8),
  #                     init = list(list(alpha = alpha[i] %>% as.array(),
  #                                 beta = beta[i] %>% as.array(),
  #                                 lambda = lambda[i] %>% as.array(),
  #                                 epsilon = epsilon[i] %>% as.array())))
  
  
  samples <- sampling(shoji_cusp,
                      stan_data,
                      iter = 1000, 
                      chains = 1, 
                      control = list(adapt_delta = 0.8))
  
  samples
})


## Results ******************************************* ####
results_df_with_r <- lapply(1:N, function(i) {
  
  lapply(1:length(cusp_parameters),
         function(x) get_stan_results(shoji_cusp_with_r_samples[[i]],
                                      paste0("theta\\[.+,", x,"\\]"))) %>%
    do.call(rbind, .) %>% 
    cbind(parameter = cusp_parameters, 
          series = i, 
          true = c(alpha[i], beta[i], lambda[i], epsilon[i], r[i]))
  
}) %>% do.call(rbind, .)

results_df <- lapply(1:N, function(i) {
  
  lapply(1:length(cusp_parameters[1:4]),
         function(x) get_stan_results(shoji_cusp_samples[[i]],
                                      paste0("theta\\[.+,", x,"\\]"))) %>%
    do.call(rbind, .) %>% 
    cbind(parameter = cusp_parameters[1:4], 
          series = i, 
          true = c(alpha[i], beta[i], lambda[i], epsilon[i]))
  
}) %>% do.call(rbind, .)
  
  
## Plot estimates ************************************ ####


results_df_with_r %>% 
  ggplot() +
  geom_errorbar(aes(x = as.factor(series), ymin = lower2.5, ymax = upper97.5)) + 
  geom_point(aes(x = as.factor(series), y = true)) + 
  facet_wrap(~parameter)


results_df %>% 
  ggplot() +
  geom_errorbar(aes(x = as.factor(series), ymin = lower2.5, ymax = upper97.5)) + 
  geom_point(aes(x = as.factor(series), y = true)) + 
  facet_wrap(~parameter)


## Draw samples, simulate, stationary density ******** ####

n_sim <- 25


simulated_densities_with_r <- lapply(1:N, function(i) {
  
  # number of samples
  n_samples <- shoji_cusp_with_r_samples[[i]]@stan_args[[1]]$iter - shoji_cusp_with_r_samples[[i]]@stan_args[[1]]$warmup
  
  # samples indices
  sample_ind <- sample(1:n_samples, n_sim)
  
  par_samples <- rstan::extract(shoji_cusp_with_r_samples[[i]], "theta")[[1]][sample_ind, , ]
  
  
  
  # simulate series
  ticks <- 0.1
  series_sim <- lapply(1:n_sim, function(s) {
    
      shoji_generator(y0 = rnorm(1, 5, 1), times = seq(from=0, to=n_points, by = ticks),
                      # r = 1,
                      r = par_samples[s, 5],
                      alpha = par_samples[s, 1],
                      beta = par_samples[s, 2],
                      lambda = par_samples[s, 3],
                      epsilon = par_samples[s, 4],
                      seed = sample(1:1000, 1)) %>% 
      as.data.frame() %>%
      cbind(time = seq(from=0, to=n_points, by = ticks)) %>%
      mutate(sample = s)
  }) %>% 
    do.call(rbind, .) %>% 
    cbind(series = i) %>% 
    set_colnames(c("value", "time", "sample", "series"))
  
  
  
  # thin: dt = 1
  series_sim <- thin(series_sim, time = "time")
  
  return(series_sim)
  
}) %>% do.call(rbind, .)

simulated_densities <- lapply(1:N, function(i) {
  
  # number of samples
  n_samples <- shoji_cusp_samples[[i]]@stan_args[[1]]$iter - shoji_cusp_samples[[i]]@stan_args[[1]]$warmup
  
  # samples indices
  sample_ind <- sample(1:n_samples, n_sim)
  
  par_samples <- rstan::extract(shoji_cusp_samples[[i]], "theta")[[1]][sample_ind, , ]
  
  
  
  # simulate series
  ticks <- 0.1
  series_sim <- lapply(1:n_sim, function(s) {
    
    shoji_generator(y0 = rnorm(1, 5, 1), times = seq(from=0, to=n_points, by = ticks),
                    r = 1,
                    # r = par_samples[j, 5],
                    alpha = par_samples[s, 1],
                    beta = par_samples[s, 2],
                    lambda = par_samples[, 3],
                    epsilon = par_samples[s, 4],
                    seed = sample(1:1000, 1)) %>% 
      as.data.frame() %>%
      cbind(time = seq(from=0, to=n_points, by = ticks)) %>%
      mutate(sample = s)
  }) %>% 
    do.call(rbind, .) %>% 
    cbind(series = i) %>% 
    set_colnames(c("value", "time", "sample", "series"))
  
  
  
  # thin: dt = 1
  series_sim <- thin(series_sim, time = "time")
  
  return(series_sim)
  
}) %>% do.call(rbind, .)

# Plot
samples_p <- ggplot() +
  stat_density(data = simulated_densities,
               aes(x = value, group = sample),
               geom = "line",
               position = "identity",
               color = "red",
               alpha = 0.5) +
  stat_density(data =  df_shoji %>%
                 set_colnames(c(1:N %>% as.character, "time")) %>%
                 melt(id.vars = "time") %>% mutate(series = variable),
               aes(x = value),
               geom = "line", position = "identity") +
  facet_wrap(~series, scales = "free") +
  labs(subtitle = "David et al data (CLR) and posterior samples; 100 time points, WITHOUT r parameter") 



samples_p_with_r <- ggplot() +
  stat_density(data = simulated_densities_with_r,
               aes(x = value, group = sample),
               geom = "line",
               position = "identity",
               color = "red",
               alpha = 0.5) +
  stat_density(data =  df_shoji %>%
                 set_colnames(c(1:N %>% as.character, "time")) %>%
                 melt(id.vars = "time") %>% mutate(series = variable),
               aes(x = value),
               geom = "line", position = "identity") +
  facet_wrap(~series, scales = "free") +
  labs(subtitle = "David et al data (CLR) and posterior samples; 100 time points, WITH r parameter") 


plot_grid(samples_p, samples_p_with_r, ncol = 1)