# Cusp + observation model


## Cusp
N <- 9
puolet <- 5

r <- runif(N, .5, 2)
alpha <- c(0, 0, 0, 0, 0,
           -0.5, -0.5, -1, 1)
beta <- c(2, -2, 1.5, -1.5, -1,
          1, -1,  0, 0)
lambda <- rnorm(N, 0, .25)
log_epsilon <- rnorm(N, 0, 0.5)
epsilon <- exp(log_epsilon)


grid <- seq(from = 0, to = 100, by = 0.1)
seed = sample(1:1000, 1)


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
  labs(title = "Latent cusp process", subtitle = "100 observations, dt = 1")


p




## Softmax for probabilities 
## Softmax works for negative values as well
theta_df <- lapply(1:nrow(df_shoji), function(i) {
  
  df <- df_shoji[i, ]
  
  x <- df[1:(length(df) - 1)] %>% 
    softMax
  
  cbind(x, df[length(df)])
  
}) %>% do.call(rbind, .)

p <- theta_df %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value)) +
  geom_line() +
  facet_wrap(~variable) +
  labs(title = "Cusp + softmax", subtitle = "100 observations, dt = 1")


p

softmax_id <- function(x) {

  x <- c(x, 0)
  
  softMax(x);
}


## Multinomial sampling


obs_df <- apply(theta_df[, -ncol(theta_df)], 1, function(i) {
  rmultinom(1, 1e5, i)
}) %>% t %>%
  as.data.frame() %>%
  cbind(x = theta_df$x) %>%
  set_colnames(c(paste0("y", 1:N), "x"))


p <- obs_df %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value)) +
  geom_line() +
  facet_wrap(~variable) +
  labs(title = "Cusp + softmax + multinomial", subtitle = "100 observations, dt = 1")


p


## Inference with Stan ******************************** ####

## Uninformative priors ******************************
stan_data <- list(N_obs = nrow(obs_df),
                  N_OTUs = N,
                  x = obs_df$x,
                  y = t(select(obs_df, -x)))


# cusp_plus_obsevations <- stan_model("cusp/cusp_plus_observation.stan")
cusp_plus_obsevations_with_r <- stan_model("cusp/cusp_plus_observation_with_r.stan")

# cusp_plus_obsevations_samples <- sampling(cusp_plus_obsevations, stan_data,
#                                           iter = 1000, 
#                                           chains = 1, 
#                                           control = list(adapt_delta = 0.9))

cusp_plus_obsevations_with_r_samples <- sampling(cusp_plus_obsevations_with_r, 
                                                 stan_data,
                                                 iter = 1000, 
                                                 chains = 1, 
                                                 control = list(adapt_delta = 0.9, 
                                                                max_treedepth = 11))


## Informative priors *********************************

# MLE estimates for stationary densities
# 
# for(i in 1:N) {
#   
#   x <- 
#   cusp_mLL <- function(par) {
#     
#     # alpha <- par[1]
#     # beta <- par[2]
#     # lambda <- par[3]
#     # epsilon <- par[4]
#     
#     # Normalizing factor
#     M <- integrate(
#       function(x, alpha = par[1], beta = par[2], lambda = par[3], epsilon = par[4]) {
#         exp((alpha*(x - lambda) + .5*beta*(x - lambda)^2 - .25*(x - lambda)^4 ) / epsilon)
#       },
#       -Inf, Inf)$value
#     
#     # Log likelihood
#     LL <- (par[1]*(x - par[3]) + .5*par[2]*(x - par[3])^2 - .25*(x - par[3])^4 ) / par[4] - log(M)
#     
#     return((-sum(LL)))
#   }
#   ml_est <- optim(c(0, 0, 0, 1),
#                   cusp_mLL, method = "L-BFGS-B",
#                   upper = c(20, 20, 20, 20),
#                   lower = c(-20, -20, -20, 0.025),
#                   control = list(maxit = 1000))
#   
#   
# }
# 
# 
# 
# 
# 
# 
# 
# prior_parameters <- matrix(NA, 4, 2)


## Results ******************************************** ####

# # get results
# results_df <- lapply(1:4, function(x) get_stan_results(cusp_plus_obsevations_samples, paste0("theta\\[.+,", x,"\\]"))) %>% do.call(rbind, .) %>% 
#   cbind(parameter = rep(c("alpha", "beta", "lambda", "epsilon"), each = N),
#         index = factor(rep(1:N, 4)),
#         true = c(alpha, beta, lambda, epsilon))

results_with_r_df <- lapply(1:5, function(x) get_stan_results(cusp_plus_obsevations_with_r_samples, paste0("theta\\[.+,", x,"\\]"))) %>% do.call(rbind, .) %>% 
  cbind(parameter = rep(c("alpha", "beta", "lambda", "epsilon", "r"), each = N),
        index = factor(rep(1:N, 5)),
        true = c(alpha, beta, lambda, epsilon, r))


# plot estimates and true values
# results_df %>% 
#   ggplot() + 
#   geom_point(aes(x = index, y = true), shape = 4) + 
#   geom_errorbar(aes(x = index, ymax = upper97.5, ymin = lower2.5), width = .2) +
#   facet_wrap(~parameter)
# 

results_with_r_df %>% 
  ggplot() + 
  geom_point(aes(x = index, y = true), shape = 4) + 
  geom_errorbar(aes(x = index, ymax = upper97.5, ymin = lower2.5), width = .2) +
  facet_wrap(~parameter)


# compare true, estimated and data densities

x <- seq(-5, 5, length.out = 100)
true_d <- lapply(1:N, function(i) {
  
  d <- cc_density(x, r[i], alpha[i], beta[i], lambda[i], epsilon[i])
  
  data.frame(x, d, r[i], alpha[i], beta[i], lambda[i], epsilon[i], "true", i) %>% 
    set_colnames(c("x", "y", "r", "alpha", "beta", "lambda", "epsilon", "type", "index"))
  
}) %>% do.call(rbind, .)


# map_d <- lapply(1:N, function(i) {
#   
#   rr <- results_df %>%
#     filter(parameter == "r", index == i) %>% 
#     pull(mode)
#   
#   aa <- results_df %>%
#     filter(parameter == "alpha", index == i) %>% 
#     pull(mode)
#   
#   bb <- results_df %>%
#     filter(parameter == "beta", index == i) %>% 
#     pull(mode)
#   
#   ee <- results_df %>%
#     filter(parameter == "epsilon", index == i) %>% 
#     pull(mode)
#   
#   ll <- results_df %>%
#     filter(parameter == "lambda", index == i) %>% 
#     pull(mode)
#   
#   d <- cc_density(x, rr, aa, bb, ll, ee)
#   
#   data.frame(x, d, r[i], alpha[i], beta[i], lambda[i], epsilon[i], "true", i) %>% 
#     set_colnames(c("x", "y", "r", "alpha", "beta", "lambda", "epsilon", "type", "index"))
#   
# }) %>% do.call(rbind, .)
map_with_r_d <- lapply(1:N, function(i) {
  
  rr <- results_with_r_df %>%
    filter(parameter == "r", index == i) %>% 
    pull(mode)
  
  aa <- results_with_r_df %>%
    filter(parameter == "alpha", index == i) %>% 
    pull(mode)
  
  bb <- results_with_r_df %>%
    filter(parameter == "beta", index == i) %>% 
    pull(mode)
  
  ee <- results_with_r_df %>%
    filter(parameter == "epsilon", index == i) %>% 
    pull(mode)
  
  ll <- results_with_r_df %>%
    filter(parameter == "lambda", index == i) %>% 
    pull(mode)
  
  d <- cc_density(x, rr, aa, bb, ll, ee)
  
  data.frame(x, d, r[i], alpha[i], beta[i], lambda[i], epsilon[i], "true", i) %>% 
    set_colnames(c("x", "y", "r", "alpha", "beta", "lambda", "epsilon", "type", "index"))
  
}) %>% do.call(rbind, .)


n_sim <- 20

# simulated_densities <- lapply(1, function(i) {
#   
#   # number of samples
#   n_samples <- cusp_plus_obsevations_samples@stan_args[[1]]$iter - cusp_plus_obsevations_samples@stan_args[[1]]$warmup
#   
#   # samples indices
#   sample_ind <- sample(1:n_samples, n_sim)
#   
#   par_samples <- rstan::extract(cusp_plus_obsevations_samples, "theta")[[1]][sample_ind, , ]
#   
#   
#   
#   
#   # simulate series
#   series_sim <- lapply(1:n_sim, function(s) {
#     lapply(1:N, function(j) {
#       shoji_generator(y0 = rnorm(1, 0, 1), times = seq(from=stan_data[['x']][1], to=stan_data[['x']][length(stan_data[['x']])], by = 0.1),
#                       r = 1,
#                       # r = par_samples[1, j, 5],
#                       alpha = par_samples[1, j, 1],
#                       beta = par_samples[1, j, 2],
#                       lambda = par_samples[1, j, 3],
#                       epsilon = par_samples[1, j, 4],
#                       seed = sample(1:5000, 1))
#     }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame() %>%
#       cbind(time = seq(from=stan_data[['x']][1], to=stan_data[['x']][length(stan_data[['x']])], by = 0.1)) %>%
#       set_colnames(c(paste0("y", 1:N), "x")) %>% 
#       mutate(sample = s)
#   }) %>% 
#     do.call(rbind, .)
#   
#   
#   
#   # thin: dt = 1
#   series_sim <- thin(series_sim)
#   series_sim_no_time <- series_sim %>% select(-c(x, sample))
#   
#   
#   # Softmax at each time step
#   series_sim_softmax <- lapply(1:nrow(series_sim), function(r) {
#     softMax(series_sim_no_time[r, ])
#   }) %>%
#     do.call(rbind, .) %>% 
#     cbind(., series_sim %>%
#             select(c(x, sample)))
#   
#   
#   series_sim_softmax_no_time <- series_sim_softmax %>% select(-c(x, sample))
#   
#   
#   # Multinomial samping; read counts from data
#   
#   series_sim_multin <- lapply(1:nrow(series_sim_softmax_no_time), function(r) {
#     rmultinom(n = 1, size = 50000, prob = series_sim_softmax_no_time[r, ]) %>% 
#       t
#   }) %>% 
#     do.call(rbind, .) %>% 
#     cbind(., series_sim %>%
#             select(c(x, sample)))
#   
#   
#   rm(series_sim_softmax_no_time, series_sim_no_time)
#   
#   return(list(cusp = series_sim, softmax = series_sim_softmax, multin = series_sim_multin))
#   
# })
# map_simulated_densities <- lapply(1, function(i) {
#   
#   
#   # simulate series
#   series_sim <- lapply(1:10, function(s) {
#     lapply(1:N, function(j) {
#       
#       map_pars <- results_df %>%
#         filter(index == j) %>%
#         select(mode, parameter)
#       
#       shoji_generator(y0 = rnorm(1, 0, 1), times = seq(from=0, to=100, by = 0.1),
#                       r = 1,
#                       # map_pars[5, 1],
#                       alpha = map_pars[1, 1],
#                       beta = map_pars[2, 1],
#                       lambda = map_pars[3, 1],
#                       epsilon = map_pars[4, 1],
#                       seed = sample(1:1000, 1))
#     }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame() %>%
#       cbind(time = seq(from=0, to=100, by = 0.1)) %>%
#       set_colnames(c(paste0("y", 1:N), "x")) %>% 
#       mutate(sample = s)
#   }) %>% 
#     do.call(rbind, .)
#   
#   
#   
#   # thin: dt = 1
#   series_sim <- thin(series_sim)
#   series_sim_no_time <- series_sim %>% select(-c(x, sample))
#   
#   
#   # Softmax at each time step
#   series_sim_softmax <- lapply(1:nrow(series_sim), function(r) {
#     softMax(series_sim_no_time[r, ])
#   }) %>%
#     do.call(rbind, .) %>% 
#     cbind(., series_sim %>%
#             select(c(x, sample)))
#   
#   
#   series_sim_softmax_no_time <- series_sim_softmax %>% select(-c(x, sample))
#   
#   
#   # Multinomial samping; read counts from data
#   series_sim_multin <- lapply(1:nrow(series_sim_softmax_no_time), function(r) {
#     rmultinom(n = 1, size = 50000, prob = series_sim_softmax_no_time[r, ]) %>% 
#       t
#   }) %>% 
#     do.call(rbind, .) %>% 
#     cbind(., series_sim %>%
#             select(c(x, sample)))
#   
#   
#   rm(series_sim_softmax_no_time, series_sim_no_time)
#   
#   return(list(cusp = series_sim, softmax = series_sim_softmax, multin = series_sim_multin))
#   
# })


simulated_densities_with_r <- lapply(1, function(i) {
  
  # number of samples
  n_samples <- cusp_plus_obsevations_with_r_samples@stan_args[[1]]$iter - cusp_plus_obsevations_with_r_samples@stan_args[[1]]$warmup
  
  # samples indices
  sample_ind <- sample(1:n_samples, n_sim)
  
  par_samples <- rstan::extract(cusp_plus_obsevations_with_r_samples, "theta")[[1]][sample_ind, , ]
  
  
  
  
  # simulate series
  series_sim <- lapply(1:n_sim, function(s) {
    lapply(1:N, function(j) {
      shoji_generator(y0 = rnorm(1, 0, 1), times = seq(from=stan_data[['x']][1], to=stan_data[['x']][length(stan_data[['x']])], by = 0.1),
                      # r = 1,
                      r = par_samples[1, j, 5],
                      alpha = par_samples[1, j, 1],
                      beta = par_samples[1, j, 2],
                      lambda = par_samples[1, j, 3],
                      epsilon = par_samples[1, j, 4],
                      seed = sample(1:5000, 1))
    }) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      cbind(time = seq(from=stan_data[['x']][1], to=stan_data[['x']][length(stan_data[['x']])], by = 0.1)) %>%
      set_colnames(c(paste0("y", 1:N), "x")) %>% 
      mutate(sample = s)
  }) %>% 
    do.call(rbind, .)
  
  
  
  # thin: dt = 1
  series_sim <- thin(series_sim)
  series_sim_no_time <- series_sim %>% select(-c(x, sample))
  
  
  # Softmax at each time step
  series_sim_softmax <- lapply(1:nrow(series_sim), function(r) {
    softMax(series_sim_no_time[r, ])
  }) %>%
    do.call(rbind, .) %>% 
    cbind(., series_sim %>%
            select(c(x, sample)))
  
  
  series_sim_softmax_no_time <- series_sim_softmax %>% select(-c(x, sample))
  
  
  # Multinomial samping; read counts from data
  
  series_sim_multin <- lapply(1:nrow(series_sim_softmax_no_time), function(r) {
    rmultinom(n = 1, size = 50000, prob = series_sim_softmax_no_time[r, ]) %>% 
      t
  }) %>% 
    do.call(rbind, .) %>% 
    cbind(., series_sim %>%
            select(c(x, sample)))
  
  
  rm(series_sim_softmax_no_time, series_sim_no_time)
  
  return(list(cusp = series_sim, softmax = series_sim_softmax, multin = series_sim_multin))
  
})
map_simulated_densities_with_r <- lapply(1, function(i) {
  
  
  # simulate series
  series_sim <- lapply(1:10, function(s) {
    lapply(1:N, function(j) {
      
      map_pars <- results_with_r_df %>%
        filter(index == j) %>%
        select(mode, parameter)
      
      shoji_generator(y0 = rnorm(1, 0, 1), times = seq(from=0, to=100, by = 0.1),
                      # r = 1,
                      map_pars[5, 1],
                      alpha = map_pars[1, 1],
                      beta = map_pars[2, 1],
                      lambda = map_pars[3, 1],
                      epsilon = map_pars[4, 1],
                      seed = sample(1:1000, 1))
    }) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      cbind(time = seq(from=0, to=100, by = 0.1)) %>%
      set_colnames(c(paste0("y", 1:N), "x")) %>% 
      mutate(sample = s)
  }) %>% 
    do.call(rbind, .)
  
  
  
  # thin: dt = 1
  series_sim <- thin(series_sim)
  series_sim_no_time <- series_sim %>% select(-c(x, sample))
  
  
  # Softmax at each time step
  series_sim_softmax <- lapply(1:nrow(series_sim), function(r) {
    softMax(series_sim_no_time[r, ])
  }) %>%
    do.call(rbind, .) %>% 
    cbind(., series_sim %>%
            select(c(x, sample)))
  
  
  series_sim_softmax_no_time <- series_sim_softmax %>% select(-c(x, sample))
  
  
  # Multinomial samping; read counts from data
  series_sim_multin <- lapply(1:nrow(series_sim_softmax_no_time), function(r) {
    rmultinom(n = 1, size = 50000, prob = series_sim_softmax_no_time[r, ]) %>% 
      t
  }) %>% 
    do.call(rbind, .) %>% 
    cbind(., series_sim %>%
            select(c(x, sample)))
  
  
  rm(series_sim_softmax_no_time, series_sim_no_time)
  
  return(list(cusp = series_sim, softmax = series_sim_softmax, multin = series_sim_multin))
  
})



# Transform simulated data to compositional
# simulated_multin <- simulated_densities[[1]][["multin"]] %>% 
#   select(-c(x, sample))
# 
# simulated_multin_compositional <- lapply(1:nrow(simulated_multin), function(i) {
#   
#   simulated_multin[i, ]/rowSums(simulated_multin)[i]
#   
# }) %>%
#   do.call(rbind, .) %>% 
#   cbind(simulated_densities[[1]][["multin"]][, c("x", "sample")])
# 
# map_simulated_multin <- map_simulated_densities[[1]][["multin"]] %>% 
#   select(-c(x, sample))
# 
# map_simulated_multin_compositional <- lapply(1:nrow(map_simulated_multin), function(i) {
#   
#   map_simulated_multin[i, ]/rowSums(map_simulated_multin)[i]
#   
# }) %>%
#   do.call(rbind, .) %>% 
#   cbind(map_simulated_densities[[1]][["multin"]][, c("x", "sample")])
# 
# 
# 
# 


simulated_multin_with_r <- simulated_densities_with_r[[1]][["multin"]] %>% 
  select(-c(x, sample))

simulated_multin_compositional_with_r <- lapply(1:nrow(simulated_multin_with_r), function(i) {
  
  simulated_multin_with_r[i, ]/rowSums(simulated_multin_with_r)[i]
  
}) %>%
  do.call(rbind, .) %>% 
  cbind(simulated_densities_with_r[[1]][["multin"]][, c("x", "sample")])

map_simulated_multin_with_r <- map_simulated_densities_with_r[[1]][["multin"]] %>% 
  select(-c(x, sample))

map_simulated_multin_compositional_with_r <- lapply(1:nrow(map_simulated_multin_with_r), function(i) {
  
  map_simulated_multin_with_r[i, ]/rowSums(map_simulated_multin_with_r)[i]
  
}) %>%
  do.call(rbind, .) %>% 
  cbind(map_simulated_densities_with_r[[1]][["multin"]][, c("x", "sample")])




# Plot
# p_samples_vs_data <- ggplot() + 
#   stat_density(data =  simulated_multin_compositional %>%
#                  melt(id.vars = c("x", "sample"))  %>%
#                  mutate(variable = gsub("y", "", variable)),
#                aes(x = value, group = sample),
#                geom = "line", position = "identity", alpha = 0.25, color = "red") +
#   
#   
#   stat_density(data =  map_simulated_multin_compositional %>%
#                  melt(id.vars = c("x", "sample")) %>% 
#                  mutate(variable = gsub("y", "", variable)),
#                aes(x = value),
#                geom = "line", position = "identity", color = "blue", size = 1) +
#   
#   
#   stat_density(data = theta_df %>% select(-x) %>% set_colnames(as.character(1:N)) %>% melt(),
#                aes(x = value), geom = "line", size = 1) +
#   
#   facet_wrap(~variable, scales = "free") +
#   labs(title = "Simulated; compositional", subtitle = paste0("Black = comp. data; Blue = MAP, simulated, Red = posterior samples, simulated;",  dim(df_shoji)[1],  "time points"))
# 
# 
# p_samples_vs_data_map <- ggplot() +
#   geom_line(data = true_d, aes(x = x, y = y)) +
#   geom_line(data = map_d, aes(x = x, y = y), color = "red") +
#   stat_density(data =df_shoji %>%
#                  melt(id.vars = "x") %>%
#                  mutate(index = rep(1:N, each=nrow(df_shoji))),
#                aes(x = value), geom = "line", linetype = "dashed") +
#   facet_wrap(~index) +
#   labs(title = "Latent", subtitle = "Solid = true; dashed = data; red = MAP estimate")
# 
# 
# plot_grid(p_samples_vs_data, p_samples_vs_data_map)





p_samples_vs_data_with_r <- ggplot() + 
  stat_density(data =  simulated_multin_compositional_with_r %>%
                 melt(id.vars = c("x", "sample"))  %>%
                 mutate(variable = gsub("y", "", variable)),
               aes(x = value, group = sample),
               geom = "line", position = "identity", alpha = 0.25, color = "red") +
  
  
  stat_density(data =  map_simulated_multin_compositional_with_r %>%
                 melt(id.vars = c("x", "sample")) %>% 
                 mutate(variable = gsub("y", "", variable)),
               aes(x = value),
               geom = "line", position = "identity", color = "blue", size = 1) +
  
  
  stat_density(data = theta_df %>% select(-x) %>% set_colnames(as.character(1:N)) %>% melt(),
               aes(x = value), geom = "line", size = 1) +
  
  
  facet_wrap(~variable, scales = "free") +
  labs(title = "Simulated; compositional", subtitle = paste0("Black = comp. data; Blue = MAP, simulated, Red = posterior samples, simulated;",  dim(df_shoji)[1],  "time points"))



p_samples_vs_data_map_with_r <- ggplot() +
  geom_line(data = true_d, aes(x = x, y = y)) +
  geom_line(data = map_with_r_d, aes(x = x, y = y), color = "blue") +
  stat_density(data = df_shoji %>%
                 melt(id.vars = "x") %>%
                 mutate(index = rep(1:N, each=nrow(df_shoji))),
               aes(x = value), geom = "line", linetype = "dashed") +
  facet_wrap(~index, scales = "free") +
  labs(title = "Latent", subtitle = "Solid = true; dashed = data; red = MAP estimate")

# p_with_error <- plot_grid(p_samples_vs_data_with_r, p_samples_vs_data_map_with_r)