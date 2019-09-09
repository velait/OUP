oup_model <- stan_model("stan_models/oup_single_transition.stan")
ar1_model <- stan_model("stan_models/ar1.stan")
ar1_missing_model <- stan_model("stan_models/ar1_missing.stan")


times <- seq(from = 0, to = 100, by = 1)

## AR(1) ******************** ####

sigma <- .1
lambda <- .75
epsilon <- 0


ar1_series <- ar1_generator(times = times, lambda = lambda, sigma = sigma, epsilon = epsilon)


ar1_series %>% plot


## AR(1) with NA ************ ####

## Remove randomly p%
na_prop <- .5
na_index <- sample(1:length(ar1_series), na_prop*length(ar1_series), replace = FALSE)
obs_index <- which(!(1:length(ar1_series) %in% na_index))

ar1_missing_samples <- sampling(ar1_missing_model,
                        list(N_obs = length(obs_index),
                             N_mis = length(na_index), 
                             ii_obs = obs_index,
                             ii_mis = na_index,
                             Y_obs = ar1_series[obs_index]), 
                        iter = 2000, chains = 1)

ar1_missing_samples




## OUP

oup_series <- oup_generator(N = length(times), times = times, y0 = 0, mu = 0, lambda = lambda, sigma = sigma)



oup_samples <- sampling(oup_model, list(T = length(times), time = times, Y = oup_series), 
                        iter = 2000, chains = 1)

oup_samples