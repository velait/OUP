# 
# Test gaussian processes on ews data 
# 


## Data **************************** ####

set.seed(2424)
times <- seq(from = 0, to = 100, by = 0.1)
N_series <- 2

# Simulation parameters
r <- 1
K <- 10
cs <- 1
h <- 1
sigma <- 0.1

## Simulate data with CSD
ews_set <- lapply(1:N_series, function(j) {
  
  y <- rep(NA, length(times))
  
  # Initial value
  y[1] <- 8
  
  cs <- 1 + (2.6771 - 1)*(times/max(times))
  
  
  for(i in 2:length(times)) {
    
    dt <- times[(i-1):i]
    
    finer_grid <- seq(from = dt[1], to = dt[2], by = 0.01)
    
    y[i] <- ews_generator(y[i-1], finer_grid, c = cs[i], milstein = T)[length(finer_grid)]
    
    
  }
  
  y
  # return(log(y + 1))
  
  
}) %>%
  do.call(cbind, .) %>%
  cbind(x = times) %>%
  as.data.frame()

ews_residuals <- lapply(1:N_series, function(i) {
  
  series <- ews_set[, i]
  smoothed <- smth(series, window = .1, tails = TRUE)
  res <- series - smoothed
  return(res)

  }) %>% do.call(cbind, .) %>%
  cbind(., ews_set$x) %>% 
  as.data.frame() %>% 
  set_colnames(c(paste0("V", 1:N_series), "x"))



## Observations
non_holes <- times[sample(which((times %% 1) == 0), 50)]
holes <- times[!(times %in% non_holes)]

## Stan Models ********************* ####

oup_model <- stan_model("stan_models/oup_fitter.stan")
se_model <- stan_model("stan_models/square_exp_fitter.stan")
matern_1.5_model <- stan_model("stan_models/matern_1.5_fitter.stan")
matern_2.5_model <- stan_model("stan_models/matern_2.5_fitter.stan")
ar1_model <- stan_model("stan_models/ar1.stan")





## Stan **************************** ####

stan_data <- list(N = length(non_holes), 
                 N_pred = length(holes), 
                 y = ews_residuals %>% filter(x %in% non_holes) %>% pull(V1), 
                 x = ews_residuals %>% filter(x %in% non_holes) %>% pull(x),
                 x_pred = ews_residuals %>% filter(x %in% holes) %>% pull(x))


## Fit models
fits <- lapply(c("se", "oup", paste0("matern_", 1:2, ".5")), function(k) {
  
  stan_data <- list(N = length(non_holes), 
                   N_pred = length(holes), 
                   y = ews_residuals %>% filter(x %in% non_holes) %>% pull(V1), 
                   x = ews_residuals %>% filter(x %in% non_holes) %>% pull(x),
                   x_pred = ews_residuals %>% filter(x %in% holes) %>% pull(x))
  
  
  samples <- sampling(get(paste0(k, "_model")), data = stan_data,
                      chains = chains, iter = iter, control = list(adapt_delta = 0.95))
  
  get_samples <- rstan::extract(samples)
  post_pred <- data.frame(x = stan_data$x_pred,
                          pred_mu = colMeans(get_samples$f_pred), kernel = k) 
  plt_df_rt <- data.frame(x = stan_data$x_pred, f = t(get_samples$f_pred), kernel = k)
  par_est <- data.frame(rho = get_samples$rho,
                        epsilon = get_samples$epsilon, 
                        sigma = get_samples$sigma,
                        lp = get_samples$lp__,
                        kernel = k)
  
  list("post_pred" = post_pred, "plt_df_rt" = plt_df_rt, "par_est" = par_est)
  
})


## Results

post_pred <- lapply(1:length(fits), function(i) {
  fits[[i]][["post_pred"]]
}) %>% do.call(rbind, .)

plt_df_rt <- lapply(1:length(fits), function(i) {
  fits[[i]][["plt_df_rt"]]
}) %>% do.call(rbind, .)

par_est <- lapply(1:length(fits), function(i) {
  fits[[i]][["par_est"]]
}) %>% do.call(rbind, .)



plt_df_rt_melt <- plt_df_rt %>% melt(id.vars = c("x", "kernel"))



# Plot process
# process_p <- 
ggplot(data = ews_residuals %>% filter(x %in% non_holes), aes(x=x, y=V1)) +
  geom_line(data = plt_df_rt_melt, aes(x = x, y = value, group = variable), color = "blue", alpha = 0.25) +
  geom_line(data = ews_residuals, aes(x=x, y=V1, colour = 'True process')) +
  geom_line(data = post_pred, aes(x = x, y = pred_mu, colour = 'Posterior mean function'), size = 1, color = "green") +
  geom_point(size = 2) +
  geom_point(size = 1, color = "white") +
  theme_bw(20) +
  facet_wrap(~kernel, ncol = 2)


# Plot estimates


estimate_p <- par_est %>% 
  melt(id.vars = "kernel") %>%
  filter(variable != "lp") %>% 
  ggplot() + 
  stat_density(aes(x = value, color = kernel), geom = "line", position = "identity") + 
  facet_wrap(~variable, scales = "free")


par_est %>% 
  ggplot(aes(x = kernel, y = lp)) +
  geom_boxplot() + 
  geom_point()








