# Single series scan


# simulate underlying process for the quantiles
stan_simulator <- stan_model(file='fiddling/simu_gauss_dgp.stan')

## Data ****************************** ####

rho_grid <- 1:10
alpha_grid <- 1:10
sigma <- .25
n_observations <- 100
n_total_samples <- 501
x_total <- 100 * (0:(N_total - 1)) / (N_total - 1)
seed <- 11235

# Simulate data (including hold out data)
single_data_set <- lapply(rho_grid, function(r) {
  
  res <- lapply(alpha_grid, function(a) {
    
    gp_simulator_data <- list(N = n_total_samples,
                              x = x_total,
                              alpha = a,
                              rho = r,
                              sigma =  sigma_true)
    
    samples <- sampling(oup_simulator,
                        data = gp_simulator_data, 
                        iter=1,
                        chains=1,
                        seed=seed,
                        algorithm="Fixed_param")
    
    samples
      
  }) %>% set_names(alpha_grid %>% as.character())

  
  res
  
}) %>% set_names(rho_grid %>% as.character())


# collect simulated data
single_simulated_data <- lapply(rho_grid, function(r) {
  
  res <- lapply(alpha_grid, function(a) {
    
    stan_model <- single_data_set[[r]][[a]]
    
    fs <- grep("f", rownames(summary(stan_model)$summary), value = T)
    ys <- grep("y", rownames(summary(stan_model)$summary), value = T)
    
    
    df <- data.frame(f = summary(stan_model)$summary[fs, "mean"],
                     y = summary(stan_model)$summary[ys, "mean"],
                     rho = r,
                     alpha = a,
                     x = x_total)
    
  }) %>%  set_names(alpha_grid %>% as.character())
  
  res
  
}) %>% set_names(alpha_grid %>% as.character())
  

# Get the underlying data generating process
single_data_process_plot_data <- lapply(rho_grid, function(r) { 
  
  lapply(alpha_grid, function(a) {
  
  
  
  f_data <- list(sigma=sigma,
                 N=n_total_samples,
                 f=single_simulated_data[[r]][[a]] %>% pull(f))
  
  
  dgp_fit <- sampling(stan_simulator, data=f_data, iter=1000, warmup=0,
                      chains=1, seed=5838298, refresh=1000, algorithm="Fixed_param")
  
  df <- summary(dgp_fit)$summary[, c("mean", "2.5%", "25%", "50%", "75%", "97.5%")] %>% 
    as.data.frame()
  
  df <- df[1:(nrow(df) - 1), ]
  
  df <- df %>% mutate(x = x_total)
  return(df)
}) %>% 
  set_names(alpha_grid %>% as.character)

}) %>% 
  set_names(rho_grid %>% as.character)

# Neat
neat_single_data_process_plot_data <- lapply(rho_grid, function(r) {
  lapply(alpha_grid, function(a) {
    
    df <- cbind(single_data_process_plot_data[[r]][[a]], 
                rho = r,
                alpha = a)
    df
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  set_colnames(c("mean", "low2.5", "low25", "median", "high75", "high97.5", "x", "rho", "alpha"))


# Plot the process
single_process_plots <- lapply(rho_grid, function(r) {
  lapply(alpha_grid, function(a) {
    
   p <-  ggplot(neat_single_data_process_plot_data %>% filter(rho == r, alpha == a)) +
      geom_ribbon(aes(x = x, ymin = low2.5, ymax = high97.5), fill = "darkgreen") +
      geom_ribbon(aes(x = x, ymin = low25, ymax = high75), fill = "chartreuse") +
      geom_line(aes(x = x, y = mean))
    
    return(p)
    
  })
})

# Plot some in a panel
single_process_plots_panel <- plot_grid(single_process_plots[[1]][[1]] + ggtitle("rho = 1; alpha = 1"),
                                        single_process_plots[[1]][[10]]+ ggtitle("rho = 1; alpha = 10"),
                                        single_process_plots[[10]][[1]]+ ggtitle("rho = 10; alpha = 1"),
                                        single_process_plots[[10]][[10]]+ ggtitle("rho = 10; alpha = 10"),
                                        single_process_plots[[5]][[5]]+ ggtitle("rho = 5; alpha = 5"),
                                        single_process_plots[[10]][[5]]+ ggtitle("rho = 10; alpha = 5"), 
                                        ncol = 2)




## Fit ******************************* ####
single_process_results <- lapply(rho_grid, function(r) {
  
  lapply(alpha_grid, function(a) {
    
    dat <- single_simulated_data[[r]][[a]] %>% 
      filter((x%%1==0))
    
    dat <- list(N = nrow(dat), y = dat$y, x = dat$x)
    
    
    samples <- sampling(oup_fitter,
            dat,
            iter=iter,
            chains=chains,
            init=1)
    
    summary(samples)$summary[c("rho", "alpha", "sigma"), c("2.5%", "50%", "97.5%")] %>% 
      as.data.frame() %>% 
      cbind(parameter = c("rho", "alpha", "sigma"), true_value = c(r, a, sigma)) %>% 
      set_colnames(c("lower2.5", "mode", "upper97.5", "parameter", "true_value"))
    
  })
  
})

# save(single_process_results, file = "fiddling/single_series_scan_results.Rdata")
# load(file = "fiddling/single_series_scan_results.Rdata")


## Plots ****************************** ####
single_process_results <- lapply(single_process_results, function(x) {
  x %>% do.call(rbind, .)
}) %>%
  do.call(rbind, .) %>% 
  cbind(index = rep(1:100, each = 3))



df <- single_process_results %>% 
  filter(parameter == "alpha")


  ggplot() +
  geom_point(data = df, aes(x = index, y = mode)) +
  geom_point(data = df, aes(x = index, y = true_value), color = "red")
  
  
  ggplot(data = df) +
    geom_point(aes(x = index, y = mode)) +
    geom_point(aes(x = index, y = true_value), color = "red") +
    geom_errorbar(aes(x = index, ymin = lower2.5, ymax = upper97.5)) +
    labs(title = "alpha 95% estimates; \n red point are true values; \n index 30-40 --> rho = 3, alpha = 1-10 etc.")

  
## length scale/variance ratio plot ************************
true_ratio <- single_process_results[single_process_results$parameter == "rho", "true_value"]/single_process_results[single_process_results$parameter == "alpha", "true_value"]
  
estimate_ratio <- single_process_results[single_process_results$parameter == "rho", "mode"]/single_process_results[single_process_results$parameter == "alpha", "mode"]
  

ratio_df <- data.frame(true_ratio = true_ratio,
                       estimate_ratio = estimate_ratio)


ratio_df %>% ggplot(aes(x = true_ratio, y = estimate_ratio)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Posterior model ratio vs. simulation value ratio")

## 