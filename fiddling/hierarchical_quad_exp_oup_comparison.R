## stan models ************************************************************* ####

# generators
oup_simulator <- stan_model("fiddling/simulate_gp_OUP.stan")
gp_simulator <- stan_model("fiddling/simulate_gp_exp_quad.stan")

# simulate underlying process for the quantiles
stan_simulator <- stan_model(file='fiddling/simu_gauss_dgp.stan')

# fitters
exp_quad_fitter <- stan_model("fiddling/fit_gp.stan")
oup_fitter <- stan_model("fiddling/fit_oup.stan")
fitter_list <- list(oup = oup_fitter, exp_quad = exp_quad_fitter)


# complete, partial and non-pooled versions
pooled_gp_oup <- stan_model("fiddling/stan_models/fit_pooled_oup.stan")
non_pooled_gp_oup <- stan_model("fiddling/stan_models/fit_non_pooled_oup.stan")
hierarchical_pooled_gp_oup <- stan_model("fiddling/stan_models/fit_hierarchical_oup.stan")



## parameters ************************************************************** ####


# alpha = kappa^2
# rho = (inv_lambda)^2

rho_shape <- 15.4  # length scale
rho_rate <- 79.4

alpha_mean <- 1    # marginal variance
alpha_sd <- 0.5

sigma <- 0.01 # measurement error

parameters <- c("rho", "alpha", "sigma")

## Data ******************************************************************** ####

series_grid <- c(5, 10, 25, 50)
observation_grid <- c(5, 10)

gp_data_set <- lapply(series_grid, function(n_series) {
  print(paste0(n_series, " series"))
  
  data <- list()
  
  obs_data <- lapply(observation_grid, function(n_observations) {
    print(paste0(n_observations, " observations"))

  # Data
  seed <- n_series*n_observations + 1
  set.seed(seed)
  
  rho_true <- rinvgamma(n_series, shape = rho_shape, rate = rho_rate)       # marginal variance (~ kappa^.5)
  alpha_true <- restricted_rnorm(n_series, alpha_mean, alpha_sd, lower = 0)        # length scale (~ inv_lambda^.5)
  sigma_true <- rep(sigma, n_series)                              # measurement error
  
  dat <- genera_gp_set_stan(n_series = n_series,
                                            alpha = alpha_true,
                                            rho = rho_true, 
                                            sigma = sigma_true,
                                            intervals = 1:n_observations,
                                            stan_model = oup_simulator,
                                            seed = seed)
  
  return(dat)


  }) %>% set_names(as.character(observation_grid))
  
  data[[as.character(n_series)]] <- obs_data
  
  
}) %>% set_names(as.character(series_grid))


# Loop over series and observations ************************** ####

loop_results <- list()

for(n_series in series_grid) {
  print(n_series)
  
  results <- list()
  
  for(n_observations in observation_grid) {
    print(n_observations)
    
    stan_data <- gp_data_set[[as.character(n_series)]][[as.character(n_observations)]]
    
    # Stan samples ************************************************
    pooled_samples <- sampling(pooled_gp_oup,
                               stan_data,
                               iter=iter,
                               chains=chains,
                               init=1)
    
    non_pooled_samples <- sampling(non_pooled_gp_oup,
                                   stan_data,
                                   iter=iter,
                                   chains=chains,
                                   init=1)
    
    partially_pooled_samples <- sampling(hierarchical_pooled_gp_oup,
                                         stan_data,
                                         iter=iter,
                                         chains=chains,
                                         init=1)
    
    # models <- list(pooled_samples) %>% 
    #   set_names(c("pooled"))
    
    models <- list(pooled_samples,
                   non_pooled_samples,
                   partially_pooled_samples) %>%
      set_names(c("pooled",
                  "non_pooled",
                  "partially"))
    
    
    
    
    # Results ****************************************************
    results[[n_observations %>% as.character()]] <- lapply(names(models), function(model) {
      
      lapply(parameters, function(p) {
        
        p_values <- stan_data[[grep(p, names(stan_data))]]
        
        get_IQRs(models[[model]] , p, p_values)
        
      }) %>%
        do.call(rbind, .) %>% 
        mutate(model = model)
      
      
      
    }) %>%
      do.call(rbind, .)
    
    # hyperparameter_results[[n_observations %>% as.character()]]  <- lapply(names(models), function(model) {
    #   
    #   lapply(hyper_parameters, function(p) {
    #     
    #     get_hyperparamters(models[[model]] , p, get(p))
    #     
    #   }) %>%
    #     do.call(rbind, .) %>% 
    #     mutate(model = model)
    #   
    # }) %>%
    #   do.call(rbind, .)
    
    
  }
  
  results <- results %>% do.call(rbind, .)
  
  loop_results[[n_series %>% as.character()]] <- results
}

results <- loop_results %>% 
  do.call(rbind, .) %>% 
  set_rownames(NULL)

for(col in c("mean_sd", "sd", "simulation_value", "IQR50_lower", "IQR50_upper", "mean_IQR", "IQR_min", "mean", "median")) {
  results[, col] <- results[, col] %>% as.character() %>% as.numeric()
}

write.csv(results, file = "fiddling/results.csv")
# results <- read.csv(file = "results/results.csv")


       


# Plot ******************************************************* ####


## Modelwise scatter plots: posterior mean vs. simulation value; length scale/alpha ratios
modelwise_ratio_plots <- lapply(c("partially", "pooled", "non_pooled"), function(x) {
  
  
  df <- results %>% 
    filter(model == x)
  
  rho_values <- df %>% filter(parameter == "rho")
  rho_values <- rbind(data.frame(rho = rho_values$mean,
                                 index = rho_values$index,
                                 type = rep("mean", nrow(rho_values)),
                                 n_series = rho_values$n_series,
                                 n_observations = rho_values$n_observations),
                      data.frame(rho = rho_values$simulation_value,
                                 index = rho_values$index,
                                 type = rep("simulation_value", nrow(rho_values)),
                                 n_series = rho_values$n_series,
                                 n_observations = rho_values$n_observations))
  
  alpha_values <- df %>% filter(parameter == "alpha")
  alpha_values <- rbind(data.frame(alpha = alpha_values$mean),
                        data.frame(alpha = alpha_values$simulation_value))
  
  df <- cbind(rho_values, alpha_values)
  
  
  ggplot() +
    geom_line(data = df, aes(x = alpha, y = rho, group = index)) +
    geom_point(data = df, aes(x = alpha, y= rho, color = type)) +
    facet_wrap(c("n_observations", "n_series"), labeller = "label_both", ncol=length(series_grid)) +
    theme_bw() +
    labs(title = paste0("Pooling: ", x),
         x = "Alpha (marginal standard deviation)", 
         y = "Rho (length scale)")
  
})

modelwise_ratio_comparison_plots <- lapply(c("partially", "pooled", "non_pooled"), function(x) {
  
  
  df <- results %>% 
    filter(model == x)
  
  rho_values <- df %>% filter(parameter == "rho")
  
  rho_values <- rbind(data.frame(rho = rho_values$mean,
                                 index = rho_values$index,
                                 type = rep("mean", nrow(rho_values)),
                                 n_series = rho_values$n_series,
                                 n_observations = rho_values$n_observations),
                      data.frame(rho = rho_values$simulation_value,
                                 index = rho_values$index,
                                 type = rep("simulation_value", nrow(rho_values)),
                                 n_series = rho_values$n_series,
                                 n_observations = rho_values$n_observations))
  
  alpha_values <- df %>% filter(parameter == "alpha")
  alpha_values <- rbind(data.frame(alpha = alpha_values$mean),
                        data.frame(alpha = alpha_values$simulation_value))
  
  df <- cbind(rho_values, alpha_values)
  
  sim_df <- df %>% filter(type == "simulation_value") %>% 
    mutate(ratio = rho/alpha)
  est_df <- df %>% filter(type == "mean") %>% 
    mutate(ratio = rho/alpha)
  
  df <- data.frame(sim_ratio = sim_df$ratio,
                   est_ratio = est_df$ratio,
                   n_series = est_df$n_series, 
                   n_observations = est_df$n_observations)
  
  p <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    geom_point(data = df, aes(x = sim_ratio, y = est_ratio)) +
    facet_wrap(c("n_observations", "n_series"), labeller = "label_both", ncol=length(series_grid)) +
    theme_bw() +
    labs(title = paste0("Pooling: ", x),
         x = "True Ratio", 
         y = "Ratio of esimated means",
         subtitle = "Length scale / variation ratio; rho/alpha")
  
  p
  
}) %>% 
  set_names(c("partially", "pooled", "non_pooled"))
    
    
## Get modelwise estimate plots
modelwise_estimate_panels <- lapply(parameters, function(par) {
  
  df <- results %>% 
    filter(parameter == par) %>% 
    group_by(model, n_observations, n_series) %>% 
    mutate(rank = rank(simulation_value))
  
  p <- df %>% 
    ggplot()  + 
    geom_line(data = df, aes(x=rank, y=simulation_value)) +
    geom_errorbar(data = df %>% filter(model != "pooled"),
                  aes(x=rank, ymin=IQR50_lower, ymax=IQR50_upper, color = model, width = 0.25),
                  position = "dodge") +
    facet_wrap(c( "n_observations", "n_series"), labeller = "label_both", ncol = length(series_grid)) +
    labs(x="Series", y="Posterior estimate", title=par) +
    guides(color=guide_legend("Model")) +
    scale_color_manual(values = c("#999999", "#E69F00"))+
    theme_bw()
  
  # Add complete pooling estimate
  p <- p + geom_line(data = df %>% filter(model == "pooled"), aes(x = rank, y = median), linetype ="dashed")
  
  return(p)
  
}) %>% set_names(parameters)

## Parameterwise poesterior means for pooled model
pooled_parameters_posterior_means <- lapply(parameters, function(par) {
  
  results <- results %>%  mutate(mean = as.numeric(as.character(mean)))
  
  df <- results %>% 
    filter(parameter == par) %>% 
    group_by(model, n_observations, n_series) %>% 
    mutate(rank = rank(simulation_value), inv_mean = 1/mean)
  
  
  p <- df %>%
    filter(model == "pooled") %>%
    ggplot(aes(x = n_series, y = mean, color = as.factor(n_observations))) +
    geom_line() +
    labs(title = par) +
    theme(legend.position = "top")
  
})
plot_grid(pooled_parameters_posterior_means[[1]],
          pooled_parameters_posterior_means[[2]], 
          pooled_parameters_posterior_means[[3]], ncol = 3)

## Ratio plot


# Minimal IQR panel
minimal_IQR_panel <- lapply(parameters, function(par) {
  
  df <- results %>% 
    filter(parameter == par) %>% 
    select(index, simulation_value, IQR_min, n_series, n_observations, model)
  
  df %>% ggplot(aes(x = model, y = IQR_min)) +
    geom_point() +
    facet_wrap(c("n_series", "n_observations"), labeller = "label_both") +
    geom_line(data = df, aes(x=model, y=IQR_min, group = index)) +
    theme_bw()
  
})%>% set_names(parameters)


minimal_average_IQR_panel <- lapply(parameters, function(par) {
  
  results %>% 
    filter(parameter == par) %>% 
    ggplot(aes(x = n_observations, y = mean_IQR, color = as.factor(n_series))) +
    geom_line() +
    facet_wrap(c( "model"), nrow =1) +
    theme_bw() +
    labs(y = "Mean min IQR") 
  
}) %>% set_names(parameters)


average_sd_panel <- lapply(parameters, function(par) {
  
  results %>% 
    filter(parameter == par) %>% 
    ggplot(aes(x = n_observations, y = mean_sd, color = as.factor(n_series))) +
    geom_line() +
    facet_wrap(c( "model"), nrow =1) +
    theme_bw()
  
  
}) %>% set_names(parameters)
