# Data: 
# pooled_samples,
#  non_pooled_samples,
#   partially_pooled_samples,
#     intervals,
#       n_series,
#         lambda, mu, sigma (simulation values)

# load(file = "results/samples_5_30.Rds")


# Settings *************************************************** #### 
seed <- 11235

parameters <- c('mu', 'lambda', 'sigma')


chains <- 1
iter <- 2000


# Hyperparameters *********************
# lambda ~ inv_gamma
lambda_mean <- 20
lambda_sd <- 10

# sigma ~ Normal
sigma_mean <- .2
sigma_sd <- .1

# mu ~ Normal
mu_mean <- 0
mu_sd <- .25
#**************************************


series_grid <- c(2, 5, 10, 20)
observation_grid <- c(5, 10, 15, 20, 25)

# Stan models ************************************************ ####
pooled_student_t_oup <- stan_model("stan_models/pooled_student_t_oup.stan")
non_pooled_student_t_oup <- stan_model("stan_models/non_pooled_student_t_oup.stan")
hierarchical_student_t_oup <- stan_model("stan_models/hierarchical_student_t_oup.stan")



# Loop over series and observations ************************** ####

loop_results <- list()

for(n_series in series_grid) {
  print(n_series)
  
  results <- list()
  
  for(n_observations in observation_grid) {
    print(n_observations)
    
    # Data
    seed <- n_series*n_observations + 1
    set.seed(seed)
    lambda <- oup_invG_lambda(n_series, shape = lambda_mean, scale = lambda_sd)
    # lambda <- rep(.5, n_series)
    sigma <- rnorm(n_series, mean = sigma_mean, sd = sigma_sd)
    # sigma <- rep(.25, n_series)
    mu <- rnorm(n_series, mean = mu_mean, sd = mu_sd)
    # mu <- rep(0, n_series)
    
    diff_compare_data <- generate_student_set(n_series = n_series,
                                              student_df = 7,
                                              mu = mu,
                                              sigma =  sigma,
                                              lambda = lambda,
                                              intervals = 1:n_observations,
                                              seed = seed)
    

    
    # ************************************************************
    
    # Stan samples ************************************************
    pooled_samples <- sampling(pooled_student_t_oup,
                               diff_compare_data,
                               iter=iter,
                               chains=chains,
                               init=1)
    
    # non_pooled_samples <- sampling(non_pooled_student_t_oup,
    #                                diff_compare_data,
    #                                iter=iter,
    #                                chains=chains, 
    #                                init=1)
    
    # partially_pooled_samples <- sampling(hierarchical_student_t_oup,
    #                                      diff_compare_data,
    #                                      iter=iter,
    #                                      chains=chains,
    #                                      init=1)
    
    models <- list(pooled_samples) %>% 
      set_names(c("pooled"))
    
    # models <- list(pooled_samples, 
    #                non_pooled_samples,
    #                partially_pooled_samples) %>% 
    #   set_names(c("pooled", 
    #               "non_pooled",
    #               "partially"))
    
    
    
    
    # Results ****************************************************
    results[[n_observations %>% as.character()]] <- lapply(names(models), function(model) {
      
      lapply(parameters, function(p) {
        
        get_IQRs(models[[model]] , p, get(p))
        
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

# write.csv(results, file = "results/results.csv")
# 
# results <- read.csv(file = "results/results.csv")
# Plot ******************************************************* ####


  
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
  p <- p + geom_line(data = df %>% filter(model == "pooled"), aes(x = index, y = median), linetype ="dashed")
  
  return(p)
  
}) %>% set_names(parameters)

## Parameterwise psoterior means for pooled model
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
