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
iter <- 1000


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
    seed <- n_series*n_observations
    set.seed(seed)
    lambda <- oup_invG_lambda(n_series, shape = lambda_mean, scale = lambda_sd)
    sigma <- rnorm(n_series, mean = sigma_mean, sd = sigma_sd)
    mu <- rnorm(n_series, mean = mu_mean, sd = mu_sd)
    
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
    
    non_pooled_samples <- sampling(non_pooled_student_t_oup,
                                   diff_compare_data,
                                   iter=iter,
                                   chains=chains, 
                                   init=1)
    
    partially_pooled_samples <- sampling(hierarchical_student_t_oup,
                                         diff_compare_data,
                                         iter=iter,
                                         chains=chains,
                                         init=1)
    
    models <- list(pooled_samples, 
                   non_pooled_samples,
                   partially_pooled_samples) %>% 
      set_names(c("pooled", 
                  "non_pooled",
                  "partially"))
    
    
    
    
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

write.csv(results, file = "results/results.csv")

# Plot ******************************************************* ####


  
## Get modelwise estimate plots
modelwise_estimate_panels <- lapply(parameters, function(par) {
  
  df <- results %>% 
    filter(parameter == par) %>% 
    mutate(index = index %>% as.numeric)

  p <- df %>% 
    ggplot()  + 
    geom_line(data = df, aes(x=index, y=simulation_value)) +
    geom_errorbar(data = df %>% filter(model != "pooled"),
                  aes(x=index, ymin=IQR50_lower, ymax=IQR50_upper, color = model, width = 0.1),
                  position = "dodge") +
    facet_wrap(c( "n_observations", "n_series"), labeller = "label_both") +
    labs(x="Series", y="Posterior estimate", title=par) +
    guides(color=guide_legend("Model")) +
    scale_color_manual(values = c("#999999", "#E69F00"))+
    theme_bw()
  
  # Add complete pooling estimate
  p <- p + geom_line(data = df %>% filter(model == "pooled"), aes(x = index, y = median), linetype ="dashed")
  
  return(p)
  
})



# Minimal IQR panel
minimal_IQR_panel <- lapply(parameters, function(par) {
  
  results %>% 
    filter(parameter == par) %>% 
    ggplot(aes(x = n_observations, y = IQR_min, color = index)) +
    geom_line() +
    facet_wrap(c( "model", "n_series"), labeller = "label_both", nrow = 3) +
    theme_bw() +
    scale_color_brewer(palette = 7)
  
})


minimal_average_IQR_panel <- lapply(parameters, function(par) {
  
  results %>% 
    filter(parameter == par) %>% 
    ggplot(aes(x = n_observations, y = mean_IQR, color = as.factor(n_series))) +
    geom_line() +
    facet_wrap(c( "model"), nrow =1) +
    theme_bw() +
    scale_color_brewer(palette = 7)
  
})


average_sd_panel <- lapply(parameters, function(par) {
  
  results %>% 
    filter(parameter == par) %>% 
    ggplot(aes(x = n_observations, y = mean_sd, color = as.factor(n_series))) +
    geom_line() +
    facet_wrap(c( "model"), nrow =1) +
    theme_bw() +
    scale_color_brewer(palette = 7)
  
})