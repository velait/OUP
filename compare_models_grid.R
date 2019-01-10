# Data: 
# pooled_samples,
#  non_pooled_samples,
#   partially_pooled_samples,
#     intervals,
#       n_series,
#         lambda, mu, sigma (simulation values)

load(file = "results/samples_5_30.Rds")


# Settings *************************************************** #### 
seed <- 11235

parameters <- c('mu', 'lambda', 'sigma')

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


series_grid <- 2
observation_grid <- 5


# Loop over series and observations ************************* ####

loop_results <- list()

chains <- 2
iter <- 2000

for(n_series in series_grid) {
  print(n_series)
  
  results <- list()
  
  for(n_observations in observation_grid) {
    print(n_observations)
    
    # Data
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
      
    }) %>% do.call(rbind, .)
    
    
  }
  
  results <- results %>% do.call(rbind, .)
  
  loop_results[[n_series %>% as.character()]] <- results
}

results <- loop_results %>% 
  do.call(rbind, .) %>% 
  set_rownames(NULL)

write.csv(results, file = "results/results.csv")

# Plot ******************************************************* ####

results %>% 
  filter(n_series == 2) %>% 
  filter(parameter == "lambda") %>% 
  filter(model == "partially") %>% 
  ggplot(aes(x = n_observations, y = IQR)) +
  geom_point() +
  geom_smooth()



  




