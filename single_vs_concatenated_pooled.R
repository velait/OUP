# Compare single long vs. SAME DATA split into shorter series and input into completely pooled model. 


## Data ******************************** ####

# Set parameter values
mu <- 0
lambda <- 0.5
inv_lambda <- 1/lambda
sigma <- 0.25
intervals <- 1:150

parameters <- c("lambda", "mu", "sigma")

# Simulate data
compare_single_vs_pooled_data <- generate_student_set(n_series = 1,
                                  student_df = 5,
                                  mu = mu,
                                  lambda = lambda,
                                  sigma = sigma,
                                  intervals = intervals,
                                  seed = 11235)

compare_single_vs_pooled_data$Y <- compare_single_vs_pooled_data$Y %>% as.matrix(1, 150) %>% t

## Model ******************************* ####
pooled_student_t_oup <- stan_model("stan_models/pooled_student_t_oup.stan")

## Single ****************************** ####

# Sample
single_series_samples <- sampling(pooled_student_t_oup,
         compare_single_vs_pooled_data,
         iter=iter,
         chains=chains,
         init=1)



# Get posterior samples
single_series_posterior <- sapply(parameters,
                    function(x) rstan::extract(single_series_samples, x)) %>%
  set_names(parameters)


# Plot posterior with simulation values
single_posterior_plot <- lapply(parameters, function(x) {
  
  p <- single_series_posterior[[x]] %>% 
    as.data.frame() %>%
    ggplot(aes(x = .)) +
    geom_density() +
    geom_vline(xintercept = eval(parse(text = x)), linetype = "dashed") +
    labs(title = x)
  
  
  
}) %>% set_names(parameters)

single_series_posterior_plot_grid <- plot_grid(single_posterior_plot[[1]],
                                               single_posterior_plot[[2]],
                                               single_posterior_plot[[3]])


## Completely pooled ******************* ####

# Butcher data into 10 short series
compare_single_vs_pooled_data$Y <- matrix(compare_single_vs_pooled_data$Y, 10, 15, byrow = TRUE)
compare_single_vs_pooled_data$time <- 1:15
compare_single_vs_pooled_data$T <- 15
compare_single_vs_pooled_data$N <- 10

# Sample
pooled_series_samples <- sampling(pooled_student_t_oup,
                                   compare_single_vs_pooled_data,
                                   iter=iter,
                                   chains=chains,
                                   init=1)


# Get posterior samples
pooled_posterior <- sapply(parameters,
                    function(x) rstan::extract(pooled_series_samples, x)) %>%
  set_names(parameters)


pooled_posterior_plot <- lapply(parameters, function(x) {
  
  p <- pooled_posterior[[x]] %>% 
    as.data.frame() %>%
    ggplot(aes(x = .)) +
    geom_density() +
    geom_vline(xintercept = eval(parse(text = x)), linetype = "dashed") +
    labs(title = x)
  
  
  
}) %>% set_names(parameters)

pooled_posterior_plot_grid <- plot_grid(pooled_posterior_plot[[1]],
          pooled_posterior_plot[[2]],
          pooled_posterior_plot[[3]])
