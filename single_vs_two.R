# Compare long series: one vs. set of two (hierarchical) 

## Data ******************************** ####

# Set parameter values
mu <- 0
lambda <- 0.5
inv_lambda <- 1/lambda
sigma <- 0.25
intervals <- 1:50

parameters <- c("inv_lambda", "lambda", "mu", "sigma")

# Simulate data
single_series <- generate_student_set(n_series = 1,
                                      student_df = 5,
                                      mu = mu,
                                      lambda = lambda,
                                      sigma = sigma,
                                      intervals = intervals,
                                      seed = 11235)
single_series2 <- generate_oup(n = 50, lambda = .5, sigma = .5, mu = 0)


single_series$Y <- single_series2 %>% as.matrix(1, 50) %>% t


two_mu <- rep(mu, 2)
two_lambda <- rep(lambda, 2)
two_inv_lambda <- rep(1/lambda, 2)
two_sigma <- rep(sigma, 2)


two_series <- generate_student_set(n_series = 2,
                                   student_df = 5,
                                   mu = two_mu,
                                   lambda = two_lambda,
                                   sigma = two_sigma,
                                   intervals = intervals,
                                   seed = 11235)


## Models ****************************** ####
pooled_student_t_oup <- stan_model("stan_models/pooled_student_t_oup.stan")
single_student_t_oup <- stan_model("stan_models/single_series_oup.stan")
non_pooled_student_t_oup <- stan_model("stan_models/non_pooled_student_t_oup.stan")
hierarchical_student_t_oup <- stan_model("stan_models/hierarchical_student_t_oup.stan")


## Single ****************************** ####

# Sample
single_series_samples <- sampling(pooled_student_t_oup,
                                  single_series,
                                  iter=iter,
                                  chains=chains,
                                  init=1)


single_series_samples2 <- sampling(single_student_t_oup,
                                  single_series,
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
                                               single_posterior_plot[[3]],
                                               single_posterior_plot[[4]])



## Two ********************************* ####
# Sample
two_series_samples <- sampling(hierarchical_student_t_oup,
                                  two_series,
                                  iter=iter,
                                  chains=chains,
                                  init=1)



# Get posterior samples
two_series_posterior <- sapply(parameters,
                                  function(x) rstan::extract(two_series_samples, x)) %>%
  set_names(parameters)


two_posterior_plot <- lapply(parameters, function(x) {
  
  p <- two_series_posterior[[x]] %>% 
    as.data.frame() %>%
    melt() %>% 
    ggplot(aes(x = value, color = variable)) +
    geom_density() +
    geom_vline(xintercept = get(paste0("two_",x))[1], linetype = "dashed") +
    geom_vline(xintercept = get(paste0("two_",x))[2], linetype = "dashed") +
    labs(title = x)
  
  
  
}) %>% set_names(parameters)

