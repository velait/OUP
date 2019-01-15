
# ************************************************************************
# ******************** Single series example *****************************
# ************************************************************************
 

library(rstan)
library(magrittr)
library(shinystan)
library(mvtnorm)
library(cowplot)

source("OU.functions.R")

# Compile model
oup_model <- stan_model(file = "stan_models/alt_single_series_oup.stan")




# Set parameter values
mu <- 0
lambda <- 0.5
sigma <- 0.25
intervals <- 1:100

# Simulate data
stan_data <- generate_student_set(n_series = 1,
                                  student_df = 5,
                                  mu = mu,
                                  lambda = lambda,
                                  sigma = sigma,
                                  intervals = intervals,
                                  seed = 11235)


# Plot data
stan_data$Y %>% plot(x = intervals, type = "l")





# Sample from posterior
samples <- sampling(oup_model,
                    stan_data,
                    chains = 1,
                    iter = 1000)


# Get posterior samples
posterior <- sapply(c("lambda", "mu", "sigma"),
                    function(x) rstan::extract(samples, x)) %>%
  set_names(c("lambda", "mu", "sigma"))


# Plot posterior with simulation values
posterior_plot <- lapply(c("lambda", "mu", "sigma"), function(x) {
  
  p <- posterior[[x]] %>% 
    as.data.frame() %>%
    ggplot(aes(x = .)) +
    geom_density() +
    geom_vline(xintercept = eval(parse(text = x)), linetype = "dashed") +
    labs(title = x) +
    scale_x_continuous(limits = c(eval(parse(text = x))-1, eval(parse(text = x))+1))
  
  
}) 


plot_grid(posterior_plot[[1]], posterior_plot[[2]], posterior_plot[[3]], nrow = 1)



# ShinyStan
samples %>% launch_shinystan()




