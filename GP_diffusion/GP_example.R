source("setup.R")


# Following is an example of non parametric regression with GPs using Stan

# Generate data

x <- runif(n = 50, min = -1, max = 1)
y <- sin(3*x) + rnorm(length(x), mean = 0, sd = .25)

data <- data.frame(x, y)

# Check how it looks
plot(data)

# Compile Stan program
gp_model <- stan_model("GP_diffusion/GP.stan")

# Grid for posterior predictions
x_pred <- seq(from = -1.5, to = 1.5, by = 0.01)

# MCMC, sample the posterior with stan. See "GP.stan" for details on input data parameters.
gp_samples <- sampling(gp_model, 
                       data = list(N = nrow(data), 
                                   N_pred = length(x_pred), 
                                   y = data$y, 
                                   x = data$x, 
                                   x_pred = x_pred, 
                                   length_scale = 1, # lambda
                                   stat_var = 1,     # alpha^2
                                   error = .5),      
                       iter = 2000, chains = 2       # MCMC iterations and chains
                       )


# Get posterior of latent "true" signal f
f_pred <- rstan::extract(gp_samples, "f_pred")[[1]] 

# Transpose and add x_pred
f_pred <- f_pred %>% 
  t %>%
  as.data.frame() %>% 
  mutate(x = x_pred)


# Get posterior summary of f
summary <- summary(gp_samples)$summary
f_summary <- summary[grep("f_pred", rownames(summary)), ] %>% 
  as.data.frame() %>% 
  mutate(x = x_pred)

# Edit column names 
f_summary <- f_summary[, c("mean", "2.5%", "97.5%", "x")] %>% 
  set_colnames(c("mean", "lower_2.5", "upper_2.5", "x"))

# Plot

ggplot() + 
  # posterior mean and 95% interval
  geom_line(data = f_summary, 
            aes(x = x, y = mean), color = "blue", size = 1) +
  geom_ribbon(data = f_summary, 
              aes(x = x, ymin = lower_2.5, ymax = upper_2.5), 
              fill = "blue", alpha = 0.5) +
  # Posterior samples, take 25
  geom_line(data = f_pred[, sample(1:(ncol(f_pred) - 1), 25)] %>% 
              cbind(x = f_pred$x) %>% 
              melt("x"), 
            aes(x = x, y = value, group = variable), size = 0.1) +
  # Original signal
  geom_line(data = data.frame(x = x_pred, 
                              y = sin(3*x_pred)), 
            aes(x = x, y = y), color = "red", size = 1) + 
  # Data
  geom_point(data = data, aes(x = x, y = y))
  
