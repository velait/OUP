# Mixture of experts: two OUPs

# x values
n_obs <- 101
x_val <- seq(from = 0, to = 100, length.out = n_obs)

## Functions

# sigmoid function as the latent process
sigmoid_function <- function(x) {
  
  return(1/(1 + exp(-.5*x + 25)))
}

# Step function
stepper <- function(x) {
  
  ifelse(x > 33 & x < 66, 0, 1)
  
}

normie <- function(x) {
  1 - (dnorm(x, 50, 5)/dnorm(50, 50, 5))
}


line <- function(x) {
  
  x/100
  
}

func <- sigmoid_function

# oup parameters
mu <- c(-1, 1)
sigma <- c(.5, 1.2)
lambda <- c(.5, 2)

parameters <- c("mu", "sigma", "lambda")

## generate an observation given the previous one

# generator
oup_mixture <- function(y, delta_t, weight, mu, sigma, lambda, seed) {
  
  # set.seed(seed)
  
  mean <- mu - (mu - y)*exp(-lambda*delta_t)
  sd <- (sigma^2/(2*lambda))*(1 - exp(-2*lambda*delta_t))
  
  oup <- c(rnorm(1, mean = mean[1], sd = sd[1]), 
           rnorm(1, mean = mean[2], sd = sd[2]))
  
  choose <- sample(1:2, 1, prob = c(1-weight, weight))
  
  return(oup[choose])
}


# bimodal series
seed <- sample(1:1000, 1)
set.seed(seed)
y <- 0
for(x in 2:length(x_val)) {
  print(x)
  delta_t <- x_val[x] - x_val[x-1]
  
  obs <- oup_mixture(y[x-1], delta_t, func(x_val[x]), mu, sigma, lambda, seed = seed)
  
  y <- c(y, obs)
}

df <- data.frame(x = x_val, y = y, latent = func(x_val))


# Sliding window variance
# variation <- data.frame(x = x_val[1:(.9*length(x_val))], variance = NA)
# for(x in 1:(.9*length(x_val))) {
#   
#   x1 <- x
#   x2 <- x + 2
#   
#   variation[x, "variance"] <- var(y[x1:x2])
#   
#   
# }


# Plot
mixture_plot <- ggplot() +
  theme_bw() +
  geom_line(data = df, aes(x = x, y = latent), color = "red") +
  geom_line(data = df, aes(x = x, y = y)) +
  labs(subtitle = "Red = mixture weight, Black = mixture of experts, two OUPs")
  
mixture_plot
  


## Stan
fit_oup_mixture <- stan_model("mixture_of_experts/fit_oup_mixture.stan")

oup_mixture_samples <- sampling(fit_oup_mixture,
                                list(N = n_obs, y = y, x = x_val),
                                chains = 1,
                                iter = 1000, 
                                init = -10)


# Plot the inferred latent function

fit_weights <- summary(oup_mixture_samples)$summary[grep("weight", rownames(summary(oup_mixture_samples)$summary)), c("mean","25%", "50%", "75%")]

# Plot mixture weight
mixture_weights <- cbind(x = x_val, fit_weights) %>% 
  set_colnames(c("x", "mean", "lower", "mode", "upper")) %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_line(data =  data.frame(x = x_val, latent = func(x_val)), aes(x = x, y = latent), color = "red") +
  geom_line(aes(x = x, y = mode), color = "blue") +
  geom_ribbon(aes(x = x, ymin = lower, ymax = upper), fill = "grey", alpha = .25) +
  geom_line(data = df, aes(x = x, y = y)) +
  theme_bw() +
  labs(subtitle = "Red = mixture weight; Blue = MAP weight")

mixture_weights

# Estimates
mixture_estimates <- summary(oup_mixture_samples)$summary[grep(paste(parameters, collapse = "|"), rownames(summary(oup_mixture_samples)$summary)), 
                                     c("25%", "50%", "75%", "2.5%", "97.5%")] %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  set_colnames(c("parameter", "lower25", "mode", "upper75", "lower2.5", "upper97.5")) %>% 
  cbind(true_value = c(mu, sigma, lambda)) %>% 
  ggplot() +
  geom_point(aes(x = parameter, y = mode), shape = 3, size = 2) + 
  geom_errorbar(aes(x = parameter, ymin = lower25, ymax = upper75), width = .05) +
  geom_point(aes(x = parameter, y = true_value), color = "red", alpha = 0.75) +
  labs(y = "Estimate", subtitle = "Red = true values, black = 50% credible interval")

mixture_estimates

# Save plots

# png("mixture_of_experts/figures/weights_normal.png", height = 500, width = 1000)
# print(mixture_weights)
# dev.off()
# 
# 
# png("mixture_of_experts/figures/estimates_normal.png", height = 750, width = 750)
# print(mixture_estimates)
# dev.off()
# 

