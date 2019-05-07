# Fit cusp catastrophy distribution with MLE
library(stats4) 

### 1 series **************************************** #####
## Data **************************************** #####

r <- 1
alpha <- -0.3
beta <- 2
lambda <- 0
epsilon <- 0.1

grid <- seq(from = 0, to = 100, by = 0.1)
# seed = sample(1:1000, 1)
# set.seed(1)

inits <- c(-2)

df_shoji <- lapply(inits, function(y0) shoji_generator(y0 = y0, times = grid,
                                                       r = r, alpha = alpha, beta = beta, lambda = lambda, epsilon = epsilon,
                                                       seed = 2)) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  cbind(time = grid) %>%
  # set_colnames(c("y", "x"))
  set_colnames(c("shoji", "x"))

p <- df_shoji %>%
  ggplot(aes(x = x, y = shoji)) +
  geom_line() +
  labs(x ="time", y = "value", subtitle = "r = 1, epsilon = 1")


p



## Compare with euler discretization. 
# df_euler <- lapply(inits, function(y0) em_generator(y0 = y0, grid)) %>%
#   do.call(cbind, .) %>%
#   as.data.frame() %>%
#   cbind(time = grid) %>%
#   # set_colnames(c("y", "x"))
#   set_colnames(c("euler", "x"))


# df <- full_join(df_shoji, df_euler, "x") %>% 
#   melt(id.vars = "x")


## MLE for cross-sectional ********************* #####   
x <- df_shoji$shoji
cusp_mLL <- function(par) {
  
  # alpha <- par[1]
  # beta <- par[2]
  # lambda <- par[3]
  # epsilon <- par[4]
  
  # Normalizing factor
  M <- integrate(
    function(x, alpha = par[1], beta = par[2], lambda = par[3], epsilon = par[4]) {
      exp((alpha*(x - lambda) + .5*beta*(x - lambda)^2 - .25*(x - lambda)^4 ) / epsilon)
    },
    -Inf, Inf)$value
  
  # Log likelihood
  LL <- (par[1]*(x - par[3]) + .5*par[2]*(x - par[3])^2 - .25*(x - par[3])^4 ) / par[4] - log(M)
  
  return((-sum(LL)))
}
ml_est <- optim(c(0, 0, 0, 1),
      cusp_mLL, method = "L-BFGS-B",
      upper = c(20, 20, 20, 20),
      lower = c(-20, -20, -20, 0.025),
      control = list(maxit = 1000))


# Compare true and estimated densities

x <- seq(-5, 5, length.out = 100)
compare_density_df <- data.frame(x = x,
                                 y_true = cc_density(x, r, alpha, beta, lambda, epsilon),
                                 y_est = cc_density(x, r, alpha = ml_est$par[1],
                                                    beta = ml_est$par[2], lambda = ml_est$par[3], epsilon = ml_est$par[4])) %>%
  melt(id.vars = "x")

inv_density_plot <- ggplot() +
  stat_density(data = df_shoji, aes(x = shoji), linetype = "dashed", geom="line") +
  geom_line(data = compare_density_df,aes(x = x, y = value, color = variable)) +
  labs(subtitle = "Dashed black = data")

inv_density_plot


## Fit time series with Stan ******************* ####

fit_shoji_cusp <- stan_model(file = "cusp/shoji_cusp.stan")

thin_shoji_df <- df_shoji[(df_shoji$x %% 1) == 0, ]

shoji_data <- list(N = nrow(thin_shoji_df), x = thin_shoji_df$x, y = thin_shoji_df$shoji)

shoji_samples <- sampling(fit_shoji_cusp, shoji_data,
                          iter = 4000, chains = 4, 
                          control = list(max_treedepth = 12))



### Many series ************************************* #####

## Data **************************************** #####

r <- 1
alpha <- -0.25
beta <- 2
lambda <- 0
epsilon <- 1


theta <- c(alpha, beta, lambda, epsilon)

grid <- seq(from = 0, to = 25, by = 0.01)
seed = sample(1:1000, 1)

N <- 50

inits <- rnorm(N, lambda, 1)

df_shoji <- lapply(inits, function(y0) shoji_generator(y0 = y0, times = grid,
                                                       r = r, alpha = alpha, beta = beta, lambda = lambda, epsilon = epsilon,
                                                       seed = seed)) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  cbind(time = grid) %>%
  set_colnames(c(paste0("y", 1:N), "x"))


## MLE for cross-sectional ********************* ####

df_shoji_vec <- df_shoji %>%
  filter((x %% 1) == 0) %>% 
  select(-x) %>%
  as.matrix %>% as.vector()

x <- df_shoji_vec

cusp_mLL <- function(par) {
  
  # alpha <- par[1]
  # beta <- par[2]
  # lambda <- par[3]
  # epsilon <- par[4]
  
  # Normalizing factor
  M <- integrate(
    function(x, alpha = par[1], beta = par[2], lambda = par[3], epsilon = par[4]) {
      exp(1*(alpha*(x - lambda) + .5*beta*(x - lambda)^2 - .25*(x - lambda)^4 ) / epsilon)
    },
    -Inf, Inf)$value
  
  # Log likelihood
  LL <- 1*(par[1]*(x - par[3]) + .5*par[2]*(x - par[3])^2 - .25*(x - par[3])^4 ) / par[4] - log(M)
  
  return(log(-sum(LL)))
}
ml_est <- optim(c(0, 0, 0, 1),
                cusp_mLL, method = "L-BFGS-B",
                upper = c(20, 20, 20, 20),
                lower = c(-20, -20, -20, 0.2),
                control = list(maxit = 1000))

# Bootstrap for confidence intervals

# boot_mle_est <- sapply(1:100, FUN = function(i) {
#   
#   x <- sample(df_shoji_vec, size = length(df_shoji_vec), replace = TRUE)
#   cusp_mLL <- function(par) {
#     
#     # alpha <- par[1]
#     # beta <- par[2]
#     # lambda <- par[3]
#     # epsilon <- par[4]
#     
#     # Normalizing factor
#     M <- integrate(
#       function(x, alpha = par[1], beta = par[2], lambda = par[3], epsilon = par[4]) {
#         exp((alpha*(x - lambda) + .5*beta*(x - lambda)^2 - .25*(x - lambda)^4 ) / epsilon)
#       },
#       -Inf, Inf)$value
#     
#     # Log likelihood
#     LL <- (par[1]*(x - par[3]) + .5*par[2]*(x - par[3])^2 - .25*(x - par[3])^4 ) / par[4] - log(M)
#     
#     return((-sum(LL)))
#   }
#   est <- optim(c(0, 0, 0, 1),
#           cusp_mLL, method = "L-BFGS-B",
#           upper = c(20, 20, 20, 20),
#           lower = c(-20, -20, -20, 0.025),
#           control = list(maxit = 1000))
#   
#   est$par
#   
# }) %>% t

# Compare true and estimated densities

x_vals <- seq(-5, 5, length.out = 100)
compare_density_df <- data.frame(x = x_vals,
                                 y_true = cc_density(x_vals, r = 1, alpha, beta, lambda, epsilon),
                                 y_est = cc_density(x_vals, r = 1, alpha = ml_est$par[1],
                                                    beta = ml_est$par[2], lambda = ml_est$par[3], epsilon = ml_est$par[4])) %>%
  melt(id.vars = "x")

inv_density_plot <- ggplot() +
  stat_density(data = data.frame(dens = x), aes(x = dens), linetype = "dashed", geom="line") +
  geom_line(data = compare_density_df,aes(x = x, y = value, color = variable)) +
  labs(subtitle = "Dashed black = data")

inv_density_plot

## Fit time series with Stan ******************* ####

fit_shoji_cusp <- stan_model(file = "cusp/shoji_cusp.stan")
uninformative_fit_shoji_cusp <- stan_model(file = "cusp/uninformative_shoji_cusp.stan")

thin_shoji_df <- df_shoji[(df_shoji$x %% 1) == 0, ]


shoji_estimates <- lapply(1:(ncol(thin_shoji_df) - 1), function(i) {
  
  shoji_data <- list(N = nrow(thin_shoji_df), x = thin_shoji_df$x, y = thin_shoji_df[, i])
  
  shoji_samples <- sampling(fit_shoji_cusp, shoji_data,
                            iter = 1000, chains = 2)
  
  
  res <- summary(shoji_samples)$summary %>% 
    grep_rows("theta") %>% 
    as.data.frame() %>% 
    select(c("2.5%", "50%", "97.5%")) %>% 
    mutate(index = i)
    
}) %>% do.call(rbind, .)

uninformative_shoji_estimates <- lapply(1:(ncol(thin_shoji_df) - 1), function(i) {
  
  shoji_data <- list(N = nrow(thin_shoji_df), x = thin_shoji_df$x, y = thin_shoji_df[, i])
  
  shoji_samples <- sampling(uninformative_fit_shoji_cusp, shoji_data,
                            iter = 1000, chains = 2)
  
  
  res <- summary(shoji_samples)$summary %>% 
    grep_rows("theta") %>% 
    as.data.frame() %>% 
    select(c("2.5%", "50%", "97.5%")) %>% 
    mutate(index = i)
  
}) %>% do.call(rbind, .)

shoji_estimates <- shoji_estimates %>%
  mutate(theta = rep(1:4, N)) %>% 
  set_colnames(c("lower", "mode", "upper", "index", "theta")) %>% 
  mutate(test = "informative")

uninformative_shoji_estimates <- uninformative_shoji_estimates %>%
  mutate(theta = rep(1:4, N)) %>% 
  set_colnames(c("lower", "mode", "upper", "index", "theta")) %>% 
  mutate(test = "uninformative")

both_estimates <- rbind(shoji_estimates, uninformative_shoji_estimates) %>% 
  mutate(true = NA)

for(i in 1:nrow(both_estimates)) {
  th <- both_estimates[i, "theta"] 
  both_estimates[i, "true"] <- theta[th]
  
}


theta_names <- c("1" = "alpha", "2" = "beta", "3" = "lambda", "4" = "epsilon")

both_estimates %>% 
  ggplot() +
  geom_errorbar(aes(x = index, ymin = lower, ymax = upper, color = test), position = "dodge") + 
  geom_hline(aes(yintercept = true), linetype = "dashed") + 
  facet_wrap(~theta, labeller = as_labeller(theta_names))


