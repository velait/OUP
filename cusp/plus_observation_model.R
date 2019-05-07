# Cusp + observation model


## Cusp
N <- 10
puolet <- 5

r <- rep(1, N)
alpha <- c(0, 0, 0, 0, 0, -0.5, -0.5, -0.5, -1, -1)
beta <- c(2, 2, 1.5, 1.5, -1, -1, 1, 1, 0, 0)
lambda <- rnorm(N, 0, .25)
epsilon <- rnorm(N, 1, 0.1)

grid <- seq(from = 0, to = 100, by = 0.1)
seed = sample(1:1000, 1)


inits <- rnorm(N, lambda[1], 0.1)

df_shoji <- lapply(1:N, function(i) shoji_generator(y0 = inits[i], times = grid,
                                                    r = r, alpha = alpha[i], beta = beta[i],
                                                    lambda = lambda[i], epsilon = epsilon[i],
                                                    seed = sample(1:1000, 1))) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  cbind(time = grid) %>%
  set_colnames(c(paste0("y", 1:N), "x"))


df_shoji <- thin(df_shoji)

p <- df_shoji %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value)) +
  geom_line() +
  facet_wrap(~variable) +
  labs(title = "Latent cusp process", subtitle = "100 observations, dt = 1")


p




## Softmax for probabilities 
## Softmax works for negative values as well
theta_df <- lapply(1:nrow(df_shoji), function(i) {
  
  df <- df_shoji[i, ]
  
  x <- df[1:(length(df) - 1)] %>% 
    softMax
  
  cbind(x, df[length(df)])
  
}) %>% do.call(rbind, .)

p <- theta_df %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value)) +
  geom_line() +
  facet_wrap(~variable) +
  labs(title = "Cusp + softmax", subtitle = "100 observations, dt = 1")


p

softmax_id <- function(x) {

  x <- c(x, 0)
  
  softMax(x);
}


## Multinomial sampling


obs_df <- apply(theta_df[, -ncol(theta_df)], 1, function(i) {
  rmultinom(1, 1e5, i)
}) %>% t %>%
  as.data.frame() %>%
  cbind(x = theta_df$x) %>%
  set_colnames(c(paste0("y", 1:N), "x"))


p <- obs_df %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value)) +
  geom_line() +
  facet_wrap(~variable) +
  labs(title = "Cusp + softmax + multinomial", subtitle = "100 observations, dt = 1")


p


## Inference with Stan ******************************** ####

## Uninformative priors ******************************
stan_data <- list(N_obs = nrow(obs_df),
                  N_OTUs = N,
                  x = obs_df$x,
                  y = t(select(obs_df, -x)))


cusp_plus_obsevations <- stan_model("cusp/cusp_plus_observation.stan")

cusp_plus_obsevations_samples <- sampling(cusp_plus_obsevations, stan_data,
                                          iter = 1000, 
                                          chains = 1, 
                                          control = list(adapt_delta = 0.8))


## Informative priors *********************************

# MLE estimates for stationary densities

for(i in 1:N) {
  
  x <- 
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
  
  
}







prior_parameters <- matrix(NA, 4, 2)


## Results ******************************************** ####

# get results
results_df <- lapply(1:4, function(x) get_stan_results(cusp_plus_obsevations_samples, paste0("theta\\[.+,", x,"\\]"))) %>% do.call(rbind, .) %>% 
  cbind(parameter = rep(c("alpha", "beta", "lambda", "epsilon"), each = N),
        index = factor(rep(1:N, 4)),
        true = c(alpha, beta, lambda, epsilon))


# plot estimates and true values
results_df %>% 
  ggplot() + 
  geom_point(aes(x = index, y = true), shape = 4) + 
  geom_errorbar(aes(x = index, ymax = upper97.5, ymin = lower2.5), width = .2) +
  facet_wrap(~parameter)


# compare true, estimated and data densities

x <- seq(-5, 5, length.out = 100)
true_d <- lapply(1:N, function(i) {
  
  d <- cc_density(x, r[i], alpha[i], beta[i], lambda[i], epsilon[i])
  
  data.frame(x, d, r[i], alpha[i], beta[i], lambda[i], epsilon[i], "true", i) %>% 
    set_colnames(c("x", "y", "r", "alpha", "beta", "lambda", "epsilon", "type", "index"))
  
}) %>% do.call(rbind, .)


map_d <- lapply(1:N, function(i) {
  
  rr <- 1
  aa <- results_df %>%
    filter(parameter == "alpha", index == i) %>% 
    pull(mode)
  
  bb <- results_df %>%
    filter(parameter == "beta", index == i) %>% 
    pull(mode)
  
  ee <- results_df %>%
    filter(parameter == "epsilon", index == i) %>% 
    pull(mode)
  
  ll <- results_df %>%
    filter(parameter == "lambda", index == i) %>% 
    pull(mode)
  
  d <- cc_density(x, rr, aa, bb, ll, ee)
  
  data.frame(x, d, r[i], alpha[i], beta[i], lambda[i], epsilon[i], "true", i) %>% 
    set_colnames(c("x", "y", "r", "alpha", "beta", "lambda", "epsilon", "type", "index"))
  
}) %>% do.call(rbind, .)


p <- ggplot() +
  geom_line(data = true_d, aes(x = x, y = y)) +
  geom_line(data = map_d, aes(x = x, y = y), color = "red") +
  stat_density(data =df_shoji %>%
                 melt(id.vars = "x") %>%
                 mutate(index = rep(1:N, each=nrow(df_shoji))),
               aes(x = value), geom = "line", linetype = "dashed") +
  facet_wrap(~index) +
  labs(title = "Density estimates for 10 series", subtitle = "Solid = true; dashed = data; red = MAP estimate")

p

