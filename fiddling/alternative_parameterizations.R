# test alternative parameterizations for the hierarchical model

seed <- 11235

# n_series, n_observations
n_series <- 50
n_observations <- 10


## Data

# hyperparameters
rho_shape <- 25  # length scale
rho_rate <- 65
alpha_mean <- 2    # marginal variance
alpha_sd <- 1
sigma <- 0.25 # measurement error

# parameters
parameters <- c("rho", "alpha", "sigma")
rho_true <- rinvgamma(n_series, shape = rho_shape, rate = rho_rate)       # marginal variance (~ kappa^.5)
alpha_true <- restricted_rnorm(n_series, alpha_mean, alpha_sd, lower = 0)        # length scale (~ inv_lambda^.5)
sigma_true <- rep(sigma, n_series)                              # measurement error


hierarchical_test_data <- genera_gp_set_stan(n_series = n_series,
                                  alpha = alpha_true,
                                  rho = rho_true, 
                                  sigma = sigma_true,
                                  intervals = 1:n_observations,
                                  stan_model = oup_simulator,
                                  seed = seed)
    
true_ratio <- rho_true/alpha_true

## basic case ******************************** ####
hierarchical_pooled_gp_oup <- stan_model("fiddling/stan_models/fit_hierarchical_oup.stan")
samples0 <- sampling(hierarchical_pooled_gp_oup,
                     hierarchical_test_data,
                     iter=iter,
                     chains=chains,
                     init=1)



df0 <- lapply(c("rho", "alpha"), function(par) {
  print(par)
  cbind(summary(samples0)$summary[grep(paste0(par, "\\["), summary(samples0)$summary %>% rownames), c("25%", "75%", "mean")], 
        par = get(paste0(par, "_true"))) %>% 
    as.data.frame() %>% 
    arrange(par) %>% 
    mutate(index = 1:n_series) %>% 
    set_colnames(c("lower", "upper", "mean", par, "index"))
}) %>% set_names(c("rho", "alpha"))

ggplot(data = df0[["rho"]]) +
  geom_line(aes(x = index, y = rho)) +
  geom_errorbar(aes(x = index, ymin = lower, ymax = upper)) +
  labs(title = "rho")
  

sim_ratio0 <- df0[["rho"]]$mean/df0[["alpha"]]$mean


ggplot(data.frame(true = true_ratio, sim0 = sim_ratio0), aes(x = true, y = sim0)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")

## rho --> inv_rho *************************** ####
inv_rho_model <- stan_model("fiddling/stan_models/alt_models/fit_hierarchical_oup2.stan")

samples1 <- sampling(inv_rho_model,
         hierarchical_test_data,
         iter=iter,
         chains=chains,
         init=1)


df1 <- lapply(c("rho", "alpha"), function(par) {
  print(par)
  cbind(summary(samples1)$summary[grep(paste0("^",par, "\\["), summary(samples1)$summary %>% rownames), c("25%", "75%", "mean")], 
        par = get(paste0(par, "_true"))) %>% 
    as.data.frame() %>% 
    arrange(par) %>% 
    mutate(index = 1:n_series) %>% 
    set_colnames(c("lower", "upper", "mean", par, "index"))
}) %>% set_names(c("rho", "alpha"))

ggplot(data = df1[["rho"]]) +
  geom_line(aes(x = index, y = rho)) +
  geom_errorbar(aes(x = index, ymin = lower, ymax = upper)) +
  labs(title = "rho")


sim_ratio1 <- df1[["rho"]]$mean/df1[["alpha"]]$mean


ggplot(data.frame(true = true_ratio, sim1 = sim_ratio1), aes(x = true, y = sim1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")



df01 <- data.frame(true = true_ratio,
                   sim1 = sim_ratio1,
                   sim0 = sim_ratio0,
                   index = 1:n_series)

df01 <- rbind(
  
  cbind(df01[, c(1, 3, 4)] %>%
              mutate(type = "t0")) %>% 
        set_colnames(c("true", "sim", "index", "version")),
      
      cbind(df01[, c(1, 2, 4)] %>%
              mutate(type = "t1") %>%
              set_colnames(c("true", "sim", "index", "version")))
     )

ggplot() +
  geom_line(data = df01,aes(x = true, y = sim, group = index)) +
  geom_point(data = df01, aes(x = true, y = sim, color = version)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")

