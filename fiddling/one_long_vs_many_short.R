# one long vs many short

# Stan models

oup_fitter <- stan_model("fiddling/stan_models/fit_oup.stan")
hierarchical_pooled_gp_oup <- stan_model("fiddling/stan_models/fit_hierarchical_oup.stan")

# Simulator
oup_simulator <- stan_model("fiddling/stan_models/simulate_gp_OUP.stan")

## Data ****************************** ####

rho_true <- 8
alpha_true <- 8
sigma <- 1
n_observations <- 100
n_total_samples <- 101
x_total <- 100 * (0:(n_total_samples - 1)) / (n_total_samples - 1)
seed <- 11235

# Divide long series into this many short series
into <- 

# Simulate data (including hold out data)

gp_simulator_data <- list(N = n_total_samples,
                          x = x_total,
                          alpha = alpha_true,
                          rho = rho_true,
                          sigma =  sigma)

long_samples <- sampling(oup_simulator,
                    data = gp_simulator_data, 
                    iter=1,
                    chains=1,
                    seed=seed,
                    algorithm="Fixed_param")


fs <- grep("f", rownames(summary(long_samples)$summary), value = T)
ys <- grep("y", rownames(summary(long_samples)$summary), value = T)


df <- data.frame(f = summary(long_samples)$summary[fs, "mean"],
                 y = summary(long_samples)$summary[ys, "mean"],
                 rho = rho_true,
                 alpha = alpha_true,
                 x = x_total)

long_data <- list(y = df[x_total %% 1 == 0, "y"] %>%
                    head(n = n_observations),
                  N = n_observations,
                  x = x_total[x_total %% 1 == 0] %>%
                    head(n = n_observations))


short_data <- divide_long_series(long_data, into = into)

## Stan
# 
# long_fit <- sampling(oup_fitter,
#                      long_data,
#                      iter = 1000,
#                      chains = 1,
#                      init = 1)

short_fit <- sampling(hierarchical_pooled_gp_oup,
                      short_data,
                      iter = 1000,
                      chains = 1,
                      init = 1)

## Results


long_MAP <- summary(long_fit)$summary[c("rho", "alpha"), "50%"]  %>% t%>% as.data.frame()

short_MAP_rho <- summary(short_fit)$summary[grepl("rho\\[", rownames(summary(short_fit)$summary)), "50%"] %>% 
  as.data.frame() %>% 
  set_colnames("rho")

short_MAP_alpha <- summary(short_fit)$summary[grepl("alpha\\[", rownames(summary(short_fit)$summary)), "50%"] %>% 
  as.data.frame() %>% 
  set_colnames("alpha")


short_MAP <- cbind(short_MAP_rho, short_MAP_alpha, index = 1:10)




## Plot
# long_prior <- prior_data(y_shape = "normal", y1=0, y2=1
#                     , x_shape = "invgamma", 4, 10,
#                     x = 1:100/10, y = 1:100/10)

short_prior <- prior_data(y_shape = "normal", y1=summary(short_fit)$summary["alpha_mean", "50%"], y2=summary(short_fit)$summary["alpha_sd", "50%"],
                          x_shape = "invgamma", x1 = summary(short_fit)$summary["rho_shape", "50%"], x2 = summary(short_fit)$summary["rho_rate", "50%"],
                          x = 1:100/10, y = 1:100/10)


ggplot() +
  geom_abline(intercept = 0, slope = alpha_true/rho_true, linetype = "dashed") +
  geom_contour(data = long_prior, aes(x = x, y = y, z = value), color = "blue") +
  geom_contour(data = short_prior, aes(x = x, y = y, z = value), color = "red") +
  geom_point(data = long_MAP, aes(x = rho, y = alpha), color = "blue") +
  geom_point(data = short_MAP, aes(x = rho, y = alpha), color = "red") +
  geom_point(data = data.frame(rho = rho_true, alpha = alpha_true), aes(x = rho, y = alpha), size = 3) +
  coord_fixed(xlim = c(0, 10), ylim = c(0, 10)) +
  labs(x = "Rho", y = "Alpha", subtitle = paste0("Contours: fixed/learned prior (blue/red) \n 1*100 samples  vs.  ", into, "*", 100/into, "samples"))
  
  








