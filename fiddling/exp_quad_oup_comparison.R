## stan models ************************************************************* ####

# generators
oup_simulator <- stan_model("fiddling/stan_models/simulate_gp_OUP.stan")
gp_simulator <- stan_model("fiddling/simulate_gp_exp_quad.stan")

# simulate underlying process for the quantiles
stan_simulator <- stan_model(file='fiddling/simu_gauss_dgp.stan')

# fitters
exp_quad_fitter <- stan_model("fiddling/fit_gp.stan")
oup_fitter <- stan_model("fiddling/fit_oup.stan")
fitter_list <- list(oup = oup_fitter, exp_quad = exp_quad_fitter)


## parameters ************************************************************** ####


# alpha = kappa^2
# rho = (inv_lambda)^2

N_total <- 501
x_total <- 100 * (0:(N_total - 1)) / (N_total - 1)

# alpha_true <- sqrt(kappa)        # marginal variance (~ kappa)
# rho_true <- 1/sqrt(lambda)        # length scale (~ inv_lambda)
# sigma_true <- 0.25   # measurement error

alpha_true <- 4 # marginal variance (~ kappa^0.5)
rho_true <- 2  # length scale (~ inv_lambda^0.5)
sigma_true <- 1    # measurement error

kappa <- alpha_true^2
lambda <- 1/rho_true^2
paste("lambda:",  round(1/rho_true^2, 3))
paste("sigma:",  round(sqrt(2*kappa*lambda), 3))


gp_simulator_data <- list(N = N_total, x = x_total, alpha = alpha_true, rho = rho_true, sigma =  sigma_true)

## Data ******************************************************************** ####

# Generate exp_quad and oup data
seed <- runif(1, 1, 1000)
simulate_oup_and_gp <- lapply(list(oup_simulator, gp_simulator), function(model) {
  samples <- sampling(model,
                         data = gp_simulator_data, 
                         iter=1,
                         chains=1,
                         seed=seed,
                         algorithm="Fixed_param")
  
  samples

  
}) %>% 
  set_names(c("oup", "exp_quad"))


# collect simulated data
simulated_data <- lapply(names(simulate_oup_and_gp), function(i) {
  
  fs <- grep("f", rownames(summary(simulate_oup_and_gp[[i]])$summary), value = T)
  ys <- grep("y", rownames(summary(simulate_oup_and_gp[[i]])$summary), value = T)
  
  
  df <- data.frame(f = summary(simulate_oup_and_gp[[i]])$summary[fs, "mean"],
                   y = summary(simulate_oup_and_gp[[i]])$summary[ys, "mean"],
                   model = i,
                   x = x_total)
  
}) %>% do.call(rbind, .)
  
 
# 
# ggplot() +
#   geom_line(data = simulated_data, aes(x = x, y = f), color ="grey", size = 1) +
#   geom_point(data = simulated_data, aes(x = x, y = y), color = "sienna2", size = 0.5) +
#   geom_point(data = simulated_data %>% filter((x %% 1) == 0), aes(x = x, y = y), color = "black") +
#   facet_wrap(~model, nrow = 2) +
#   theme_bw()
  

# Get the underlying data generating process
data_generating_process <- lapply(c("oup", "exp_quad") , function(x) {
  
  
  
  f_data <- list(sigma=sigma_true,
                 N=N_total,
                 f=simulated_data[simulated_data$model == x, ] %>% pull(f))
  
  
  dgp_fit <- sampling(stan_simulator, data=f_data, iter=1000, warmup=0,
                  chains=1, seed=5838298, refresh=1000, algorithm="Fixed_param")
  
  df <- summary(dgp_fit)$summary[, c("mean", "2.5%", "25%", "50%", "75%", "97.5%")] %>% 
    as.data.frame()
  
  df <- df[1:(nrow(df) - 1), ]
  
  df %>% mutate(model = x, x = x_total)
  return(df)
}) %>% 
  set_names(c("oup", "exp_quad"))

# Plotter for the underlying process
plot_posterior_error <- function(df, x_values = x_total, obs) {
  

  df <- df %>%
    set_colnames(c("mean", "low2.5", "low25", "median", "high75", "high97.5"))
  
  p <- ggplot() +
    geom_ribbon(data = df, aes(x = x_values, ymin = low2.5, ymax = high97.5), fill = "darksalmon") +
    geom_ribbon(data = df, aes(x = x_values, ymin = low25, ymax = high75), fill = "brown3") +
    geom_point(data = obs, aes(x =x, y = y))
    
  p
}


#*** ugly stuff ******
obs <- simulated_data %>% 
  filter((x %% 1) == 0) %>%
  filter(model == "oup") %>%
  select(x, y)

underlying_process_observations_plot <- list()

underlying_process_observations_plot[["oup"]] <- plot_posterior_error(data_generating_process[["oup"]], x_total, obs) + ggtitle("OUP")
  
obs <- simulated_data %>% 
  filter((x %% 1) == 0) %>%
  filter(model == "exp_quad") %>%
  select(x, y)

underlying_process_observations_plot[["exp_quad"]] <- plot_posterior_error(data_generating_process[["exp_quad"]], x_total, obs)+ ggtitle("Squared exponential")



underlying_process_observations_plot <- plot_grid(underlying_process_observations_plot[["oup"]],
                                                  underlying_process_observations_plot[["exp_quad"]],
                                                  ncol = 1)
underlying_process_observations_plot

#************************************

## Fit data **************************************************************** ####

# Fit
gp_fits <- lapply(c("oup", "exp_quad"), function(x) {
  print(x)
  observations <- simulated_data[simulated_data$model == x, ] %>% 
    filter((x %% 1) ==0) %>% pull(y)
  
  samples <- sampling(fitter_list[[x]],
           list(N = 100 + 1, x = 1:(100 + 1), y = observations),
           iter = 1000, 
           chains =1,
           init =1)
  

  
}) %>% set_names(c("oup", "exp_quad"))

gp_results <- lapply(c("oup", "exp_quad"), function(x) {
  samples <- gp_fits[[x]]
  summary(samples)$summary[c("rho", "alpha", "sigma"), c("2.5%", "50%", "97.5%")] %>% 
    as.data.frame() %>% 
    cbind(model = rep(x, 3), parameter = c("rho", "alpha", "sigma")) %>% 
    set_colnames(c("lower", "mode", "upper", "model", "parameter"))
  
}) %>% do.call(rbind, .)

# Compare results plot 
ggplot() +
    geom_point(data = gp_results, aes(x = parameter, y = mode, color = model)) +
    geom_errorbar(data = gp_results, aes(x = parameter, ymin = lower, ymax = upper, color = model), width =.1) +
    geom_point(data = data.frame(parameter = c("alpha", "rho", "sigma"), value = c(alpha_true, rho_true, sigma_true)), aes(x = parameter, y = value))
  
  
  
  
  