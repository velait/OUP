## Gaussian process + generic ews comparison. 
## No hierarchical models
## Data from ecological model
## Sliding windows


## Data ***************************************************** ####
set.seed(242425)
end <- 1000
times <- seq(from = 0, to = end, by = 0.1)
N_series <- 1
n_obs <- 900




# Simulation parameters
r <- 1
K <- 10
cs <- 1
h <- 1
sigma <- 0.01
obs_error_sd <- .1

## Simulate data with CSD
ews_set <- lapply(1:N_series, function(j) {
  
  y <- rep(NA, length(times))
  
  # Initial value
  y[1] <- 8
  
  cs <- 1 + (2.6771 - 1)*(times/max(times))
  
  
  for(i in 2:length(times)) {
    
    dt <- times[(i-1):i]
    
    finer_grid <- seq(from = dt[1], to = dt[2], by = 0.01)
    
    y[i] <- ews_generator(y[i-1], finer_grid, c = cs[i], milstein = T)[length(finer_grid)]
    
    
  }
  
  
  # Add measurement error
  y <- y + rnorm(length(y), 0, obs_error_sd)
  
  y
  
}) %>%
  do.call(cbind, .) %>%
  cbind(x = times) %>%
  as.data.frame()

ews_set <- ews_set[1:(which(((1 + (2.6771 - 1)*(times/max(times))) < 2.604)) %>% max()), ]




## Residual data
ews_residuals <- lapply(1:N_series, function(i) {
  
  series <- ews_set[, i]
  smoothed <- smth(series, window = .1, tails = TRUE)
  res <- series - smoothed
  return(res)
}) %>% do.call(cbind, .) %>% cbind(., ews_set$x) %>% 
  as.data.frame %>% 
  set_colnames(c(paste0("V", 1:N_series), "x"))


obs_index <- sample((length(ews_residuals$x))/10 %>% floor, n_obs)*10 + sample(-4:5, n_obs, replace = T)


ews_residuals_obs <- ews_residuals[obs_index, ] %>%
# ews_residuals_obs <- ews_residuals %>%
  as.data.frame() %>%
  set_colnames(c(paste0("V", 1:N_series), "x")) %>% 
  arrange(x)
  


## Stan: Separate ******************************************* ####
oup_model <- stan_model("stan_models/oup_fitter.stan")
se_model <- stan_model("stan_models/square_exp_fitter.stan")
matern_1.5_model <- stan_model("stan_models/matern_1.5_fitter.stan")
matern_2.5_model <- stan_model("stan_models/matern_2.5_fitter.stan")
ar1_model <- stan_model("stan_models/ar1.stan")
ar1_missing_model <- stan_model("stan_models/ar1_missing.stan")


iter <- 500
chains <- 1



# CDS set
# ar1_ews_res <- ar1_ews_set(ews_residuals_obs, iter, chains, window_prop = window_prop)

# generic_ews_res <- generic_ews_set(ews_residuals_obs, window_prop = .05)


ar1_ews_res <- ar1_missing_ews_set(ews_residuals_obs, iter = 500, chains, window_length = 50)

oup_ews_res <- GP_ews_set(ews_residuals_obs, iter = 500, chains = 1, window_length = 50, kernel = "OUP")
se_ews_res <- GP_ews_set(ews_residuals_obs, iter = 500, chains = 1, window_length = 50, kernel = "se")
matern_1.5_ews_res <- GP_ews_set(ews_residuals_obs, iter = 500, chains = 1, window_length = 50, kernel = "matern_1.5")
matern_2.5_ews_res <- GP_ews_set(ews_residuals_obs, iter = 500, chains = 1, window_length = 50, kernel = "matern_2.5")





# ews_res <- list("ar1" = ar1_ews_res, "oup" = oup_ews_res, "se" = se_ews_res, "generic" = generic_ews_res)
ews_res <- list("ar1" = ar1_ews_res,
                "oup" = oup_ews_res,
                "se" = se_ews_res, 
                "matern_1.5" = matern_1.5_ews_res, 
                "matern_2.5" = matern_2.5_ews_res)


# saveRDS(object = ews_res, file = "EWS/res/GP_EWS_one_long_1000_NA.RDS")
# ews_res <- readRDS(file = "EWS/res/GP_EWS_one_long_1000_NA.RDS")
# ews_res <- readRDS(file = "EWS/res/.RDS")

## Results ************************************************** ####

## Check parameter values are OK
par_res <- lapply(1:N_series, function(s) {
  
  my_res <- ews_res
  
  # # # Generic parameters
  # generic_df <- my_res[["generic"]] %>%
  #   filter(series == s) %>%
  #   dplyr::select(ar1, sd, sk, kurt, cv, densratio, timeindex) %>%
  #   set_colnames(c(paste0("generic_", c("ar1", "sd", "sk", "kurt", "cv", "densratio")), "timeindex"))%>% 
  #   melt(id.vars = "timeindex")
  
  # AR(1)
  ar1_df <- my_res[["ar1"]] %>% 
    filter(series == s) %>% 
    dplyr::select(lambda_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("ar1_", c("rho", "sigma")), "timeindex"))%>% 
    melt(id.vars = "timeindex")
  
  # OUP 
  oup_df <- my_res[["oup"]] %>% 
    filter(series == s) %>%
    dplyr::select(rho_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("OUP_", c("rho", "sigma")), "timeindex")) %>% 
    melt(id.vars = "timeindex")
  
  # SE 
  se_df <- my_res[["se"]] %>%
    filter(series == s) %>%
    dplyr::select(rho_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("se_", c("rho", "sigma")), "timeindex")) %>% 
    melt(id.vars = "timeindex")
  
  # Matern 1.5
  matern_1.5_df <- my_res[["matern_1.5"]] %>%
    filter(series == s) %>%
    dplyr::select(rho_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("matern_1.5_", c("rho", "sigma")), "timeindex")) %>% 
    melt(id.vars = "timeindex")
  
  # Matern 2.5
  matern_2.5_df <- my_res[["matern_2.5"]] %>%
    filter(series == s) %>%
    dplyr::select(rho_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("matern_2.5_", c("rho", "sigma")), "timeindex")) %>% 
    melt(id.vars = "timeindex")
  
  return(rbind(se_df, ar1_df, oup_df, matern_1.5_df, matern_2.5_df) %>% mutate(series = s))   
  
  
}) %>% do.call(rbind, .) %>% as.data.frame()

par_res %>%
  ggplot(aes(x = timeindex, y = value, color = as.factor(series))) + 
  geom_smooth() +
  facet_wrap(~variable, scales = "free")



par_res %>%
  ggplot(aes(x = timeindex, y = value, color = variable)) + 
  geom_line(aes(group = series)) + 
  facet_wrap(~variable, scales = "free")



par_res_scaled <- lapply(unique(par_res$variable), function(par) {
  
  lapply(1:N_series, function(i) {
    
    df <- par_res %>%
      filter(variable == par, series == i)
    
    df <- df %>%
      mutate(scaled = value/max(abs(df$value)))
    
  }) %>%
    do.call(rbind, .)
  
  
}) %>% do.call(rbind, .)


par_res_scaled %>%
  ggplot(aes(x = timeindex, y = scaled, color = variable)) + 
  geom_smooth() + 
  facet_wrap(~variable, scales = "free")

# Kendall's Tau of all variables. For the Bayesian posteriors mode is used.
tau_res <- lapply(1:N_series, function(s) {
  
  my_res <- ews_res
  
  # # Generic parameters
  # generic_df <- my_res[["generic"]] %>%
  #   filter(series == s) %>%
  #   dplyr::select(ar1, sd, sk, kurt, cv, densratio, timeindex) %>%
  #   set_colnames(c(paste0("generic_", c("ar1", "sd", "sk", "kurt", "cv", "densratio")), "timeindex"))
  # 
  # AR(1)
  ar1_df <- my_res[["ar1"]] %>% 
    filter(series == s) %>% 
    dplyr::select(lambda_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("ar1_", c("lambda", "sigma")), "timeindex"))%>% 
    melt(id.vars = "timeindex")
  
  # OUP 
  oup_df <- my_res[["oup"]] %>% 
    filter(series == s) %>%
    dplyr::select(rho_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("OUP_", c("rho", "sigma")), "timeindex")) %>% 
    melt(id.vars = "timeindex")
  
  # SE 
  se_df <- my_res[["se"]] %>%
    filter(series == s) %>%
    dplyr::select(rho_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("se_", c("rho", "sigma")), "timeindex")) %>% 
    melt(id.vars = "timeindex")
  
  # Matern 1.5
  matern_1.5_df <- my_res[["matern_1.5"]] %>%
    filter(series == s) %>%
    dplyr::select(rho_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("matern_1.5_", c("rho", "sigma")), "timeindex")) %>% 
    melt(id.vars = "timeindex")
  
  # Matern 2.5
  matern_2.5_df <- my_res[["matern_2.5"]] %>%
    filter(series == s) %>%
    dplyr::select(rho_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("matern_2.5_", c("rho", "sigma")), "timeindex")) %>% 
    melt(id.vars = "timeindex")
  


  cor_res <- lapply(list(ar1_df, oup_df, se_df, matern_1.5_df, matern_2.5_df), function(df) {
    
    
    lapply(unique(df$variable), function(i) {
      
      sub_df <- df %>% filter(variable == i)
      
      cor(sub_df %>% pull(value), sub_df$timeindex, method = "kendall", use = "complete.obs")
      
    }) %>% set_names(unique(df$variable)) %>% unlist
    
    
  }) %>% unlist %>% as.data.frame() %>% t %>% cbind(series = s)
  
  
  return(cor_res)   
  
  
}) %>% do.call(rbind, .) %>% as.data.frame()

tau_res <- tau_res %>%
  mutate(matern1.5_rho = matern_1.5_rho, 
         matern1.5_sigma = matern_1.5_sigma,
         matern2.5_rho = matern_2.5_rho, 
         matern2.5_sigma = matern_2.5_sigma, 
         ar1_rho = ar1_lambda) %>%
  select(-c(matern_1.5_rho, matern_1.5_sigma, matern_2.5_rho, matern_2.5_sigma, ar1_lambda))


tau_res %>% 
  melt_pos_sort(id.vars = "series") %>% 
  separate(col = "variable", sep = "_", into = c("model", "variable"), extra = "merge") %>% 
  ggplot() +
  geom_boxplot(aes(x = variable, y = value, color =model)) + 
  coord_cartesian(ylim = -1:1) +
  geom_hline(yintercept = c(-1, 0,  1), linetype = "dashed")





## Separate ****************************

## Kendall's Tau of all variables. For the Bayesian posteriors mode is used
separate_tau_res <- lapply(unique(ews_res[["oup"]]$series), function(s) {
  
  # Generic parameters
  generic_df <- ews_res[["generic"]] %>% 
    filter(series == s) %>%
    dplyr::select(ar1, sd, sk, kurt, cv, densratio, timeindex) %>% 
    set_colnames(c(paste0("generic_", c("ar1", "sd", "sk", "kurt", "cv", "densratio")), "timeindex"))
  
  # # DDJ
  # ddj_df <- ews_res[["ddj"]] %>% 
  #   filter(series == s) %>%
  #   dplyr::select(conditional_var, total_var, diffusion, jump_intensity, timeindex) %>% 
  #   set_colnames(c(paste0("ddj_", c("conditional_var", "total_var", "diffusion", "jump_intensity")), "timeindex"))
  
  # AR(1)
  ar1_df <- ews_res[["ar1"]] %>% 
    filter(series == s) %>% 
    dplyr::select(lambda_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("ar1_", c("lambda", "sigma")), "timeindex"))
  
  # OUP 
  oup_df <- ews_res[["oup"]] %>% 
    filter(series == s) %>%
    dplyr::select(length_scale_mode, sigma_mode, stat_var_mode, timeindex) %>% 
    set_colnames(c(paste0("OUP_", c("length_scale", "sigma", "stat_var")), "timeindex"))
  
  
  # se 
  se_df <- ews_res[["se"]] %>% 
    filter(series == s) %>%
    dplyr::select(length_scale_mode, sigma_mode, stat_var_mode, timeindex) %>% 
    set_colnames(c(paste0("se_", c("length_scale", "sigma", "stat_var")), "timeindex"))
  
  
  cor_res <- lapply(list(generic_df, se_df, ar1_df, oup_df), function(df) {
    
    timeindex <- df$timeindex
    df <- df %>% dplyr::select(-timeindex)
    
    lapply(colnames(df), function(i) {
      
      cor(df[, i], timeindex, method = "kendall", use = "complete.obs")
      
    }) %>% set_names(colnames(df)) %>% unlist
    
    
  }) %>% unlist %>% as.data.frame() %>% t %>% cbind(series = s)
  
  
  return(cor_res)   
  
  
}) %>% do.call(rbind, .) %>% as.data.frame()

