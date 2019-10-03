## Gaussian process + generic ews comparison. 
## No hierarchical models
## Data from ecological model
## Sliding windows


## Non-stationary GP data *********************************** ####
non_stat_square_exp_generator_input_process <- stan_model("stan_models/non_stat_square_exp_generator_input_process.stan")


non_stat_oup_generator_input_process <- stan_model("stan_models/non_stat_oup_generator_input_process.stan")


n <- 10
x <- seq(from = 0, to = 100, length.out = 1001)
sigma <- seq(from = .5, to = .75, length.out = length(x))
rho <- seq(from = 1, to = 30, length.out = length(x))
epsilon <- .01

# set.seed(11)
# ews_set_exp <- square_exp_input_set(n, x, rho, sigma, epsilon)
# set.seed(11)
ews_set_oup <- oup_input_set(n, x, rho, sigma, epsilon)

(ews_set_oup)[, i] %>% plot(type = "l")

i <- 2
(ews_set_exp)[, i] %>% generic_ews()
(ews_set_oup)[, i] %>% generic_ews()




# ews_residuals <- get_residuals(ews_set, window = .1)






# residual_plot <- full_join(ews_set %>% set_colnames(c("original", "x")),
#   ews_residuals %>% set_colnames(c("residual", "x"))) %>%
#   mutate(diff = original - residual) %>% 
#   melt(id.vars = "x") %>% 
#   ggplot(aes(x = x, y = value, color = variable)) + 
#   geom_line() + 
#   scale_color_aaas()
# 
# residual_plot


ggplot() +
  geom_line(data = ews_set %>%
              melt(id.vars = "x"), aes(x = x, y = value, group = variable)) 
  # geom_line(data = ews_residuals %>%
  #             melt(id.vars = "x")%>% filter(variable == "V1"), aes(x = x, y = value), color = "blue")




res <- generic_ews_set(ews_set, detrending = "gaussian", window_prop = .5)


res %>% 
  ggplot(aes(x = timeindex, y = ar1)) +
  geom_line(aes(group = series)) + 
  geom_smooth()


## Data ***************************************************** ####
# set.seed(2424)
end <- 100
times <- seq(from = 0, to = end, by = 0.1)
N_series <- 2
n_obs <- 100




# Simulation parameters
r <- 1
K <- 10
cs <- 1
h <- 1
sigma <- 0.03
obs_error_sd <- .01

## Simulate data with CSD ***********************************
ews_set <- lapply(1:N_series, function(j) {
  
  y <- rep(NA, length(times))
  
  # Initial value
  y[1] <- 8
  
  # cs <- 1 + (2.6771 - 1)*(times/max(times))
  cs <- 2 + (2.6771 - 2)*(times/max(times))
  
  
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


# obs_index <- sample((length(ews_residuals$x))/10 %>% floor, n_obs)*10 + sample(-4:5, n_obs, replace = T)


# ews_residuals_obs <- ews_residuals[obs_index, ] %>%
ews_residuals_obs <- ews_set %>%
  as.data.frame() %>%
  set_colnames(c(paste0("V", 1:N_series), "x")) %>% 
  arrange(x)

plot_grid(ews_residuals_obs %>% ggplot(aes(x = x, y = V1)) + geom_line(), ews_residuals %>% ggplot(aes(x = x, y = V1)) + geom_line(), ncol= 1)





## Vostok data ********************************************** ####
vostok <- readRDS(file = "EWS/data/vostok_pre.RDS")

vostok_residuals <- function(df, inter = FALSE) {
  
  mydf <- df
  
  ## Join observations
  
  filtered_dfs <- lapply(1:4, function(i) {
    
    series <- mydf %>%
      filter(Era == i, interpolation == inter) %>% 
      select(Delta, Age) %>% 
      set_colnames(c(paste0("Era", i), "x"))
    
    # Smooth
    smoothed <- smth(series[, 1], window = .1, tails = TRUE)
    series[, 1] <- series[, 1] - smoothed
    
    series
  })
  
  joined <- full_join(filtered_dfs[[1]], filtered_dfs[[2]], "x")
  joined <- full_join(joined,  filtered_dfs[[3]], "x")
  joined <- full_join(joined,  filtered_dfs[[4]], "x")
  
  joined <- joined %>% 
    select(Era1, Era2, Era3, Era4, x)
 
  return(joined)
   
}

vostok_data <- vostok_residuals(vostok, inter = FALSE)
imputed_vostok_data <- vostok_residuals(vostok, inter = TRUE)


## Stan: Separate ******************************************* ####
oup_model <- stan_model("stan_models/oup_fitter.stan")
se_model <- stan_model("stan_models/square_exp_fitter.stan")
non_stat_se_model <- stan_model("stan_models/non_stat_square_exp_fitter.stan")
non_stat_oup_model <- stan_model("stan_models/non_stat_oup_fitter.stan")
matern_1.5_model <- stan_model("stan_models/matern_1.5_fitter.stan")
matern_2.5_model <- stan_model("stan_models/matern_2.5_fitter.stan")
ar1_model <- stan_model("stan_models/ar1.stan")
ar1_missing_model <- stan_model("stan_models/ar1_missing.stan")


iter <- 500
chains <- 1



# CDS set
ar1_ews_res <- ar1_ews_set(imputed_vostok_data %>% select(-c(Era1)), iter, chains, window_length = 25)

generic_ews_res <- generic_ews_set(ews_residuals, window_prop = .5)


# ar1_ews_res <- ar1_missing_ews_set(ews_residuals_obs, iter = 500, chains, window_length = 50)

oup_ews_res <- GP_ews_set(vostok_data %>% select(-c(Era1)), iter = 500, chains = 1, window_length = 25, kernel = "OUP")
se_ews_res <- GP_ews_set(vostok_data %>% select(-c(Era1)), iter = 500, chains = 1, window_length = 25, kernel = "se")
matern_1.5_ews_res <- GP_ews_set(vostok_data %>% select(-c(Era1)), iter = 500, chains = 1, window_length = 25, kernel = "matern_1.5")
matern_2.5_ews_res <- GP_ews_set(vostok_data %>% select(-c(Era1)), iter = 500, chains = 1, window_length = 25, kernel = "matern_2.5")





# ews_res <- list("ar1" = ar1_ews_res, "oup" = oup_ews_res, "se" = se_ews_res, "generic" = generic_ews_res)
ews_res <- list("ar1" = ar1_ews_res,
                "oup" = oup_ews_res,
                "se" = se_ews_res, 
                "matern_1.5" = matern_1.5_ews_res, 
                "matern_2.5" = matern_2.5_ews_res)


# saveRDS(object = ews_res, file = "EWS/res/vostok_results.RDS")
# ews_res <- readRDS(file = "EWS/res/GP_EWS_one_long_1000_NA.RDS")
ews_res <- readRDS(file = "EWS/res/vostok_results.RDS")

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
  
  ar1_df <- ar1_df %>% mutate(timeindex = timeindex - min(timeindex))
  
  
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
  mutate(series = as.factor(series)) %>% 
  ggplot(aes(x = timeindex, y = value, color = variable)) + 
  geom_line(aes(group = series, linetype = series)) + 
  facet_wrap(~variable, scales = "free")


par_res %>%
  filter(series == 3) %>%
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
  geom_hline(yintercept = c(-1, 0,  1), linetype = "dashed") + 
  facet_wrap(~series)





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


## Non-stationary ******************************************* ####

non_stat_res <- ns_GP_ews_set(vostok %>% filter(Era == 3, interpolation == FALSE) %>% select(Delta, Age) %>% arrange(Age) %>% set_colnames(c("Delta", "x")), iter = iter, chains = 1, n_samples = 50)
# non_stat_res <- ns_GP_ews_set(vostok_data[, c(3,5)], iter = iter, chains = 1, n_samples = 50)
# non_stat_res <- ns_GP_ews_set(vostok_data[, c(2,5)], iter = iter, chains = 1, n_samples = 50)

non_stat_res_summary <- non_stat_res[["summary"]]
non_stat_res_samples <- non_stat_res[["samples"]]



# Plot results
p <- ggplot()+
geom_line(data = non_stat_res_samples %>% 
              melt(id.vars = c("x", "parameter", "series")) %>%
              mutate(variable = gsub("V", "", variable) %>% as.numeric()), 
            aes(x = x, y = value, group = variable), color = "blue", alpha = .25) + 
  geom_line(data = non_stat_res_summary, aes(x = x, y = mode)) +
facet_wrap(parameter~series, scales = "free")



non_stat_taus <- lapply(c("rho", "sigma"), function(par) {
  
  mydf <- non_stat_res_samples %>% 
    filter(parameter == par) %>% 
    select(-c(series, parameter))
  
  x <- mydf %>% pull(x)
  mydf <- mydf %>% select(-x)
  
  cors <- lapply(1:ncol(mydf), function(i) {
    
    cor(mydf[, i], x)
    
  }) %>% unlist
  
  return(cors)
}) %>% do.call(cbind, .) %>% set_colnames(c("rho", "sigma"))


q <- non_stat_taus %>% 
  melt %>% 
  ggplot(aes(x = value)) + 
  geom_histogram() + 
  geom_density() +
  facet_wrap(~Var2)


plot_grid(pp, pp_res, p, q, ncol = 1)




mode_bacf <- lapply(1, function(dum) {
  
  rho <- non_stat_res_summary %>% 
    filter(parameter == "rho") %>% 
    pull(mode)
  
  sigma <- non_stat_res_summary %>% 
    filter(parameter == "sigma") %>% 
    pull(mode)
  
  pars <- data.frame(rho, sigma)
  
  x_grid <- seq(from = 0, to = 20, length.out = 1000)
  
  mydf <- lapply(1:nrow(pars), function(i) {
    
    exp(-x_grid/pars[i, 1])
    
  }) %>% do.call(cbind, .) %>%
    set_colnames(paste0("T", 1:120)) %>% 
    cbind(x_grid)
  
  return(mydf)
})[[1]]

mode_bacf %>%
  as.data.frame() %>%
  melt(id.vars = "x_grid") %>%
  mutate(time = gsub("T", "", variable) %>% as.numeric()) %>% 
  ggplot(aes(y = x_grid, x = time, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = .9)



bacf <- lapply(1, function(i) {
  
  rho_df <- non_stat_res_samples %>% 
    filter(parameter == "rho") %>% 
    select(-c(x, parameter, series))
  
  sigma_df <- non_stat_res_samples %>% 
    filter(parameter == "sigma") %>% 
    select(-c(x, parameter, series))
  
  exp(-1/rho_df)
  
}) %>% do.call(cbind, .) %>% as.data.frame()

