## Hierarchical OUP/AR1 sliding window ews

## Functions ************************************************ ####
melt_pos_sort <- function(df, id.vars) {
  
  my_df <- df
  
  cmeans <- my_df %>%
    select(-one_of(id.vars)) %>%
    colMeans()
  
  my_df <- lapply(1:length(cmeans), function(i) {
    
    x <- (my_df %>% select(-one_of(id.vars)))[, i]
    
    if(cmeans[i] < 0) {
      return(-x)
    } else {
      return(x)
    }
    
  }) %>% 
    do.call(cbind, .) %>% 
    cbind(my_df[, id.vars]) %>% 
    set_colnames(colnames(my_df))
  
  
  melt_df <- my_df %>%
    as.data.frame() %>% 
    melt(id.vars = id.vars) %>% 
    mutate(variable = factor(variable, levels = names(sort(abs(cmeans)))))
  
  
  return(melt_df)
  
}
list_bind <- function(df_list) {
  
  element_length <- sapply(1:length(df_list), function(i) {
    df_list[[i]] %>% length()
  })
  
  
  if(length(unique(element_length)) != 1) {
    stop("List elements are not of equal length")
  }
  
  
  res <- lapply(1:(element_length)[1], function(i)  {
    lapply(1:length(df_list), function(j) {
      
      df_list[[j]][[i]]
      
    }) %>% 
      do.call(rbind, .)
  })
  
  return(res)
  
}
get_population_distribution <- function(df, par, model = NULL, distribution, padding = 2.5) {
  
  my_df <- population_par_res %>% filter(parameter == par)
  
  if(!is.null(model)) {
    
    my_df <- my_df %>% filter(model == model)
    
  }
  
  # Get distributions
  if(distribution == "normal") {
    
    dist_df <- lapply(1:nrow(my_df), function(i) {
      
      p <- my_df[i, c("a_mode", "b_mode")] %>% as.matrix()
      
      x <- seq(from = p[1] - padding*p[2], to = p[1] + padding*p[2], by = .1)
      y <- dnorm(x, p[1], p[2])
      
      data.frame(x = x, y = y, window = i, parameter = par)
    }) %>% do.call(rbind, .)
    
  } else if(distribution == "gamma") {
    
    dist_df <- lapply(1:nrow(my_df), function(i) {
      
      p <- my_df[i, c("a_mode", "b_mode")] %>% as.matrix()
      
      x <- seq(from = 0, to = 5, by = .1)
      y <- dgamma(x, p[1], p[2])
      
      data.frame(x = x, y = y, window = i, parameter = par)
    }) %>% do.call(rbind, .)
    
  }
  
  return(dist_df)
}

## Data ***************************************************** ####
set.seed(2)
times <- seq(from = 0, to = 100, by = 1)
N_series <- 20

# Simulation parameters
r <- 1
K <- 10
cs <- 1
h <- 1
sigma <- 0.03

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
  
  y
  # return(log(y + 1))
  
  
}) %>%
  do.call(cbind, .) %>%
  cbind(x = times) %>%
  as.data.frame()

ews_set <- ews_set[1:(which(((1 + (2.6771 - 1)*(times/max(times))) < 2.604)) %>% max()), ]

## Null data
N_null_series <- 10
null_set <- lapply(1:N_null_series, function(j) {
  
  r <- 1
  K <- 10
  cs <- 1
  h <- 1
  sigma <- 0.03
  
  
  y <- rep(NA, length(times))
  y[1] <- 8
  
  # cs <- 1 + (2.6771 - 1)*(times/max(times))
  
  
  for(i in 2:length(times)) {
    
    dt <- times[(i-1):i]
    
    finer_grid <- seq(from = dt[1], to = dt[2], by = 0.01)
    
    y[i] <- ews_generator(y[i-1], finer_grid, c = cs, milstein = T)[length(finer_grid)]
    
  }
  
  y
  # return(log(y + 1))
  
  
}) %>%
  do.call(cbind, .) %>%
  cbind(x = times) %>%
  as.data.frame()


null_set <- null_set[1:(which(((1 + (2.6771 - 1)*(times/max(times))) < 2.604)) %>% max()), ]


## Residual data
ews_residuals <- lapply(1:N_series, function(i) {
  
  series <- ews_set[, i]
  smoothed <- smth(series, window = .1, tails = TRUE)
  res <- series - smoothed
  return(res)
}) %>% do.call(cbind, .) %>% cbind(., ews_set$x)
null_residuals <- lapply(1:N_null_series, function(i) {
  
  series <- null_set[, i]
  smoothed <- smth(series, window = .1, tails = TRUE)
  res <- series - smoothed
  return(res)
}) %>% do.call(cbind, .) %>% cbind(., ews_set$x)



## Stan: Hierarchical *************************************** ####
ar1_hierarchical <- stan_model("stan_models/ar1_hierarchical.stan")
oup_hierarchical <- stan_model("stan_models/oup_hierarchical_transition.stan")


ar1_hierarchical_res <- ar1_hierarchical_ews(ews_residuals[, 1:N_series],
                     window_prop = .5,
                     iter = 500, chains = 1)

oup_hierarchical_res <- oup_hierarchical_ews(ews_residuals[, 1:N_series],
                                             window_prop = .5,
                                             iter = 500, chains = 1)


hierarchical_res <- list("ar1" = ar1_hierarchical_res,
                         "oup" = oup_hierarchical_res)

saveRDS(object = hierarchical_res, file = "EWS/short100_hierarchical_results.rds")

## Stan: Separate ******************************************* ####
oup_model <- stan_model("stan_models/oup_single_transition.stan")
ar1_model <- stan_model("stan_models/ar1.stan")


iter <- 500
chains <- 1

# CDS set
ar1_ews_res <- ar1_ews_set(ews_set, iter, chains)
oup_ews_res <- oup_ews_set(ews_set, iter, chains)
generic_ews_res <- generic_ews_set(ews_set)
ddj_ews_res <- ddj_ews_set(ews_residuals)

ews_res <- list("ar1" = ar1_ews_res, "oup" = oup_ews_res, "generic" = generic_ews_res, "ddj" = ddj_ews_res)

saveRDS(object = ews_res, file = "EWS/short100_EWS_results_10.rds")
# ews_res <- readRDS( file = "EWS/short100_EWS_results_20.rds")


# NULL set
ar1_null_res <- ar1_ews_set(ews_set, iter, chains)
oup_null_res <- oup_ews_set(ews_set, iter, chains)
generic_null_res <- generic_ews_set(ews_set)
ddj_null_res <- ddj_ews_set(null_residuals)

null_res <- list("ar1" = ar1_null_res, "oup" = oup_null_res, "generic" = generic_null_res, "ddj" = ddj_null_res)

# saveRDS(object = null_res, file = "EWS/short100_NULL_results_10.rds")
# null_res <- readRDS(file = "EWS/short100_NULL_results_30.rds")


## Results ************************************************** ####

## Hierarchical ************************

# Kendall's Tau of all variables. For the Bayesian posteriors mode is used.
hierarchical_tau_res <- lapply(1:N_series, function(s) {
  
  my_res <- hierarchical_res
  
  # # Generic parameters
  # generic_df <- my_res[["generic"]] %>%
  #   filter(series == s) %>%
  #   dplyr::select(ar1, sd, sk, kurt, cv, densratio, timeindex) %>%
  #   set_colnames(c(paste0("generic_", c("ar1", "sd", "sk", "kurt", "cv", "densratio")), "timeindex"))
  # 
  # # DDJ
  # ddj_df <- my_res[["ddj"]] %>%
  #   filter(series == s) %>%
  #   dplyr::select(conditional_var, total_var, diffusion, jump_intensity, timeindex) %>%
  #   set_colnames(c(paste0("ddj_", c("conditional_var", "total_var", "diffusion", "jump_intensity")), "timeindex"))
  # 
   
  # AR(1)
  ar1_df <- my_res[["ar1"]][["parameters"]] %>% 
    filter(series == s) %>% 
    dplyr::select(lambda_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("ar1_", c("lambda", "sigma")), "timeindex"))
  
  # OUP 
  oup_df <- my_res[["oup"]][["parameters"]] %>% 
    filter(series == s) %>%
    dplyr::select(length_scale_mode, sigma_mode, stat_var_mode, timeindex) %>% 
    set_colnames(c(paste0("OUP_", c("length_scale", "sigma", "stat_var")), "timeindex"))
  
  
  cor_res <- lapply(list(ar1_df, oup_df), function(df) {
    
    timeindex <- df$timeindex
    df <- df %>% dplyr::select(-timeindex)
    
    lapply(colnames(df), function(i) {
      
      cor(df[, i], timeindex, method = "kendall", use = "complete.obs")
      
    }) %>% set_names(colnames(df)) %>% unlist
    
    
  }) %>% unlist %>% as.data.frame() %>% t %>% cbind(series = s)
  
  
  return(cor_res)   
  
  
}) %>% do.call(rbind, .) %>% as.data.frame()

hierarchical_tau_res %>% 
  melt_pos_sort(id.vars = "series") %>% 
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot() + 
  labs(subtitle = "Hierarchical models, Kendall tau for posterior modes. \n 100 obs, 20 series, dt = 1.") + 
  coord_cartesian(ylim = -1:1) +
  geom_hline(yintercept = c(-1, 0,  1), linetype = "dashed")



## Population distribution parameters
population_par_res <- lapply(1, function(i) {

  oup_pars <- c("length_scale", "stat_var", "sigma")
  oup_population_pars <- lapply(oup_pars, function(p) {
    # print(p)
    my_df <- hierarchical_res[["oup"]][["hyperparameters"]]

    res <- my_df[, c(grep(p, colnames(my_df), value = T), "window")]

    # Remove parameter from column names
    colnames(res) <- gsub(paste0(p, "_"), "", colnames(res))

    # Add paramemeter name
    res <- res %>% mutate(parameter = p, model = "oup")

    return(res)

  }) %>% do.call(rbind, .)
  ar1_population_pars <- lapply(oup_pars, function(p) {

    my_df <- hierarchical_res[["ar1"]][["hyperparameters"]]

    res <- my_df[, c(grep(p, colnames(my_df), value = T), "window")]

    # Remove parameter from column names
    colnames(res) <- gsub(paste0(p, "_"), "", colnames(res))

    # Add paramemeter name
    res <- res %>% mutate(parameter = p, model = "ar1")

    return(res)

  }) %>% do.call(rbind, .)


  res <- rbind(ar1_population_pars, oup_population_pars)

  return(res)

})[[1]]



## Separate ****************************

## Kendall's Tau of all variables. For the Bayesian posteriors mode is used
separate_tau_res <- lapply(unique(ews_res[["oup"]]$series), function(s) {
  
  # Generic parameters
  generic_df <- ews_res[["generic"]] %>% 
    filter(series == s) %>%
    dplyr::select(ar1, sd, sk, kurt, cv, densratio, timeindex) %>% 
    set_colnames(c(paste0("generic_", c("ar1", "sd", "sk", "kurt", "cv", "densratio")), "timeindex"))
  
  # DDJ
  ddj_df <- ews_res[["ddj"]] %>% 
    filter(series == s) %>%
    dplyr::select(conditional_var, total_var, diffusion, jump_intensity, timeindex) %>% 
    set_colnames(c(paste0("ddj_", c("conditional_var", "total_var", "diffusion", "jump_intensity")), "timeindex"))
  
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
  
  
  cor_res <- lapply(list(generic_df, ddj_df, ar1_df, oup_df), function(df) {
    
    timeindex <- df$timeindex
    df <- df %>% dplyr::select(-timeindex)
    
    lapply(colnames(df), function(i) {
      
      cor(df[, i], timeindex, method = "kendall", use = "complete.obs")
      
    }) %>% set_names(colnames(df)) %>% unlist
    
    
  }) %>% unlist %>% as.data.frame() %>% t %>% cbind(series = s)
  
  
  return(cor_res)   
  
  
}) %>% do.call(rbind, .) %>% as.data.frame()





## Figures ************************************************** ####

pos_sort <- function(df, id.vars) {
  
  my_df <- df
  
  cmeans <- my_df %>%
    select(-one_of(id.vars)) %>%
    colMeans()
  
  my_df <- lapply(1:length(cmeans), function(i) {
    
    x <- (my_df %>% select(-one_of(id.vars)))[, i]
    
    if(cmeans[i] < 0) {
      return(-x)
    } else {
      return(x)
    }
    
  }) %>% 
    do.call(cbind, .) %>% 
    cbind(my_df[, id.vars]) %>% 
    set_colnames(colnames(my_df))
  
}



## Tau *********************

# hierachical and separate oup/ar1 parameter comparison
rbind(cbind(hierarchical_tau_res,
            model = "hierarchical") %>% pos_sort(id.vars = c("series", "model")),
      cbind(separate_tau_res[ ,colnames(hierarchical_tau_res)],
            model = "separate")  %>% pos_sort(id.vars = c("series", "model"))) %>% 
  melt(id.vars = c("series", "model")) %>% 
  ggplot() +
  geom_boxplot(aes(x = variable, y = value, color = model), position = "dodge") +
  labs(x = "", y = expression(~tau), subtitle = "Kendall tau in hierarchical and individual models") +
  scale_color_startrek() + 
  geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed")





## Parameter values and evolution
# OUP
oup_pars <- c("length_scale", "stat_var")
oup_par_comparison_df <- lapply(oup_pars, function(p){
  
  
  hier_df <- hierarchical_res[["oup"]][["parameters"]] %>%
    select(timeindex, series, paste0(rep(p, each = 3),
                                     c("_mode", "_lower2.5", "_upper97.5")))
  colnames(hier_df) <- hier_df %>% colnames %>% gsub(paste0(p, "_"), "", .)
  
  
  sep_df <- ews_res[["oup"]] %>%
    select(timeindex,  series,paste0(rep(p, each = 3),
                                     c("_mode", "_lower2.5", "_upper97.5")))
  colnames(sep_df) <- sep_df %>% colnames %>% gsub(paste0(p, "_"), "", .)
  
  
  df <- rbind(hier_df %>% mutate(model = "hierarchical"), 
  sep_df %>% mutate(model = "separate")) %>% 
    mutate(parameter = p)
  
  
  return(df)
}) %>% do.call(rbind, .)
oup_par_comparison_df %>% 
  ggplot() + 
  geom_ribbon(aes(x = timeindex, ymin = lower2.5, ymax = upper97.5, fill = model), alpha = .25) +
  geom_line(aes(x = timeindex, y = mode, group = model)) +
  facet_wrap(c("series", "parameter"), scale = "free", ncol = 4)



ar1_pars <- c("lambda", "sigma")
ar1_par_comparison_df <- lapply(ar1_pars, function(p){
  
  
  hier_df <- hierarchical_res[["ar1"]][["parameters"]] %>%
    select(timeindex, series, paste0(rep(p, each = 3),
                                     c("_mode", "_lower2.5", "_upper97.5")))
  colnames(hier_df) <- hier_df %>% colnames %>% gsub(paste0(p, "_"), "", .)
  
  
  sep_df <- ews_res[["ar1"]] %>%
    select(timeindex,  series,paste0(rep(p, each = 3),
                                     c("_mode", "_lower2.5", "_upper97.5")))
  colnames(sep_df) <- sep_df %>% colnames %>% gsub(paste0(p, "_"), "", .)
  
  
  df <- rbind(hier_df %>% mutate(model = "hierarchical"), 
              sep_df %>% mutate(model = "separate")) %>% 
    mutate(parameter = p)
  
  
  return(df)
}) %>% do.call(rbind, .)
ar1_par_comparison_df %>% 
  ggplot() + 
  geom_ribbon(aes(x = timeindex, ymin = lower2.5, ymax = upper97.5, fill = model), alpha = .25) + 
  geom_line(aes(x = timeindex, y = mode, group = model)) +
  facet_wrap(c("series", "parameter"), scale = "free", ncol = 4)


rbind(oup_par_comparison_df %>% mutate(parameter = paste0(parameter, "_oup")),
      ar1_par_comparison_df %>% mutate(parameter = paste0(parameter, "_ar1_bayes"))) %>%
  ggplot(aes(x = parameter, y = mode, color = model), position = "dodge") +
  geom_boxplot() + 
  labs(subtitle = "Parameter values from ews time series")




## Population parameters
population_par_res %>% 
  ggplot() +
  geom_errorbar(aes(x = window, ymin = a_lower2.5, ymax = a_upper97.5, color = parameter)) + 
  facet_wrap(~model, scales = "free")

