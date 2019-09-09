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
## Data ***************************************************** ####
set.seed(2424)
times <- seq(from = 0, to = 1000, by = 1)
N_series <- 40

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


## Stan: Separate ******************************************* ####
oup_model <- stan_model("stan_models/oup_single_transition.stan")
ar1_model <- stan_model("stan_models/ar1.stan")
# se_model <- stan_model("se_single_covariance.stan")


iter <- 500
chains <- 1

# CDS set
ar1_ews_res <- ar1_ews_set(ews_residuals, iter, chains)
oup_ews_res <- oup_ews_set(ews_residuals, iter, chains)
# se_ews_res <- se_ews_set(ews_residuals, iter, chains)
generic_ews_res <- generic_ews_set(ews_residuals)


# ar1_store <- ar1_ews_res
# oup_store <- oup_ews_res
# generic_store <- generic_ews_res


# ews_res <- list("ar1" = ar1_ews_res, "oup" = oup_ews_res, "se" = se_ews_res, "generic" = generic_ews_res)
ews_res <- list("ar1" = ar1_ews_res, "oup" = oup_ews_res, "generic" = generic_ews_res)


# saveRDS(object = ews_res, file = "EWS/res/compare_oup_ar1_se_ews_1000_40.RDS")
# ews_res <- readRDS(file = "EWS/res/compare_oup_ar1_ews_1000.RDS")

# # NULL set
# ar1_null_res <- ar1_ews_set(ews_set, iter, chains)
# oup_null_res <- oup_ews_set(ews_set, iter, chains)
# generic_null_res <- generic_ews_set(ews_set)

# null_res <- list("ar1" = ar1_null_res, "oup" = oup_null_res, "generic" = generic_null_res, "ddj" = ddj_null_res)

## Results ************************************************** ####

## Check parameter values are OK
par_res <- lapply(1:N_series, function(s) {
  
  my_res <- ews_res
  
  # # Generic parameters
  generic_df <- my_res[["generic"]] %>%
    filter(series == s) %>%
    dplyr::select(ar1, sd, sk, kurt, cv, densratio, timeindex) %>%
    set_colnames(c(paste0("generic_", c("ar1", "sd", "sk", "kurt", "cv", "densratio")), "timeindex"))%>% 
    melt(id.vars = "timeindex")
  
  # AR(1)
  ar1_df <- my_res[["ar1"]] %>% 
    filter(series == s) %>% 
    dplyr::select(lambda_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("ar1_", c("lambda", "sigma")), "timeindex"))%>% 
    melt(id.vars = "timeindex")
  
  # OUP 
  oup_df <- my_res[["oup"]] %>% 
    filter(series == s) %>%
    dplyr::select(length_scale_mode, sigma_mode, stat_var_mode, timeindex) %>% 
    set_colnames(c(paste0("OUP_", c("length_scale", "sigma", "stat_var")), "timeindex")) %>% 
    melt(id.vars = "timeindex")
  
  
  # SE 
  # se_df <- my_res[["se"]] %>% 
  #   filter(series == s) %>%
  #   dplyr::select(length_scale_mode, stat_var_mode,  timeindex) %>% 
  #   set_colnames(c(paste0("se_", c("length_scale", "stat_var")), "timeindex")) %>% 
  #   melt(id.vars = "timeindex")
  
  return(rbind(generic_df, ar1_df, oup_df) %>% mutate(series = s))   
  
  
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
  generic_df <- my_res[["generic"]] %>%
    filter(series == s) %>%
    dplyr::select(ar1, sd, sk, kurt, cv, densratio, timeindex) %>%
    set_colnames(c(paste0("generic_", c("ar1", "sd", "sk", "kurt", "cv", "densratio")), "timeindex"))

  # AR(1)
  ar1_df <- my_res[["ar1"]] %>% 
    filter(series == s) %>% 
    dplyr::select(lambda_mode, sigma_mode, timeindex) %>% 
    set_colnames(c(paste0("ar1_", c("lambda", "sigma")), "timeindex"))
  
  # OUP 
  oup_df <- my_res[["oup"]] %>% 
    filter(series == s) %>%
    dplyr::select(length_scale_mode, sigma_mode, stat_var_mode, timeindex) %>% 
    set_colnames(c(paste0("OUP_", c("length_scale", "sigma", "stat_var")), "timeindex"))
  
  # 
  # # SE 
  # se_df <- my_res[["se"]] %>% 
  #   filter(series == s) %>%
  #   dplyr::select(length_scale_mode, stat_var_mode, timeindex) %>% 
  #   set_colnames(c(paste0("se_", c("length_scale", "stat_var")), "timeindex"))
  
  
  cor_res <- lapply(list(ar1_df, oup_df, generic_df), function(df) {
    
    timeindex <- df$timeindex
    df <- df %>% dplyr::select(-timeindex)
    
    lapply(colnames(df), function(i) {
      
      cor(df[, i], timeindex, method = "kendall", use = "complete.obs")
      
    }) %>% set_names(colnames(df)) %>% unlist
    
    
  }) %>% unlist %>% as.data.frame() %>% t %>% cbind(series = s)
  
  
  return(cor_res)   
  
  
}) %>% do.call(rbind, .) %>% as.data.frame()

tau_res %>% 
  melt_pos_sort(id.vars = "series") %>% 
  separate(col = "variable", sep = "_", into = c("model", "variable"), extra = "merge") %>%
  ggplot() +
  geom_boxplot(aes(x = variable, y = value, color = model)) + 
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

