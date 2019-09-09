## Functions ************************************************ ####
std_edit <- function(df, vars) {
  
  if(any(!(vars %in% colnames(df)))) {
    stop("Variables not found in the data frame")
  }
  
  
  res <- lapply(vars, function(p) {
    # print(p)
    x <- df[, p]
    
    # Center and normalize
    x <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
    
    # First to zero
    x <- x - x[1]
    
    # Last non NA to positive
    if(x[last(which(!is.na(x)))] < 0) {
      x <- -x
    }
    
    return(x)
    
  }) %>% 
    do.call(cbind, .) %>% 
    set_colnames(paste0(vars, "_edit"))
  
  
  res %>% data.frame()
}


oup_model <- stan_model("stan_models/oup_single_transition.stan")
ar1_model <- stan_model("stan_models/ar1.stan")


iter <- 500
chains <- 1

## Data ***************************************************** ####
set.seed(2)
times <- seq(from = 0, to = 100, by = 1)
N_series <- 10

# Simulation parameters
r <- 1
K <- 10
cs <- 1
h <- 1
sigma <- 0.03

## Simulate data with CSD
ews_set <- lapply(1:N_series, function(j) {
  
  
  
  
  y <- rep(NA, length(times))
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


## EWS ****************************************************** ####

# CDS set
ar1_ews_res <- ar1_ews_set(ews_set, iter, chains)
oup_ews_res <- oup_ews_set(ews_set, iter, chains)
generic_ews_res <- generic_ews_set(ews_set)
ddj_ews_res <- ddj_ews_set(ews_residuals)

ews_res <- list("ar1" = ar1_ews_res, "oup" = oup_ews_res, "generic" = generic_ews_res, "ddj" = ddj_ews_res)

# saveRDS(object = ews_res, file = "EWS/short100_EWS_results_10.rds")
# ews_res <- readRDS( file = "EWS/short100_EWS_results.rds")


# NULL set
ar1_null_res <- ar1_ews_set(ews_set, iter, chains)
oup_null_res <- oup_ews_set(ews_set, iter, chains)
generic_null_res <- generic_ews_set(ews_set)
ddj_null_res <- ddj_ews_set(null_residuals)

null_res <- list("ar1" = ar1_null_res, "oup" = oup_null_res, "generic" = generic_null_res, "ddj" = ddj_null_res)

# saveRDS(object = null_res, file = "EWS/short100_NULL_results_10.rds")
# null_res <- readRDS(file = "EWS/short100_NULL_results_10.rds")



# # EWS for NULL set
# compare_null_df <- lapply(1:N_null_series, function(ind) {
#   
#   print(ind)
#   
#   ## Extract data *************************************
#   series <- null_residuals[, ind]
#   
#   ## Generic ******************************************
#   generic <- generic_ews(series)
#   dev.off()
#   
#   generic_pars <- c("ar1", "sd", "sk", "kurt", "cv", "densratio")
#   generic <- generic[, c("timeindex", generic_pars)]
#   
#   
#   
#   # Standardize
#   # generic <- lapply(generic_pars, function(p) {
#   #   
#   #   (generic[, p] - mean(generic[, p]))/sd(generic[, p])
#   #   
#   # }) %>% 
#   #   do.call(cbind, .) %>% 
#   #   set_colnames(paste0(generic_pars, "_std")) %>% 
#   #   cbind(generic, .)
#   #   
#   
#   
#   
#   # Edit further: all start at zero and end in a >0 value
#   # edit_generic <- generic[, paste0(generic_pars, "_std")]
#   # for(i in 1:ncol(edit_generic)) {
#   #   
#   #   x <- edit_generic[, i]
#   #   
#   #   x <- x - x[1]
#   #   
#   #   if(x[length(x)] < 0 ) {
#   #     
#   #     x <- -x
#   #   }
#   #   
#   #   edit_generic[, i] <- x
#   # }
#   # 
#   # # Add time
#   # edit_generic <- data.frame(edit_generic, timeindex = generic$timeindex)
#   
#   # edit_generic %>%
#   #   melt(id.vars = "timeindex") %>%
#   #   ggplot(aes(x = timeindex, y = value, color = variable)) +
#   #   geom_line() +
#   #   scale_color_locuszoom()
#   # 
#   
#   # DDJ ********************************************************
#   # ddj <- ddjnonparam_ews(data.frame(ews_set$x, series))
#   # 
#   # 
#   # # Take timeindex and measures dependent on time
#   # ddj <- ddj[names(ddj)[6:10]] %>%
#   #   do.call(cbind, .) %>%
#   #   as.data.frame() %>% 
#   #   set_colnames(c("timeindex", "conditional_var", "total_var", "diffusion", "jump_intensity"))
#   #
#   # dev.off()
#   # 
#   # ddj_pars <- c("conditional_var", "total_var", "diffusion", "jump_intensity")
# 
#   # Standardize
#   # ddj <- lapply(ddj_pars, function(p) {
#   # 
#   #   (ddj[, p] - mean(ddj[, p]))/sd(ddj[, p])
#   # 
#   # }) %>%
#   #   do.call(cbind, .) %>%
#   #   cbind(ddj, .) %>%
#   #   set_colnames(c("timeindex", ddj_pars, paste0(ddj_pars, "_std")))
# 
#   #
#   # # Edit more
#   # edit_ddj <- ddj[, paste0(ddj_pars, "_std")]
#   # for(i in 1:ncol(edit_ddj)) {
#   #
#   #   x <- edit_ddj[, i]
#   #
#   #   x <- x - x[1]
#   #
#   #   if(x[length(x)] < 0 ) {
#   #
#   #     x <- -x
#   #   }
#   #
#   #   edit_ddj[, i] <- x
#   # }
#   #
#   # edit_ddj <- data.frame(edit_ddj, timeindex = ews_set$x)
# 
#   # edit_ddj %>%
#   #   melt(id.vars = "timeindex") %>%
#   #   ggplot(aes(x = timeindex, y = value, color = variable)) +
#   #   geom_line(alpha = .75) +
#   #   scale_color_locuszoom()
#   
#   
#   ## BDS **********************************************************
#   #
#   # (Skip for now, the function doesn't work)
#   # 
#   # bds <- bdstest_ews(data.frame(ews_set$x, series))
#   
#   
#   ## Sensitivity **************************************************
#   # sensi <- sensitivity_ews(data.frame(ews_set$x, series))
#   
#   ## Surrogates
#   # surrogates <- surrogates_ews(data.frame(ews_set$x, series))
#   
#   
#   ## Conditional heteroscedasticity
#   # ch <- ch_ews(data.frame(ews_set$x, series))
#   
#   
#   ## Bayes AR(1) *************************************************
#   ar1 <- ar1_ews(series)
#   
#   ar1 <- ar1 %>% mutate(timeindex = ar1$to)
#   
#   
#   # ar1 <- ar1 %>% 
#   #   mutate(stat_var_ar1 = sigma_ar1^2/(1 + lambda_ar1), 
#   #   ) %>% mutate(timeindex = ar1$to)
#   
#   
#   ## OUP *********************************************************
#   
#   oup <- oup_ews(series)
#   
#   oup <- oup %>% mutate(timeindex = oup$to)
#   
#   
#   # oup <- oup %>% 
#   #   mutate(stat_var = sigma^2/(2*lambda), 
#   #   ) %>% mutate(timeindex = oup$to)
#   
#   
#   # oup_edit <- std_edit(oup, c("likelihood", "lambda", "mu", "sigma", "resilience", "iqr")) %>% mutate(timeindex = oup$to)
#   
#   # oup_edit %>%
#   #   melt(id.vars = "timeindex") %>% 
#   #   ggplot(aes(x = timeindex, y = value, color = variable)) + 
#   #   geom_line() + 
#   #   scale_color_locuszoom()
#   
#   
#   ## Combine *****************************************************
#   
#   # oup %>% dim
#   # edit_ddj %>% dim
#   # edit_generic %>% dim
#   
#   
#   
#   # combined_df <- full_join(oup_edit, edit_generic, "timeindex") %>%
#   #   mutate(series = ind)
#   
#   
#   combined_df <- full_join(oup, generic, "timeindex")
#   combined_df <- full_join(combined_df, ar1, c("timeindex", "from", "to", "window")) %>% 
#     mutate(series = ind)
#   
#   
#   ## Remove .x and .y from colnames
#   colnames(combined_df)[grep("\\.y", colnames(combined_df))] <- grep("\\.y", colnames(combined_df), value = T) %>% gsub("\\.y", "_ar1", .)
#   colnames(combined_df)[grep("\\.x", colnames(combined_df))] <- grep("\\.x", colnames(combined_df), value = T) %>% gsub("\\.x", "", .)
#   
#   
#   return(combined_df)
#   
# }) %>% do.call(rbind, .) 
# 
# # compare_null_df <- rbind(compare_null_df, xx)
# # xx <- compare_null_df
# 
# compare_null_df <- compare_null_df %>% 
#   mutate(data = "null")
# 
# saveRDS(object = compare_null_df, file = "EWS/EWS_null_results_20")
# 
# 
# # EWS for CSD set
# compare_ews_df <- lapply(1:N_series, function(ind) {
#   
#   print(ind)
#   
#   ## Extract data *************************************
#   series <- ews_residuals[, ind]
#   
#   ## Generic ******************************************
#   generic <- generic_ews(series)
#   dev.off()
#   
#   generic_pars <- c("ar1", "sd", "sk", "kurt", "cv", "densratio")
#   generic <- generic[, c("timeindex", generic_pars)]
#   
#   
#   
#   # Standardize
#   # generic <- lapply(generic_pars, function(p) {
#   #   
#   #   (generic[, p] - mean(generic[, p]))/sd(generic[, p])
#   #   
#   # }) %>% 
#   #   do.call(cbind, .) %>% 
#   #   set_colnames(paste0(generic_pars, "_std")) %>% 
#   #   cbind(generic, .)
#   #   
#   
#   
#   
#   # Edit further: all start at zero and end in a >0 value
#   # edit_generic <- generic[, paste0(generic_pars, "_std")]
#   # for(i in 1:ncol(edit_generic)) {
#   #   
#   #   x <- edit_generic[, i]
#   #   
#   #   x <- x - x[1]
#   #   
#   #   if(x[length(x)] < 0 ) {
#   #     
#   #     x <- -x
#   #   }
#   #   
#   #   edit_generic[, i] <- x
#   # }
#   # 
#   # # Add time
#   # edit_generic <- data.frame(edit_generic, timeindex = generic$timeindex)
#   
#   # edit_generic %>%
#   #   melt(id.vars = "timeindex") %>%
#   #   ggplot(aes(x = timeindex, y = value, color = variable)) +
#   #   geom_line() +
#   #   scale_color_locuszoom()
#   # 
#   
#   ## DDJ ********************************************************
#   # ddj <- ddjnonparam_ews(data.frame(ews_set$x, series))
#   # 
#   # 
#   # # Take timeindex and measures dependent on time
#   # ddj <- ddj[names(ddj)[6:10]] %>% 
#   #   do.call(cbind, .) %>% 
#   #   set_colnames(c("timeindex", "conditional_var", "total_var", "diffusion", "jump_intensity"))
#   # 
#   # dev.off()
#   # 
#   # ddj_pars <- c("conditional_var", "total_var", "diffusion", "jump_intensity")
#   # 
#   # # Standardize
#   # ddj <- lapply(ddj_pars, function(p) {
#   #   
#   #   (ddj[, p] - mean(ddj[, p]))/sd(ddj[, p])
#   #   
#   # }) %>% 
#   #   do.call(cbind, .) %>% 
#   #   cbind(ddj, .) %>% 
#   #   set_colnames(c("timeindex", ddj_pars, paste0(ddj_pars, "_std")))
#   # 
#   # 
#   # # Edit more
#   # edit_ddj <- ddj[, paste0(ddj_pars, "_std")]
#   # for(i in 1:ncol(edit_ddj)) {
#   #   
#   #   x <- edit_ddj[, i]
#   #   
#   #   x <- x - x[1]
#   #   
#   #   if(x[length(x)] < 0 ) {
#   #     
#   #     x <- -x
#   #   }
#   #   
#   #   edit_ddj[, i] <- x
#   # }
#   # 
#   # edit_ddj <- data.frame(edit_ddj, timeindex = ews_set$x)
#   
#   # edit_ddj %>%
#   #   melt(id.vars = "timeindex") %>% 
#   #   ggplot(aes(x = timeindex, y = value, color = variable)) + 
#   #   geom_line(alpha = .75) + 
#   #   scale_color_locuszoom()
#   
#   
#   ## BDS **********************************************************
#   #
#   # (Skip for now, the function doesn't work)
#   # 
#   # bds <- bdstest_ews(data.frame(ews_set$x, series))
#   
#   
#   ## Sensitivity **************************************************
#   # sensi <- sensitivity_ews(data.frame(ews_set$x, series))
#   
#   ## Surrogates
#   # surrogates <- surrogates_ews(data.frame(ews_set$x, series))
#   
#   
#   ## Conditional heteroscedasticity
#   # ch <- ch_ews(data.frame(ews_set$x, series))
#   
#   
#   ## Bayes AR(1) *************************************************
#   ar1 <- ar1_ews(series)
#   
#   ar1 <- ar1 %>% mutate(timeindex = ar1$to)
#   
#   # ar1 <- ar1 %>% 
#   #   mutate(stat_var_ar1 = sigma_ar1^2/(1 + lambda_ar1), 
#   #   ) %>% mutate(timeindex = ar1$to)
#   
#   
#   ## OUP *********************************************************
#   
#   oup <- oup_ews(series)
#   
#   oup <- oup %>% mutate(timeindex = oup$to)
#   
#   # oup <- oup %>% 
#   #   mutate(stat_var = sigma^2/(2*lambda), 
#   #   ) %>% mutate(timeindex = oup$to)
#   
#   
#   # oup_edit <- std_edit(oup, c("likelihood", "lambda", "mu", "sigma", "resilience", "iqr")) %>% mutate(timeindex = oup$to)
#   
#   # oup_edit %>%
#   #   melt(id.vars = "timeindex") %>% 
#   #   ggplot(aes(x = timeindex, y = value, color = variable)) + 
#   #   geom_line() + 
#   #   scale_color_locuszoom()
#   
#   
#   ## Combine *****************************************************
#   
#   # oup %>% dim
#   # edit_ddj %>% dim
#   # edit_generic %>% dim
#   
#   
#   
#   # combined_df <- full_join(oup_edit, edit_generic, "timeindex") %>%
#   #   mutate(series = ind)
#   
#   
#   combined_df <- full_join(oup, generic, "timeindex")
#   combined_df <- full_join(combined_df, ar1, c("timeindex", "from", "to", "window")) %>% 
#     mutate(series = ind)
#   
#   # Remove .x and .y from column names
#   colnames(combined_df)[grep("\\.y", colnames(combined_df))] <- grep("\\.y", colnames(combined_df), value = T) %>% gsub("\\.y", "_ar1", .)
#   colnames(combined_df)[grep("\\.x", colnames(combined_df))] <- grep("\\.x", colnames(combined_df), value = T) %>% gsub("\\.x", "", .)
#   
#   
#   
#   
#   return(combined_df)
#   
# }) %>% do.call(rbind, .)
# 
# # compare_ews_df <- rbind(compare_ews_df, xx)
# 
# # xx <- compare_ews_df
# 
# compare_ews_df <- compare_ews_df %>% 
#   mutate(data = "ews")
# 
# saveRDS(object = compare_ews_df, file = "EWS/EWS_ews_results_10")
# 





## Results ************************************************** ####

# Table with times at which the 95% CIs parted for the first time
CI_table <- lapply(1:N_series, function(s) {
  
  ## OUP Lambda and sigma ****************************************
  oup_df <- lapply(c("lambda", "unit_variance", "sigma"), function(par) {
    
    # CSD
    my_ews <- ews_res[["oup"]] %>% 
      filter(series == s) %>%
      select(timeindex,
             mode = paste0(par, "_mode"),
             lower2.5 = paste0(par, "_lower2.5"),
             upper97.5 = paste0(par, "_upper97.5"))
    # NULL
    my_null <- null_res[["oup"]] %>% 
      select(timeindex = timeindex,
             mode = paste0(par, "_mode"),
             lower2.5 = paste0(par, "_lower2.5"),
             upper97.5 = paste0(par, "_upper97.5"),
             series)
    
    
    # Compare CSD and null seriers 95% CIs 
    my_df <- lapply(unique(my_null$series), function(i) {
      
      null_limit <- my_null %>% 
        filter(series == i) %>% 
        pull(lower2.5)
      
      # last index where null model is above
      last_ind <- last(which(null_limit < my_ews$upper97.5))
      last_ind <- ifelse(is.na(last_ind), 1, last_ind)
      
      
      data.frame(null_series = i, ews_series = s,
                 last_time = my_null$timeindex[last_ind],
                 parameter = paste0(par, "_oup"))
      
    }) %>% do.call(rbind, .)
    
    
    my_df
    
  }) %>% do.call(rbind, .)
  
  
  ## AR(1) Lambda ****************************************
  ar1_df <- lapply(c("lambda", "sigma"), function(par) {
    
    
    # CSD
    my_ews <- ews_res[["ar1"]] %>% 
      filter(series == s) %>% 
      select(timeindex,
             mode = paste0(par, "_mode"),
             lower2.5 = paste0(par, "_lower2.5"),
             upper97.5 = paste0(par, "_upper97.5"))
    # NULL
    my_null <- null_res[["ar1"]] %>% 
      select(timeindex,
             mode = paste0(par, "_mode"),
             lower2.5 = paste0(par, "_lower2.5"),
             upper97.5 = paste0(par, "_upper97.5"),
             series)
    
    
    
    my_df <- lapply(unique(my_null$series), function(i) {
      
      null_limit <- my_null %>% 
        filter(series == i) %>% 
        pull(lower2.5)
      
      # last index where null model is above
      last_ind <- last(which(null_limit < my_ews$upper97.5))
      last_ind <- ifelse(is.na(last_ind), 1, last_ind)
      
      
      data.frame(null_series = i, ews_series = s, last_time = my_null$timeindex[last_ind], parameter = paste0(par, "_ar1"))
      
    }) %>% do.call(rbind, .)
    
    return(my_df)
    
  }) %>% do.call(rbind, .)
  
  
  ## Combine
  rbind(oup_df, ar1_df)
  
  
}) %>% do.call(rbind, .)

## Kendall's Tau of all variables. For the Bayesian posteriors mode is used
tau_res <- lapply(unique(ews_res[["oup"]]$series), function(s) {
  
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
    dplyr::select(lambda_mode, sigma_mode, unit_variance_mode, timeindex) %>% 
    set_colnames(c(paste0("OUP_", c("lambda", "sigma", "unit_variance")), "timeindex"))
  
  
  cor_res <- lapply(list(generic_df, ddj_df, ar1_df, oup_df), function(df) {
    
    timeindex <- df$timeindex
    df <- df %>% dplyr::select(-timeindex)
    
    lapply(colnames(df), function(i) {
      
      cor(df[, i], timeindex, method = "kendall", use = "complete.obs")
      
    }) %>% set_names(colnames(df)) %>% unlist
    
    
  }) %>% unlist %>% as.data.frame() %>% t %>% cbind(series = s)
  
  
  return(cor_res)   
  
  
}) %>% do.call(rbind, .) %>% as.data.frame()


## Plot ***************************************************** ####

pos_tau_res <- lapply(colnames(tau_res[, 1:(ncol(tau_res) - 1)]), function(i) {
  
  x <- tau_res[, i]
  m <- mean(x)
  
  
  if(m < 0){
    x <- -x
  }
  
  
  return(x)
}) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>% 
  set_colnames(colnames(tau_res[, 1:(ncol(tau_res) - 1)])) %>% 
  cbind(., series = tau_res$series)

## Kendall's tau
pos_tau_res %>% 
  melt(id.vars = "series") %>% 
  mutate(variable = factor(variable, levels = sort(colMeans(pos_tau_res %>% select(-series)) %>% abs) %>% names)) %>%
  ggplot() +
  geom_boxplot(aes(x = (variable), y = value), position = "dodge") + 
  theme(legend.position="top")




## CI Results
# Box plot: lambda and variance
p1 <- CI_table %>%
  ggplot(aes(x = factor(parameter, levels = c("lambda_oup", "lambda_ar1", "unit_variance_oup", "sigma_ar1", "sigma_oup")), y = last_time)) +
  geom_boxplot() +
  labs(x = "Parameter", y = "Time point", subtitle = "Distributions of last time indices after which the 95% CIs no longer intersect") + 
  geom_signif(comparisons = list(c("lambda_oup", "lambda_ar1"), c("unit_variance_oup", "sigma_ar1")), )


i <- 3

# Example: lambda
p2 <- ggplot() + 
  geom_line(data = null_res[["oup"]] %>% filter(series == i),
            aes(x = timeindex, y = lambda_mode),
            size = .5) + 
  geom_ribbon(data = null_res[["oup"]] %>% filter(series == i),
              aes(x = timeindex, ymin = lambda_lower2.5, ymax = lambda_upper97.5),
              alpha = 0.1) + 
  geom_line(data = ews_res[["oup"]] %>% filter(series == i),
            aes(x = timeindex, y = lambda_mode),
            size = .5,
            color = "blue") + 
  geom_ribbon(data = ews_res[["oup"]] %>% filter(series == i),
              aes(x = timeindex, ymin = lambda_lower2.5, ymax = lambda_upper97.5),
              alpha = 0.2, fill = "royalblue") + 
  geom_vline(xintercept = CI_table %>%
               filter(null_series == i, ews_series == i, parameter == "lambda_oup") %>%
               pull(last_time),
             linetype = "dashed") + 
  labs(x = "Time", y = "Value",
       subtitle = "OUP: lambda parameter; grey = null model; blue = sliding window")


# Example: AR1 lambda
p3 <- ggplot() + 
  geom_line(data = null_res[["ar1"]] %>% filter(series == i),
            aes(x = timeindex, y = lambda_mode),
            size = .5) +
  geom_ribbon(data = null_res[["ar1"]] %>% filter(series == i),
              aes(x = timeindex, ymin = lambda_lower2.5, ymax = lambda_upper97.5),
              alpha = 0.1) +
  geom_line(data = ews_res[["ar1"]] %>% filter(series == i),
            aes(x = timeindex, y = lambda_mode),
            size = .5, color = "red") +
  geom_ribbon(data = ews_res[["ar1"]] %>% filter(series == i),
              aes(x = timeindex, ymin = lambda_lower2.5, ymax = lambda_upper97.5),
              alpha = 0.2, fill = "salmon") + 
  geom_vline(xintercept = CI_table %>% filter(null_series == i, ews_series == i, parameter == "lambda_ar1") %>% pull(last_time),
             linetype = "dashed") + 
  labs(x = "Time", y = "Value", subtitle = "AR(1): lambda parameter; grey = null model; red = sliding window")

# OUP sigma
p4 <- ggplot() + 
  geom_line(data = null_res[["oup"]] %>% filter(series == i), aes(x = timeindex, y = sigma_mode), size = .5) + 
  geom_ribbon(data = null_res[["oup"]] %>% filter(series == i), aes(x = timeindex, ymin = sigma_lower2.5, ymax = sigma_upper97.5), alpha = 0.1) + 
  geom_line(data = ews_res[["oup"]] %>% filter(series == i), aes(x = timeindex, y = sigma_mode), size = .5, color = "blue") + 
  geom_ribbon(data = ews_res[["oup"]] %>% filter(series == i), aes(x = timeindex, ymin = sigma_lower2.5, ymax = sigma_upper97.5), alpha = 0.2, fill = "royalblue") + 
  geom_vline(xintercept = CI_table %>% filter(null_series == i, ews_series == i, parameter == "sigma_oup") %>% pull(last_time), linetype = "dashed") + 
  labs(x = "Time", y = "Value", subtitle = "OUP: sigma parameter; grey = null model; blue = sliding window")


# OUP unit variance
p4 <- ggplot() + 
  geom_line(data = null_res[["oup"]] %>% filter(series == i), aes(x = timeindex, y = unit_variance_mode), size = .5) + 
  geom_ribbon(data = null_res[["oup"]] %>% filter(series == i), aes(x = timeindex, ymin = unit_variance_lower2.5, ymax = unit_variance_upper97.5), alpha = 0.1) + 
  geom_line(data = ews_res[["oup"]] %>% filter(series == i), aes(x = timeindex, y = unit_variance_mode), size = .5, color = "blue") + 
  geom_ribbon(data = ews_res[["oup"]] %>% filter(series == i), aes(x = timeindex, ymin = unit_variance_lower2.5, ymax = unit_variance_upper97.5), alpha = 0.2, fill = "royalblue") + 
  geom_vline(xintercept = CI_table %>% filter(null_series == i, ews_series == i, parameter == "unit_variance_oup") %>% pull(last_time), linetype = "dashed") + 
  labs(x = "Time", y = "Value", subtitle = "OUP: unit variance; grey = null model; blue = sliding window")


# AR1 unit variance
p5 <- ggplot() + 
  geom_line(data = null_res[["ar1"]] %>% filter(series == i), aes(x = timeindex, y = sigma_mode), size = .5) + 
  geom_ribbon(data = null_res[["ar1"]] %>% filter(series == i), aes(x = timeindex, ymin = sigma_lower2.5, ymax = sigma_upper97.5), alpha = 0.1) + 
  geom_line(data = ews_res[["ar1"]] %>% filter(series == i), aes(x = timeindex, y = sigma_mode), size = .5, color = "red") + 
  geom_ribbon(data = ews_res[["ar1"]] %>% filter(series == i), aes(x = timeindex, ymin = sigma_lower2.5, ymax = sigma_upper97.5), alpha = 0.2, fill = "salmon") + 
  geom_vline(xintercept = CI_table %>% filter(null_series == i, ews_series == i, parameter == "sigma_ar1") %>% pull(last_time), linetype = "dashed") + 
  labs(x = "Time", y = "Value", subtitle = "AR(1): sigma parameter; grey = null model; red = sliding window")




plot_grid(p2, p3, p4, p5, ncol = 1)




## Results ************************************************** ####


ggplot() + 
  geom_line(data = data.frame(ews_set[, c("V2", "x")], resids = resids) %>%
              melt(id.vars = "x") %>% mutate(residuals = ifelse(variable == "V2", FALSE, TRUE)), aes(x = x, y = value)) + 
  geom_smooth(data = compare_ews_df[, c("iqr", "timeindex", "series", "residuals")] %>% 
                # filter(series >= 89) %>% 
                melt(id.vars = c("timeindex", "series", "residuals")), 
              aes(x = timeindex, y = value), color = "black") + 
  # geom_line(data = compare_ews_df[, c("ar1_std", "sd_std", "sk_std", "kurt_std", "cv_std", "densratio_std", "timeindex", "series")] %>% 
  # filter(series >= 89) %>% 
  #   melt(id.vars = c("timeindex", "series")), 
  # aes(x = timeindex, y = value, color = variable)) + 
  # facet_wrap(~series, scales = "free") +
  geom_errorbar(data = compare_ews_df, aes(x = timeindex, ymin = pred_mean - pred_sd, ymax = pred_mean + pred_sd), color = "royalblue", alpha = 0.2) + 
  labs(subtitle = "Black = iqr") + 
  facet_wrap(~residuals) +
  scale_colour_locuszoom()




ggplot() + 
  geom_line(data = compare_ews_df[, c("likelihood_edit", "timeindex", "series")] %>% 
              # filter(series >= 89) %>% 
              melt(id.vars = c("timeindex", "series")), 
            aes(x = timeindex, y = value)) +
  geom_line(data = compare_ews_df[, c("lambda_edit", "mu_edit", "sigma_edit", "resilience_edit", "timeindex", "series")] %>% 
              # filter(series >= 89) %>% 
              melt(id.vars = c("timeindex", "series")), 
            aes(x = timeindex, y = value, color = variable)) + 
  # facet_wrap(~series, scales = "free") +
  labs(subtitle = "Black = likelihood") + 
  scale_colour_locuszoom()


#### Kendall Tau
# tau_res <- lapply(c(TRUE, FALSE), function(i) {
#   
#   compare_ews_df %>% 
#     filter(residuals == i) %>% 
#     dplyr::select(iqr, lambda, sigma, resilience, ar1, sd, sk, kurt, cv, densratio, likelihood_ar1, iqr_ar1, lambda_ar1, mu_ar1, sigma_ar1, resilience_ar1) %>% 
#     apply(X = ., MARGIN = 2, FUN = function(x) cor(x, y = compare_ews_df %>% 
#                                                      filter(residuals == TRUE) %>% pull(window), method = "kendall", use = "complete.obs")) %>% 
#     abs() %>% 
#     sort()
#   
# }) %>% 
#   do.call(cbind, .) %>% 
#   as.data.frame() %>% 
#   set_colnames(c("TRUE", "FALSE")) %>% 
#   rownames_to_column(var = "parameter") %>% 
#   melt(id.vars = "parameter") %>% 
#   mutate(residuals = as.logical(variable)) %>% 
#   dplyr::select(-variable)




#####

pars <- c("window", "from", "to", "timeindex", "residuals", "series")

generic_pars <- c("ar1", "sd", "sk", "kurt", "cv", "densratio", "acf1")
oup_pars <- c("lambda", "mu", "sigma", "resilience", "likelihood", "iqr")
ar1_pars <- c("likelihood_ar1", "iqr_ar1", "lambda_ar1", "mu_ar1", "sigma_ar1", "resilience_ar1")

ac_pars <- c("ar1", "lambda", "lambda_ar1")
sd_pars <- c("sd", "sigma", "sigma_ar1")
resilience_pars <- c("resilience", "resilience_ar1")

# compare_ews_df %>% 
#   dplyr::select(iqr, lambda, sigma, resilience, ar1, sd, sk, kurt, cv, densratio, likelihood_ar1, likelihood, iqr_ar1, lambda_ar1, mu_ar1, sigma_ar1, resilience_ar1, timeindex, residuals) %>% 
#   melt(id.vars = c("timeindex", "residuals")) %>% 
#   ggplot(aes(x = timeindex, y = value, color = residuals)) + 
#   geom_smooth() + 
#   facet_wrap(c("variable"), scales = "free")


trans_first <- function(df, vars, group) {
  
  my_df <- df %>% as.data.frame
  
  for(g in unique(my_df[, group])) {
    for(v in vars) {
      
      my_df[my_df[, group] == g, v] <- my_df[my_df[, group] == g, v] - my_df[my_df[, group] == g, ][1, v]
      
    }
  }
  
  
  my_df
}

# Set starting point to same in all

ac_trans_compare_ews_df <- lapply(unique(compare_ews_df$series), function(s) {
  
  my_df <- compare_ews_df %>% 
    dplyr::select(one_of(c(ac_pars, "timeindex", "residuals", "series"))) %>% 
    mutate(lambda = -lambda, lambda_ar1 = -lambda_ar1) %>% 
    filter(series == s)
  
  for(g in unique(my_df[, "residuals"])) {
    for(v in ac_pars) {
      
      my_df[my_df[, "residuals"] == g, v] <- my_df[my_df[, "residuals"] == g, v] - my_df[my_df[, "residuals"] == g, ][1, v]
      
    }
  }
  
  
  my_df
  
}) %>% do.call(rbind, .)
ac_trans_compare_ews_df %>% 
  filter(series < 5) %>%
  melt(id.vars = c("timeindex", "residuals", "series")) %>% 
  ggplot(aes(x = timeindex, y = value, color = variable)) + 
  geom_line() + 
  facet_wrap(c("series", "residuals"), labeller = labeller(residuals = label_both, series = label_both), ncol = 2, scales = "free") + 
  scale_color_manual(name = "Variable",
                     values=c("red", "blue", "green"), 
                     breaks=c("ar1", "lambda", "lambda_ar1"),
                     labels=c("AR(1) ac", "OUP -lambda", "Bayes AR(1) ac"))



sd_trans_compare_ews_df <- lapply(unique(compare_ews_df$series), function(s) {
  
  my_df <- compare_ews_df %>% 
    dplyr::select(one_of(c(sd_pars, "timeindex", "residuals", "series"))) %>% 
    filter(series == s)
  
  for(g in unique(my_df[, "residuals"])) {
    for(v in sd_pars) {
      
      my_df[my_df[, "residuals"] == g, v] <- my_df[my_df[, "residuals"] == g, v] - my_df[my_df[, "residuals"] == g, ][1, v]
      
    }
  }
  
  
  my_df
  
}) %>% do.call(rbind, .)

sd_trans_compare_ews_df %>% 
  filter(series < 5) %>%
  melt(id.vars = c("timeindex", "residuals", "series")) %>% 
  ggplot(aes(x = timeindex, y = value, color = variable)) + 
  geom_line() + 
  facet_wrap(c("series", "residuals"), labeller = labeller(residuals = label_both, series = label_both), ncol = 2, scales = "free") + 
  scale_color_manual(name = "Variable",
                     values=c("red", "blue", "green"), 
                     breaks=c("sd", "sigma", "sigma_ar1"),
                     labels=c("AR(1) sd", "OUP sigma", "Bayes AR(1) sd"))


#####

compare_ews_df2 <- compare_ews_df %>% 
  mutate(oup_sd = (sigma^2)/(2*lambda)*(1-exp(-2*lambda)))


sd_trans_compare_ews_df <- lapply(unique(compare_ews_df2$series), function(s) {
  
  my_df <- compare_ews_df2 %>% 
    dplyr::select(one_of(c(sd_pars, "oup_sd", "timeindex", "residuals", "series"))) %>% 
    filter(series == s)
  
  for(g in unique(my_df[, "residuals"])) {
    for(v in c(sd_pars, "oup_sd")) {
      
      my_df[my_df[, "residuals"] == g, v] <- my_df[my_df[, "residuals"] == g, v] - my_df[my_df[, "residuals"] == g, ][1, v]
      
    }
  }
  
  
  my_df
  
}) %>% do.call(rbind, .)

sd_trans_compare_ews_df %>% 
  filter(series < 5) %>%
  melt(id.vars = c("timeindex", "residuals", "series")) %>% 
  ggplot(aes(x = timeindex, y = value, color = variable)) + 
  geom_line() + 
  facet_wrap(c("series", "residuals"), labeller = labeller(residuals = label_both, series = label_both), ncol = 2, scales = "free") + 
  scale_color_manual(name = "Variable",
                     values=c("red", "blue", "green"), 
                     breaks=c("sd", "sigma", "sigma_ar1"),
                     labels=c("AR(1) sd", "OUP sigma", "Bayes AR(1) sd"))


#####

xx <- compare_ews_df[, c("timeindex", "ar1", "lambda_ar1", "lambda", "series", "residuals")]
# xx <- compare_ews_df[, c("timeindex", "sd", "sigma", "sigma_ar1", "series", "residuals")]
# xx$ar1 <- (xx$ar1/max(xx$ar1) - xx$ar1[1])
# xx$iqr <- xx$iqr/max(xx$iqr) - xx$iqr[1]

xx %>% 
  melt(id.vars = c("timeindex", "series", "residuals")) %>% 
  ggplot(aes(x = timeindex, y = value, color = variable)) + 
  geom_line(alpha = .5) + 
  geom_smooth() +
  facet_wrap(c("series", "residuals"), scales = "free")



