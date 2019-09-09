# Functions

## EWS  ********************** ####

## Hierarchical *****
ar1_hierarchical_ews <- function(series_set, window_prop = .5, iter = 1000, chains = 1) {
  
  my_data <- series_set
  
  # ar1_hierarchical <- stan_model("stan_models/ar1_hierarchical.stan")
  
  # Proportion of time series used for a single window
  window_prop <- window_prop
  window_length <- (window_prop*nrow(my_data) - 1) %>% round
  N_window <- nrow(my_data) - window_length + 1
  
  # windows_index <- which((1:N_window)%% 100 == 1)
  
  # Results
  my_res <- lapply(1:N_window, function(j) {
    
    dat <- my_data[j:(j + window_length - 1), ]
    
    
    ar1_samples <- sampling(ar1_hierarchical, list(N_times = nrow(dat),
                                                   N_series = ncol(dat),
                                                   time = 1:nrow(dat),
                                                   Y = t(dat)),
                            chains = chains, iter = iter)
    
    
    
    # Paramters
    ar1_pars <- c("lambda", "mu", "sigma", "c", "phi", "stat_variance")
    coefs <- lapply(ar1_pars, function(p) {
      
      my_df <- get_stan_results(ar1_samples, paste0("^", p, "\\["), regex = TRUE) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% 
        set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      mutate(series = 1:ncol(my_data)) %>% 
      set_rownames(NULL)
    
    
    # Hyperparameters
    ar1_hyperpars <- paste0(rep(c("lambda", "sigma", "mu"), each = 2), c("_a", "_b"))
    hyper_coefs <- lapply(ar1_hyperpars, function(p) {
      
      my_df <- get_stan_results(ar1_samples, p, regex = TRUE) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% 
        set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      mutate(window = j) %>% 
      set_rownames(NULL)
    
    x <- coefs %>% mutate(timeindex = (j + window_length - 1))
    
    list("coefs" = x, "hyper_coefs" = hyper_coefs)
    
  }) 

  # Rbind the parameter and hyperparamter tables
  par_res <- lapply(1:length(my_res), function(i) {
    
    my_res[[i]][["coefs"]]
    
  }) %>% do.call(rbind, .)
  hyper_par_res <- lapply(1:length(my_res), function(i) {
    
    my_res[[i]][["hyper_coefs"]]
    
  }) %>% do.call(rbind, .)
  
  
  
  list("parameters" = par_res, "hyperparameters" = hyper_par_res)
  
}
oup_hierarchical_ews <- function(series_set, window_prop = .5, iter = 1000, chains = 1) {
  
  my_data <- series_set
  
  # oup_hierarchical <- stan_model("stan_models/oup_hierarchical_transition.stan")
  
  # Proportion of time series used for a single window
  window_prop <- window_prop
  window_length <- (window_prop*nrow(my_data) - 1) %>% round
  N_window <- nrow(my_data) - window_length + 1
  
  # windows_index <- which((1:N_window)%% 100 == 1)
  
  # Results
  my_res <- lapply(1:N_window, function(j) {
    
    dat <- my_data[j:(j + window_length - 1), ]
    
    
    oup_samples <- sampling(oup_hierarchical, list(N_times = nrow(dat),
                                                   N_series = ncol(dat),
                                                   time = 1:nrow(dat),
                                                   Y = t(dat)),
                            chains = chains, iter = iter)
    
    
    
    # Paramters
    oup_pars <- c("length_scale", "sigma", "stat_var")
    coefs <- lapply(oup_pars, function(p) {
      
      my_df <- get_stan_results(oup_samples, paste0("^", p, "\\["), regex = TRUE) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% 
        set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      mutate(series = 1:ncol(my_data)) %>% 
      set_rownames(NULL)
    
    
    # Hyperparameters
    oup_hyperpars <- paste0(rep(c("length_scale", "stat_var"), each = 2), c("_a", "_b"))
    hyper_coefs <- lapply(oup_hyperpars, function(p) {
      
      my_df <- get_stan_results(oup_samples, p, regex = TRUE) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% 
        set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      mutate(window = j) %>% 
      set_rownames(NULL)
    
    x <- coefs %>% mutate(timeindex = (j + window_length - 1))
    
    list("coefs" = x, "hyper_coefs" = hyper_coefs)
    
  }) 
  
  # Rbind the parameter and hyperparamter tables
  par_res <- lapply(1:length(my_res), function(i) {
    
    my_res[[i]][["coefs"]]
    
  }) %>% do.call(rbind, .)
  hyper_par_res <- lapply(1:length(my_res), function(i) {
    
    my_res[[i]][["hyper_coefs"]]
    
  }) %>% do.call(rbind, .)
  
  
  
  list("parameters" = par_res, "hyperparameters" = hyper_par_res)
  
  
}
se_hierarchical_ews <- function(series_set, window_prop = .5, iter = 1000, chains = 1) {
  
  my_data <- series_set
  
  # se_hierarchical <- stan_model("stan_models/se_hierarchical_covariance.stan")
  
  # Proportion of time series used for a single window
  window_prop <- window_prop
  window_length <- (window_prop*nrow(my_data) - 1) %>% round
  N_window <- nrow(my_data) - window_length + 1
  
  # windows_index <- which((1:N_window)%% 100 == 1)
  
  # Results
  my_res <- lapply(1:N_window, function(j) {
    
    dat <- my_data[j:(j + window_length - 1), ]

    oup_samples <- sampling(se_hierarchical, list(N = nrow(dat),
                                                   N_series = ncol(dat),
                                                   x = 1:nrow(dat),
                                                   y = dat),
                            chains = chains, iter = iter)
    
    
    
    # Paramters
    oup_pars <- c("length_scale", "sigma", "stat_var")
    coefs <- lapply(oup_pars, function(p) {
      
      my_df <- get_stan_results(oup_samples, paste0("^", p, "\\["), regex = TRUE) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% 
        set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      mutate(series = 1:ncol(my_data)) %>% 
      set_rownames(NULL)
    
    
    # Hyperparameters
    oup_hyperpars <- paste0(rep(c("length_scale", "stat_var"), each = 2), c("_a", "_b"))
    hyper_coefs <- lapply(oup_hyperpars, function(p) {
      
      my_df <- get_stan_results(oup_samples, p, regex = TRUE) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% 
        set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      mutate(window = j) %>% 
      set_rownames(NULL)
    
    x <- coefs %>% mutate(timeindex = (j + window_length - 1))
    
    list("coefs" = x, "hyper_coefs" = hyper_coefs)
    
  }) 
  
  # Rbind the parameter and hyperparamter tables
  par_res <- lapply(1:length(my_res), function(i) {
    
    my_res[[i]][["coefs"]]
    
  }) %>% do.call(rbind, .)
  hyper_par_res <- lapply(1:length(my_res), function(i) {
    
    my_res[[i]][["hyper_coefs"]]
    
  }) %>% do.call(rbind, .)
  
  
  
  list("parameters" = par_res, "hyperparameters" = hyper_par_res)
  
  
}





## GP EWS

GP_ews <- function(x, y, window_length, iter = 1000, chains = 1, kernel = "OUP") {
  
  ## GET MODEL
  if(kernel == "OUP") {
    model <- oup_model
  } else if(kernel == "se") {
    model <- se_model
  } else if(kernel == "matern_1.5") {
    model <- matern_1.5_model
  } else if(kernel == "matern_2.5") {
    model <- matern_2.5_model
  }
  
  ## Number of windows needed
  N_windows <- which(x <= (x[length(x)] - window_length)) %>% length
  
  
  ## Stan
  my_res <- lapply(1:N_window, function(j) {
    
    indx <- which((x[j] + window_length) >= x)
    
    stan_data <- list(N = length(indx), 
                      N_pred = 1, 
                      y = y[indx], 
                      x = x[indx],
                      x_pred = as.array(j + window_length))
    
    
    samples <- sampling(model,
                            stan_data,
                            chains = chains, iter = iter)
    
    
    
    pars <- c("rho", "sigma", "lp")
    coefs <- lapply(pars, function(p) {
      
      my_df <- get_stan_results(samples, p) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      set_rownames(NULL)
    
   
    # Finalize results
    vec <- c(window = j, from = j, to = last(x[indx]), (coefs %>% as.matrix)[1, ])
    
    vec
    
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  
  my_res
  
  
  
}
GP_ews_set <- function(data, iter = 500, chains = 1, window_length, kernel = "OUP") {
  
  ## GET MODEL
  if(kernel == "OUP") {
    model <- oup_model
  } else if(kernel == "se") {
    model <- se_model
  } else if(kernel == "matern_1.5") {
    model <- matern_1.5_model
  } else if(kernel == "matern_2.5") {
    model <- matern_2.5_model
  }
  
  
  
  # Separate observations and time
  obs <- data[, 1:(ncol(data) - 1)] %>% data.frame()
  obs_times <- data[, ncol(data)]
  
  
  ## OUP Bayes *******************************
  print("Bayesian OUP")
  res <- lapply(1:ncol(obs), function(i) {
    
    print(paste0("OUP: ", i, " out of ", ncol(obs)))
    
    series <- obs[, i]
    
    gp <- GP_ews(x = obs_times, y = series, iter = iter, chains = chains, window_length, kernel = kernel)
    
    gp <- gp %>%
      mutate(timeindex = gp$to, series = i)
    
    return(gp)
    
  }) %>%
    do.call(rbind, .)
  
  res
  
}

# AR1 ews with missing values
ar1_missing_ews_set <- function(data, iter = 500, chains = 1, window_length) {
  
  mydata <- data %>% 
    arrange(x) %>% 
    mutate(x = round(x))
  
  na_index <- (min(mydata$x):max(mydata$x))[which(!(min(mydata$x):max(mydata$x) %in% mydata$x))]
  obs_times <- mydata$x
  
  # Separate observations and time
  obs <- mydata[, 1:(ncol(mydata) - 1)] %>% data.frame()
  
  ## AR(1) Bayes *****************************
  print("Bayesian AR(1)")
  ar1_res <- lapply(1:ncol(obs), function(i) {
    
    print(paste0("AR(1): ", i, " out of ", ncol(obs)))
    
    series <- obs[, i]
    
    ar1 <- ar1_missing_ews(series = series, na_index, obs_times,
                           iter = iter, chains = chains, window_length =  window_length)

    ar1 <- ar1 %>%
      mutate(timeindex = ar1$to, series = i)
    
    return(ar1)
    
  }) %>%
    do.call(rbind, .)
  
  return(ar1_res)
  
}
ar1_missing_ews <- function(series, na_index, obs_index, window_length, iter = 1000, chains = 1) {
  
  my_data <- series
  x <- obs_index
  
  ## Number of windows needed
  N_windows <- which(x <= (x[length(x)] - window_length)) %>% length
  
  # ar1_missing_model <- stan_model("stan_models/ar1_missing.stan")
  
  my_res <- lapply(1:N_windows, function(j) {
    
    present <- which(j:(j + window_length - 1) %in% obs_index)
    absent <- which(j:(j + window_length - 1) %in% na_index)
    
    dat <- my_data[present]
    
    ar1_samples <- sampling(ar1_missing_model,
                            list(N_obs = length(present),
                                 N_mis = length(absent),
                                 ii_obs = present,
                                 ii_mis = absent,
                                 Y_obs = dat),
                            chains = chains, iter = iter)
    
    
    
    oup_pars <- c("lambda", "sigma", "c", "phi", "stat_variance")
    coefs <- lapply(oup_pars, function(p) {
      
      my_df <- get_stan_results(ar1_samples, p, regex = FALSE) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      set_rownames(NULL)
    
    
    
    
    
    # Next time point likelihood
    # last <- my_data[j + window_length - 1]
    # present <- my_data[j + window_length]
    
    # LL <- dnorm(present,
    #             mean = coefs["c"] + coefs["phi"]*present,
    #             sd = coefs["sigma"]^2) 
    # Q <- abs(0.5 - pnorm(present, mean = coefs["c"] + coefs["phi"]*present, sd = coefs["sigma"]))*2
    
    # x <- c(likelihood_ar1 = LL, iqr_ar1 = Q, window = j, from = j, to = (j + window_length - 1), pred_mean_ar1 = unname(next_m), pred_sd_ar1 = unname(next_sd), coefs %>% set_names(paste0(oup_pars, "_ar1")))
    
    x <- c(window = j, from = j, to = (j + window_length - 1), (coefs %>% as.matrix())[1, ])
    
    x
    
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  
  my_res
  
  
  
}







# Bayesian OUP and AR(1)
oup_ews <- function(x, y, window_length, iter = 1000, chains = 1) {
  
  # oup_model <- stan_model("stan_models/oup_fitter.stan")
  
  # Proportion of time series used for a single window
  N_window <- length(y) - window_length + 1
  
  my_res <- lapply(1:N_window, function(j) {
    
    win_y <- y[j:(j + window_length - 1)]
    
    stan_data <- list(N = length(win_y), 
                      N_pred = 1, 
                      y = win_y, 
                      x = x[j:(j + window_length - 1)],
                      x_pred = as.array(j + window_length))
    
    
    oup_samples <- sampling(oup_model,
                            stan_data,
                            chains = chains, iter = iter)
    
    
    
    oup_pars <- c("rho", "sigma", "lp")
    coefs <- lapply(oup_pars, function(p) {
      
      my_df <- get_stan_results(oup_samples, p) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      set_rownames(NULL)
    
    
    
    # Next time point likelihood
    # next_m <- coefs["mu"] - (coefs["mu"] - my_data[j + window_length - 1])*exp(-coefs["lambda"])
    # next_sd <- ((coefs["sigma"]^2)/(2*coefs["lambda"]))*(1-exp(-2*coefs["lambda"]))
    # 
    # present <- my_data[j + window_length]
    # 
    # LL <- dnorm(present,
    #             mean = next_m,
    #             sd = next_sd) 
    # Q <- abs(0.5 - pnorm(present, mean = next_m, sd = next_sd))*2
    
    # x <- c(likelihood = LL, iqr = Q, window = j, from = j, to = (j + window_length - 1), pred_mean = unname(next_m), pred_sd = unname(next_sd), coefs)
    
    vec <- c(window = j, from = j, to = (j + window_length - 1), (coefs %>% as.matrix)[1, ])
    
    vec
    
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  
  my_res
  
  
  
}
ar1_ews <- function(x, y, window_prop = .5, iter = 1000, chains = 1) {
  
  
  
  # ar1_model <- stan_model("stan_models/ar1.stan")
  
  # Proportion of time series used for a single window
  # window_prop <- .5
  window_length <- (window_prop*length(y) - 1) %>% round
  N_window <- length(y) - window_length + 1
  
  # windows_index <- which((1:N_window)%% 100 == 1)
  
  my_res <- lapply(1:N_window, function(j) {
    
    dat <- y[j:(j + window_length - 1)]
    
    ar1_samples <- sampling(ar1_model, list(T = length(dat),
                                            time = x[j:(j + window_length - 1)],
                                            Y = dat),
                            chains = chains, iter = iter)
    
    
    
    oup_pars <- c("lambda", "sigma", "c", "phi", "stat_variance")
    coefs <- lapply(oup_pars, function(p) {
      
      my_df <- get_stan_results(ar1_samples, p, regex = FALSE) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      set_rownames(NULL)
    
    
    
    
    
    # Next time point likelihood
    # last <- my_data[j + window_length - 1]
    # present <- my_data[j + window_length]
    
    # LL <- dnorm(present,
    #             mean = coefs["c"] + coefs["phi"]*present,
    #             sd = coefs["sigma"]^2) 
    # Q <- abs(0.5 - pnorm(present, mean = coefs["c"] + coefs["phi"]*present, sd = coefs["sigma"]))*2
    
    # x <- c(likelihood_ar1 = LL, iqr_ar1 = Q, window = j, from = j, to = (j + window_length - 1), pred_mean_ar1 = unname(next_m), pred_sd_ar1 = unname(next_sd), coefs %>% set_names(paste0(oup_pars, "_ar1")))
    
    x <- c(window = j, from = j, to = (j + window_length - 1), (coefs %>% as.matrix())[1, ])
    
    x
    
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  
  my_res
  
  
  
}
se_ews <- function(x, y, window_prop = .5, iter = 1000, chains = 1) {
  
  # se_model <- stan_model("stan_models/square_exp_fitter.stan")
  
  # Proportion of time series used for a single window
  window_length <- (window_prop*length(y) - 1) %>% round
  N_window <- length(y) - window_length + 1
  
  my_res <- lapply(1:N_window, function(j) {
    
    win_y <- y[j:(j + window_length - 1)]
    
    stan_data <- list(N = length(win_y), 
                      N_pred = 1, 
                      y = win_y, 
                      x = x[j:(j + window_length - 1)],
                      x_pred = as.array(j + window_length))
    
    
    se_samples <- sampling(se_model,
                            stan_data,
                            chains = chains, iter = iter)
    
    
    
    se_pars <- c("rho", "sigma", "lp")
    coefs <- lapply(se_pars, function(p) {
      
      my_df <- get_stan_results(se_samples, p) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      set_rownames(NULL)
    
    
    
    # Next time point likelihood
    # next_m <- coefs["mu"] - (coefs["mu"] - my_data[j + window_length - 1])*exp(-coefs["lambda"])
    # next_sd <- ((coefs["sigma"]^2)/(2*coefs["lambda"]))*(1-exp(-2*coefs["lambda"]))
    # 
    # present <- my_data[j + window_length]
    # 
    # LL <- dnorm(present,
    #             mean = next_m,
    #             sd = next_sd) 
    # Q <- abs(0.5 - pnorm(present, mean = next_m, sd = next_sd))*2
    
    # x <- c(likelihood = LL, iqr = Q, window = j, from = j, to = (j + window_length - 1), pred_mean = unname(next_m), pred_sd = unname(next_sd), coefs)
    
    vec <- c(window = j, from = j, to = (j + window_length - 1), (coefs %>% as.matrix)[1, ])
    
    vec
    
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  
  my_res
  
  
  
}
matern_1.5_ews <- function(x, y, window_prop = .5, iter = 1000, chains = 1) {
  
  # matern_1.5_model <- stan_model("stan_models/matern_1.5_fitter.stan")
  
  # Proportion of time series used for a single window
  window_length <- (window_prop*length(y) - 1) %>% round
  N_window <- length(y) - window_length + 1
  
  my_res <- lapply(1:N_window, function(j) {
    
    win_y <- y[j:(j + window_length - 1)]
    
    stan_data <- list(N = length(win_y), 
                      N_pred = 1, 
                      y = win_y, 
                      x = x[j:(j + window_length - 1)],
                      x_pred = as.array(j + window_length))
    
    
    matern_1.5_samples <- sampling(matern_1.5_model,
                           stan_data,
                           chains = chains, iter = iter)
    
    
    
    matern_1.5_pars <- c("rho", "sigma", "lp")
    coefs <- lapply(matern_1.5_pars, function(p) {
      
      my_df <- get_stan_results(matern_1.5_samples, p) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      set_rownames(NULL)
    
    
    
    # Next time point likelihood
    # next_m <- coefs["mu"] - (coefs["mu"] - my_data[j + window_length - 1])*exp(-coefs["lambda"])
    # next_sd <- ((coefs["sigma"]^2)/(2*coefs["lambda"]))*(1-exp(-2*coefs["lambda"]))
    # 
    # present <- my_data[j + window_length]
    # 
    # LL <- dnorm(present,
    #             mean = next_m,
    #             sd = next_sd) 
    # Q <- abs(0.5 - pnorm(present, mean = next_m, sd = next_sd))*2
    
    # x <- c(likelihood = LL, iqr = Q, window = j, from = j, to = (j + window_length - 1), pred_mean = unname(next_m), pred_sd = unname(next_sd), coefs)
    
    vec <- c(window = j, from = j, to = (j + window_length - 1), (coefs %>% as.matrix)[1, ])
    
    vec
    
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  
  my_res
  
  
  
}
matern_2.5_ews <- function(x, y, window_prop = .5, iter = 1000, chains = 1) {
  
  # matern_2.5_model <- stan_model("stan_models/matern_2.5_fitter.stan")
  
  # Proportion of time series used for a single window
  window_length <- (window_prop*length(y) - 1) %>% round
  N_window <- length(y) - window_length + 1
  
  my_res <- lapply(1:N_window, function(j) {
    
    win_y <- y[j:(j + window_length - 1)]
    
    stan_data <- list(N = length(win_y), 
                      N_pred = 1, 
                      y = win_y, 
                      x = x[j:(j + window_length - 1)],
                      x_pred = as.array(j + window_length))
    
    
    matern_2.5_samples <- sampling(matern_2.5_model,
                                   stan_data,
                                   chains = chains, iter = iter)
    
    
    
    matern_2.5_pars <- c("rho", "sigma", "lp")
    coefs <- lapply(matern_2.5_pars, function(p) {
      
      my_df <- get_stan_results(matern_2.5_samples, p) %>% 
        select(mode, lower2.5, upper97.5)
      
      my_df %>% set_colnames(paste0(p, "_", colnames(my_df)))
      
    }) %>% do.call(cbind, .) %>% 
      set_rownames(NULL)
    
    
    
    # Next time point likelihood
    # next_m <- coefs["mu"] - (coefs["mu"] - my_data[j + window_length - 1])*exp(-coefs["lambda"])
    # next_sd <- ((coefs["sigma"]^2)/(2*coefs["lambda"]))*(1-exp(-2*coefs["lambda"]))
    # 
    # present <- my_data[j + window_length]
    # 
    # LL <- dnorm(present,
    #             mean = next_m,
    #             sd = next_sd) 
    # Q <- abs(0.5 - pnorm(present, mean = next_m, sd = next_sd))*2
    
    # x <- c(likelihood = LL, iqr = Q, window = j, from = j, to = (j + window_length - 1), pred_mean = unname(next_m), pred_sd = unname(next_sd), coefs)
    
    vec <- c(window = j, from = j, to = (j + window_length - 1), (coefs %>% as.matrix)[1, ])
    
    vec
    
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  
  my_res
  
  
  
}



ar1_ews_set <- function(data, iter = 500, chains = 1, window_prop = .5) {
  
  # Separate observations and time
  obs <- data[, 1:(ncol(data) - 1)]  %>% data.frame()
  obs_times <- data[, ncol(data)]
  

  ## AR(1) Bayes *****************************
  print("Bayesian AR(1)")
  ar1_res <- lapply(1:ncol(obs), function(i) {
    
    print(paste0("AR(1): ", i, " out of ", ncol(obs)))
    
    series <- obs[, i]
    
    ar1 <- ar1_ews(x = obs_times, y = series, iter = iter, chains = chains, window_prop = window_prop)
    
    ar1 <- ar1 %>%
      mutate(timeindex = ar1$to, series = i)
    
    return(ar1)
    
  }) %>%
    do.call(rbind, .)
  
  return(ar1_res)
  
}
oup_ews_set <- function(data, iter = 500, chains = 1, window_prop = .5) {
  
  # Separate observations and time
  obs <- data[, 1:(ncol(data) - 1)] %>% data.frame()
  obs_times <- data[, ncol(data)]
  
  
  ## OUP Bayes *******************************
  print("Bayesian OUP")
  oup_res <- lapply(1:ncol(obs), function(i) {
    
    print(paste0("OUP: ", i, " out of ", ncol(obs)))
    
    series <- obs[, i]
    
    oup <- oup_ews(x = obs_times, y = series, iter = iter, chains = chains, window_prop = window_prop)
    
    oup <- oup %>%
      mutate(timeindex = oup$to, series = i)
    
    return(oup)
    
  }) %>%
    do.call(rbind, .)
  
  oup_res
  
}
se_ews_set <- function(data, iter = 500, chains = 1, window_prop = .5) {
  
  # Separate observations and time
  obs <- data[, 1:(ncol(data) - 1)] %>% data.frame()
  obs_times <- data[, ncol(data)]
  
  
  ## OUP Bayes *******************************
  print("Bayesian GP with SE kernel")
  se_res <- lapply(1:ncol(obs), function(i) {
    
    print(paste0("GP: ", i, " out of ", ncol(obs)))
    
    series <- obs[, i]
    
    se <- se_ews(x = obs_times, y = series, iter = iter, chains = chains, window_prop = window_prop)
    
    se <- se %>%
      mutate(timeindex = se$to, series = i)
    
    return(se)
    
  }) %>%
    do.call(rbind, .)
  
  se_res
  
}
matern_1.5_ews_set <- function(data, iter = 500, chains = 1, window_prop = .5) {
  
  # Separate observations and time
  obs <- data[, 1:(ncol(data) - 1)] %>% data.frame()
  obs_times <- data[, ncol(data)]
  
  
  ## OUP Bayes *******************************
  print("Bayesian GP with Matern 1.5 kernel")
  matern_1.5_res <- lapply(1:ncol(obs), function(i) {
    
    print(paste0("GP: ", i, " out of ", ncol(obs)))
    
    series <- obs[, i]
    
    matern <- matern_1.5_ews(x = obs_times, y = series, iter = iter, chains = chains, window_prop = window_prop)
    
    matern <- matern %>%
      mutate(timeindex = matern$to, series = i)
    
    return(matern)
    
  }) %>%
    do.call(rbind, .)
  
  matern_1.5_res
  
}
matern_2.5_ews_set <- function(data, iter = 500, chains = 1, window_prop = .5) {
  
  # Separate observations and time
  obs <- data[, 1:(ncol(data) - 1)] %>% data.frame()
  obs_times <- data[, ncol(data)]
  
  
  ## OUP Bayes *******************************
  print("Bayesian GP with Matern 2.5 kernel")
  matern_2.5_res <- lapply(1:ncol(obs), function(i) {
    
    print(paste0("GP: ", i, " out of ", ncol(obs)))
    
    series <- obs[, i]
    
    matern <- matern_2.5_ews(x = obs_times, y = series, iter = iter, chains = chains, window_prop = window_prop)
    
    matern <- matern %>%
      mutate(timeindex = matern$to, series = i)
    
    return(matern)
    
  }) %>%
    do.call(rbind, .)
  
  matern_2.5_res
  
}
generic_ews_set <- function(data, window_prop = .5) {
  
  # Separate observations and time
  obs <- data[, 1:(ncol(data) - 1)]
  obs_times <- data[, ncol(data)]
  
  
  ## Generic EWS *****************************
  print("Generic")
  generic_res <- lapply(1:ncol(obs), function(i) {
    
    series <- obs[, i]
    
    generic <- generic_ews(series, winsize = window_prop*100)
    dev.off()
    
    generic_pars <- c("ar1", "sd", "sk", "kurt", "cv", "densratio")
    generic <- generic[, c("timeindex", generic_pars)]
    
    generic <- generic %>% mutate(series = i)
    
    return(generic)
    
  }) %>%
    do.call(rbind, .)
}
ddj_ews_set <- function(data) {
 
   # Separate observations and time
  obs <- data[, 1:(ncol(data) - 1)]
  obs_times <- data[, ncol(data)]
  
  
  ## DDJ *************************************
  print("Drift-Diffusion-Jump")
  ddj_res <- lapply(1:ncol(obs), function(i) {
    
    series <- obs[, i]
    
    ddj <- ddjnonparam_ews(data.frame(obs_times, series))
    
    
    # Take timeindex and measures dependent on time
    ddj <- ddj[names(ddj)[6:10]] %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      set_colnames(c("timeindex", "conditional_var", "total_var", "diffusion", "jump_intensity"))
    
    
    ddj <- ddj %>%
      mutate(series = i)
    
    dev.off()
    
    return(ddj)
    
    
  }) %>% 
    do.call(rbind, .)
  
  ddj_res
  
}






# Combine the above into one
ews_wrapper <- function(data, iter = 500, chains = 1) {
  
  my_data <- data
  
  ar1_res <- ar1_ews_set(my_data, iter, chains)
  oup_res <- oup_ews_set(my_data, iter, chains)
  generic_res <- generic_ews_set(my_data, iter, chains)
  ddj_res <- ddj_ews_set(my_data, iter, chains)
  

  return(list("generic" = generic_res,
              "ddj" = ddj_res, 
              "ar1" = ar1_res,
              "oup" = oup_res))
  
}




## Simulate early warning data

ews_drift <- function(x, r, K, c, h) {
  
  return(r*x*(1 - x/K) - c*x^2/(x^2 + h^2))
  
}

ews_dispersion <- function(x, sigma) {
  
  return(sigma*x)
  
}

D_ews_dispersion <- function(x, sigma) {
  
  sigma
  
}

ews_generator <- function(y0, grid, c, seed = NULL, milstein = FALSE, restrict_pos = FALSE) {
  
  # if(is.null(seed)) {
  #   set.seed(sample(1:10000, 1))
  # } else {
  #   set.seed(seed)
  # }
  # 
  
  y <- y0
  for(i in 2:length(grid)) {
    
    y_prev <- y[i-1]
    delta_t <- grid[i] - grid[i-1]
    
    a <- y_prev + ews_drift(y_prev, r, K, c, h)*delta_t + ews_dispersion(y_prev, sigma)*rnorm(1, 0, sqrt(delta_t))
    
    if(milstein) {
      a <- a + .5*ews_dispersion(y_prev, sigma)*D_ews_dispersion(y_prev, sigma)*(rnorm(1, 0, delta_t) - delta_t)
    }
    
    
    if(restrict_pos) {
      a <- ifelse(a <= 0, 0, a)
    }
    y <- c(y, a)
    
    y
    
  }
  
  return(y)
  
}



## Toy data ****************** ####

# Error function, adds stochasticity
errorer <- function(x) {
  
  sqrt(abs(x))*rnorm(1, 0, 0.05)
  
}

# Hill function
f_hill <- function(y, K, p = 2){
  
  d <- length(y)
  
  f <- sapply(1:d, function(i) {
    
    sapply(1:d, function(k) {
      ifelse(i == k, 1, 
             (K[i, k]^p)/(K[i, k]^p + y[k]^p))
    }) %>% prod()
    
  })
  
  return(f)
} 

# Toy model from "Multistability and the.... "; n species
toy_model <- function(times, y, parms) {
  
  f <- f_hill(y = y, K)
  
  dy <- y*(b*f - k*y)
  
  list(dy)
}

# Plot simulated dynamics
plot_dynamics <- function(data, time_col = "time", facet = FALSE) {
  
  p <- data %>% 
    melt(id.vars = time_col) %>% 
    ggplot(aes(x = get(time_col), y = value, color = variable)) + 
    geom_line() +
    labs(x = "Time", y = "", subtitle = "Toy model") +
    theme(legend.title=element_blank())
  
  if(isTRUE(facet)) {
    p <- p + facet_wrap(~variable)
  }
  
  p
  
}

# Get spefic species data
get_species <- function(old_data, new_data,
                        species, time_col = "time") {
  
  news <- new_data[, grep(paste0(species, collapse = "|"),
                          colnames(new_data))]
  
  if(is.null(old_data)) {
    old_data <- cbind(new_data[, time_col], news)
  } else {
    old_data <- cbind(old_data, news)
  }
  
  old_data <- old_data %>% 
    set_colnames(c("time", paste0("S", 1:(ncol(old_data) - 1))))
  
}

toy_data <- function(n_species = 3, times = seq(from = 0, to = 100, by = 0.1), inits = "random", stochastic = TRUE) {
  
  ## Parameters ********************************
  
  # death rate
  # k <- runif(n_species, min = .85, max = 1.15)
  k <- rep(1, n_species)  
  # growth rate
  b <- runif(n_species, min = .85, max = 1.15)
  
  # inhibition matrix 
  K <- matrix(rnorm(n_species^2, .75, .2) %>% round(2),
              n_species, n_species, byrow = T)
  # K <- rep(0.1, n_species^2) %>% matrix(n_species, n_species, byrow = T)
  K <- K - K*diag(n_species)
  
  
  ## Initial values ***************************
  x <- matrix(NA, length(times), n_species + 1)
  if(inits == "random") {
    x[1, ] <- c(0, sample(1:10/10, n_species, replace = T))
    
  } else {
    x[1, ] <- c(0, inits)
  }
  
  
  ## Generation ********************************
  toy_model <- function(times, y, parms) {
    
    f <- f_hill(y = y, K)
    
    dy <- y*(b*f - k*y)
    
    list(dy)
  }
  
  for(i in 2:length(times)) {
    
    out <- ode(times = times[(i-1):i],
               y = x[i-1, -1],
               func = toy_model,
               parms = NULL)
    
    x_new <- out[2, 2:ncol(out)] 
    
    if(stochastic) {
      x_new <- x_new + sapply(x_new + 0.01, errorer)
    }
    
    x_new <- ifelse(x_new < 0, 0, x_new)
    
    x[i, ] <- c(out[2, 1], x_new)
    
  }
  
  x <- x %>%
    as.data.frame() %>%
    set_colnames(c("time", paste0("S", 1:(ncol(x) - 1))))
  
  return(x)
  
}

## Euler ********************* ####
em_generator <- function(y0, grid, seed = NULL, milstein = FALSE, restrict_pos = FALSE) {
  
  if(is.null(seed)) {
    set.seed(sample(1:10000, 1))
  } else {
    set.seed(seed)
  }
  
  
  y <- y0
  for(i in 2:length(grid)) {
    
    y_prev <- y[i-1]
    delta_t <- grid[i] - grid[i-1]
    
    a <- y_prev + cc_drift(y_prev, r, alpha, beta, lambda, epsilon)*delta_t + cc_dispersion(y_prev)*rnorm(1, 0, sqrt(delta_t))
    
    # if(milstein) {
    #   a <- a + .5*dispersion(y_prev)*D_dispersion(y_prev)*(rnorm(1, 0, delta_t) - delta_t)
    # }
    
    
    if(restrict_pos) {
      a <- ifelse(a <= 0, 0, a)
    }
    y <- c(y, a)
    
    y
    
  }
  
  return(y)
  
}
## Shoji ********************* ####

shoji_cusp_transition_density <- function(x, dt, r, alpha, beta, lambda, epsilon) {
  drift <- r*(alpha + beta*(x - lambda) - (x - lambda)^3)
  
  L <- r*(beta - 3*(x - lambda)^2)
  M <- r*((epsilon/2)*(- 6*(x - lambda)))
  
  Ax <-  x + drift*(exp(L*dt) - 1)/L + M*(exp(L*dt) - 1 - L*dt)/L^2
  B <- sqrt(epsilon)*sqrt((exp(2*L*dt) - 1)/(2*L))
  
  rnorm(1, Ax, B)
  
}

shoji_generator <- function(y0, times, r, alpha, beta, lambda, epsilon, seed = 1) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  y <- y0
  
  for(i in 2:length(times)) {
    dt <- times[i] - times[i-1]
    
    y <- c(y, shoji_cusp_transition_density(y[i-1], dt = dt, 
                                            r = r, alpha = alpha, beta = beta, lambda = lambda, epsilon = epsilon))
  }
  
  y
}

## Cusp ********************** ####
cc_drift <- function(x, r, alpha, beta, lambda, epsilon) {
  
  r*(alpha + beta*(x - lambda) - (x - lambda)^3)
  
}

cc_dispersion <- function(x) {
  # sqrt(epsilon*disp_c*x)
  sqrt(epsilon)
}


cc_density <- function(x, rr, aa, bb, ll, ee) {
  
  unnormalized <- function(x, alpha1 = aa, beta1 = bb, lambda1 = ll, epsilon1 = ee, rr1 = rr) {
    exp( (alpha1*(x - lambda1) + .5*beta1*(x - lambda1)^2 - .25*(x - lambda1)^4 ) / (epsilon1/rr1))
  } 
  
  M <- integrate(unnormalized, -Inf, Inf)$value
  
  
  unnormalized(x)/M
  
}



cc_hsc_density <- function(x, rr, aa, bb, ll, ee, dd) {
  
  unnormalized <- function(x, alpha1 = aa, beta1 = bb, lambda1 = ll, epsilon1 = ee, rr1 = rr) {
    exp( (alpha1*(x - lambda1) + .5*beta1*(x - lambda1)^2 - .25*(x - lambda1)^4 ) / (epsilon1/rr1))
  } 
  
  M <- integrate(unnormalized, -Inf, Inf)$value
  
  
  unnormalized(x)/M
  
}

## Forman ******************** ####
transformator <- function(x, alpha, mix_mean, mix_sd, mix_pro) {
  
  qmixnorm(pnorm(x, 0, alpha),
           mean = mix_mean, 
           sd = mix_sd,
           pro = mix_pro)
  
}
## OUP *********************** ####
oup_generator <- function(N = 100, times = 1:100, y0 = 0, mu = 0, lambda = 0.5, sigma = .25) {
  
  y <- rep(0, N)
  
  y[1] <- y0
  
  
  for(i in 2:N) {
    dt <- (times[i] - times[i-1])
    
    transition_mu <- mu - (mu - y[i-1])*exp(-lambda*dt)
    transition_var <- (sigma^2/(2*lambda))*(1 - exp(-2*lambda*dt))
    
    y[i] <- rnorm(1, transition_mu, sqrt(transition_var))
    
  }
  
  return(y)
}

generate_gp_set <- function(N_series, times, covariance, length_scale, stat_var, error = 0,
                            seed = NULL, modules = FALSE) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  if(isTRUE(modules)) {
    my_seed <- ifelse(is_null(seed), sample.int(.Machine$integer.max, 1), seed)
  } else {
    my_seed <- sample.int(.Machine$integer.max, 1)
  }
  
  # Model
  # sim_data_model <- stan_model("stan_models/generate_gp.stan")
  
  # Set covariance
  if(covariance == "oup") {
    kernel <- 0
  } else if(covariance == "se") {
    kernel <- 1
  }
  
  
  time_grid <- seq(from = times[1], to = times[length(times)], by = 0.1)
  
  # Generate
  res <- lapply(1:N_series, function(i) {
    
    dat_list <- list(T = length(time_grid),
                     time = time_grid,
                     stat_var = stat_var[i],
                     length_scale = length_scale[i],
                     error = error[i],
                     kernel = kernel)
    
    
    draw <- sampling(sim_data_model,
                     iter=1,
                     algorithm='Fixed_param',
                     chains = 1,
                     data = dat_list, 
                     seed = my_seed)
    
    
    samps <- rstan::extract(draw)
    plt_df = with(samps, data.frame(y = y[1,], f = f[1,])) %>% 
      mutate(index = i)
    
    return(plt_df)
    
  }) %>%
    do.call(rbind, .) %>% 
    cbind(x = time_grid)
  
  return(res)
}

## AR(1) ********************* ####

ar1_generator <- function(times, lambda, sigma, mu = 0, epsilon) {
  
  phi <- -lambda
  c <- lambda*mu
  
  x <- rep(0, length(times))
  
  x[1] <- rnorm(1, c/(1 - phi), sigma/sqrt(1 - phi^2))
  
  
  for(i in 2:length(times)) {
    
    x[i] <- rnorm(1, c + phi*x[i-1], sigma)
    
  }
  
  if(epsilon != 0) {
    error <- rnorm(length(times), 0, epsilon)
  } else {
    error <- 0
  }
  
  
  return(x + error)
}

## MISC ********************** ####

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
    set_colnames(c(names(cmeans), id.vars))
  
  
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


mean_and_var <- function(x) {
  
  m <- mean(x)
  v <- var(x)
  
  c(mean = m, var = v)
}


# Random truncated normal number generator
rnorm_trunc <- function(n, mean, sd, trunc = 0) {
  
  x <- c(0, n)
  
  for(i in 1:n) {
    y <- trunc - 1 
    while(y < trunc) {
      y <- rnorm(1, mean, sd)
    }
  
    x[i] <- y
    
  }
  return(x)
}


get_stan_results <- function(samples, parameter, regex = TRUE) {
  
  df <- summary(samples)$summary
  
  if(regex) {
    parameter_regex <- str_detect(rownames(df), parameter)
    
    res <- df[parameter_regex, c("50%", "25%", "75%", "2.5%", "97.5%")] %>% 
      as.data.frame()
  } else {
    res <- df[parameter, c("50%", "25%", "75%", "2.5%", "97.5%")] %>% 
      as.data.frame()
  }
  
  
  
  if(ncol(res) == 1) {
    res <- res %>%
      t
  }
  
  res <- res %>% 
    as.data.frame() %>% 
    set_colnames(c("mode", "lower25", "upper75", "lower2.5", "upper97.5"))
  
  if(nrow(res) == 1) {
    res <- res %>% set_rownames(gsub("\\\\", "", parameter))
  }
  
  return(res) 
}


softMax <- function(x) {
  exp(x)/(sum(exp(x)))
}

thin <- function(df, time = "x", modulo = 1) {
  df[(df[, time] %% modulo) == 0, ]
}


grep_rows <- function(df, pattern) {
  
  df[grep(pattern, rownames(df)), ]
  
}

# matrix to tibble with rownames to column
m_neat <- function(x, colnames = NULL) {
  x <- x %>%
    as.data.frame() %>% 
    rownames_to_column()
  
  if(!is.null(colnames)) {
    x <- x %>% set_colnames(colnames)
  }
  
  return(x)
}

# Generate 2d oup
oup_2d <- function(length = 100, Y0, mu, B, G) {
  
  Y <- matrix(NA, nrow = length, ncol = 2)
  Y[1, ] <- Y0
  exp_mB <- expm(-B) %>% as.matrix()
  
  for(i in 2:length) {
    
    
    Mean <- mu + exp_mB%*%(Y[i-1, ] - mu)
    Cov <- G - exp_mB%*%G%*%exp_mB %>%
      as.matrix()
    
    
    
    Y[i, ] <- rmvnorm(n = 1, mean = Mean, sigma = Cov)
    
  }
  
  colnames(Y) <- c("x", "y")
  
  return(Y)
}

# Generate multi d oup
oup_multi_d <- function(length = 100, Y0, mu, B, G) {
  
  D <- nrow(B)
  
  Y <- matrix(NA, nrow = length, ncol = D)
  Y[1, ] <- Y0
  exp_mB <- expm(-B) %>% as.matrix()
  
  for(i in 2:length) {
    
    
    Mean <- mu + exp_mB%*%(Y[i-1, ] - mu)
    Cov <- G - exp_mB%*%G%*%exp_mB %>%
      as.matrix()
    
    
    
    Y[i, ] <- rmvnorm(n = 1, mean = Mean, sigma = Cov)
    
  }
  
  Y <- Y %>% as.data.frame()
  
  return(Y)
}



# Covariance matrix
covariance <- function(lambda, sigma, intervals) {
  n <- length(intervals)
  kernel <- matrix(NA, nrow=length(intervals), ncol=length(intervals))
  for(i in 1:n) {
    for(j in 1:n) {
      kernel[i,j] <- (sigma^2/(2*lambda))*exp(-lambda*abs(intervals[i]-intervals[j]))
    }
  }
  
  return(kernel)
}

# generate a student t set
generate_student_set <- function(n_series, student_df, mu, lambda, sigma, kappa = NULL, intervals, seed=1) {
  set.seed(seed)
  n <- length(intervals)
  
  if(!is.null(kappa)) {
    sigma = sqrt(2*kappa*lambda)
  }
  
  

  
  # observations
  obs <- list()
  for(i in 1:n_series) {
    
    # student or gaussian oup?
    if(is.finite(student_df)) {
      
      shape = ((student_df-2)/student_df)*covariance(lambda[i], sigma[i], intervals)
      obs[[i]] <- rmvt(1, df=student_df,  delta = rep(mu[i], n), sigma = shape, type="shifted") %>% as.matrix()
      
    } else {
      
      shape = covariance(lambda[i], sigma[i], intervals)
      obs[[i]] <- rmvt(1, df=student_df,  delta = rep(mu[i], n), sigma = shape, type="shifted") %>% as.matrix()
      
    }
    
    
  }
  
  # make matrix
  obs <- do.call(rbind, obs)
  
  if(n_series == 1) {
    obs <- obs %>% t() %>% as.vector()
  }
  
  return(list(Y=obs,
              N=n_series,
              time=intervals,
              T=n,
              student_df=ifelse(is.finite(student_df),
                                student_df,
                                -99),
              mu_values = mu,
              lambda_values = lambda,
              sigma_values = sigma,
              kappa_values = kappa))
}

restricted_rlnorm <- function(n, meanlog = 0, sdlog = 1, lower = 0, upper = 1) {
  
  vec <- c()
  for(i in 1:n) {
    
    x <- 1
    while(x >= upper | x <= lower) {
      x <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)
    }
    
    vec[i] <- x
  }
  
  return(vec)
  
}

# Is a value inside interval?
inside <- function(x, interval, or_equal = FALSE) {
  
  if(length(interval) != 2) {
    stop("Bad interval")
  }
  
  if(or_equal == TRUE) {
    is_in <- (x >= interval[1]) & (x <= interval[2])
  } else {
    is_in <- (x > interval[1]) & (x < interval[2])
  }
  
  
  return(unname(is_in))
  
}

# Get smallest IQR inside which the simulation value is
get_IQRs <- function(stanfit, parameter, parameter_values) {
  
  prob = 100:1/100
  
  posterior <- rstan::extract(stanfit, pars = parameter)[[1]] %>% 
    as.data.frame()
  
  if(dim(posterior)[2] == 1 & length(parameter_values) > 1) {
    
    posterior <- rep(posterior, length(parameter_values)) %>% do.call(cbind, .)
    
  }
  
  IQR_df <- lapply(1:length(parameter_values), function(i) {
    
    interval_membership <- sapply(1:length(prob), FUN = function(j) {
      interval <- posterior[, i] %>% quantile(probs = c(.5 - prob[j]/2, .5 + prob[j]/2))
      inside(parameter_values[i], interval)
    }) %>% set_names(prob)
    
    interval <- which(interval_membership == TRUE) %>%
      tail(n=1) %>%
      names
    
    if(sum(interval_membership == FALSE) == length(prob)) {
      interval <- 1
    }
    
    IQR50 <- quantile(posterior[, i], probs = c(.25, .75)) %>% unname()
    
    mean <- posterior[, i] %>% mean
    median <- posterior[, i] %>% median
    sd <- posterior[, i] %>% sd

    return(c(parameter = parameter,
             index = i,
             simulation_value = parameter_values[i],
             IQR_min = interval,
             IQR50_lower = IQR50[1],
             IQR50_upper = IQR50[2],
             mean = mean,
             median = median,
             sd = sd))
    
  }) %>%
    do.call(rbind, .) %>% 
    as.data.frame() %>% 
    cbind(.,
          n_series = n_series,
          n_observations = n_observations)
  
  IQR_df$IQR_min <- IQR_df$IQR_min %>% as.character() %>% as.numeric
  IQR_df$sd <- IQR_df$sd %>% as.character() %>% as.numeric
  IQR_df <- IQR_df %>%
    mutate(mean_IQR = mean(IQR_min), mean_sd = mean(sd))


  return(IQR_df)
  
}

# Get posterior median
get_posterior_mode <- function(stanfit, parameter) {
  
  posterior <- rstan::extract(stanfit, pars = parameter)[[1]] %>% 
    as.data.frame()
  
  medians <- apply(posterior, 2, FUN = function(x) quantile(x, probs = c(.5))) %>% 
    unname
  
  return(medians)
}

# Generate positive normal numbers
restricted_rnorm <- function(n, mean, sd, lower = 0) {
  vec <- c()
  for(i in 1:n) {
    
    x <- lower -1
    while(x < lower) {
      x <- rnorm(1, mean, sd)
    }
    
    vec[i] <- x
  }
  
  return(vec)
}

# Generate gaussian process set with stan
genera_gp_set_stan <- function(n_series, alpha, rho, sigma, intervals, stan_model, seed = 1) {
  set.seed(seed)
  
  x_total <- intervals
  
  N_total <- length(x_total)
  
  # Generate data with stan
  obs <- lapply(1:n_series, function(i) {
    
    gp_simulator_data <- list(N = N_total, x = x_total, alpha = alpha[i], rho = rho[i], sigma = sigma[i])
    
    seed <- runif(1, 1, 1000)
    
    samples <- sampling(stan_model,
                        data = gp_simulator_data, 
                        iter=1,
                        chains=1,
                        seed=seed,
                        algorithm="Fixed_param")
    
    fs <- grep("f", rownames(summary(samples)$summary), value = T)
    ys <- grep("y", rownames(summary(samples)$summary), value = T)
    
    
    df <- data.frame(f = summary(samples)$summary[fs, "mean"],
                     y = summary(samples)$summary[ys, "mean"]) %>% 
      set_rownames(NULL)
    
    
  })
  
  
  # observations
  y <- obs %>% lapply(function(x) {
    x[, "y"]
  }) %>% do.call(rbind, .)    
  
  # latent process
  f <- obs %>% lapply(function(x) {
    x[, "f"]
  }) %>% do.call(rbind, .)
  
  if(n_series == 1) {
    y <- y %>% t() %>% as.vector()
    f <- f %>% t() %>% as.vector()
  }
  
  return(list(y=y,
              f = f,
              S=n_series,
              x=intervals,
              N=length(intervals),
              alpha_values = alpha,
              rho_values = rho,
              sigma_values = sigma))
}


# One long into many short
divide_long_series <- function(long_data, into = 10) {
  
  if(length(long_data$y) %% into != 0) {
    stop(paste0("Number of observations not divisible with ", into))
  }
  
  y <- long_data$y %>% 
    matrix(nrow = into, byrow = TRUE)
  
  x <- 1:ncol(y)
  
  N <- ncol(y)
  
  S <- into
  
  return(list(N = N, S = S, x = x, y = y))  
}

# Prior data

prior_data <- function(y_shape = "normal", y1, y2, x_shape = "invgamma", x1, x2, x = 1:100/10, y = 1:100/10) {
  
  prior <- data.frame(x = rep(x, length(y)),
                      y = rep(y, each = length(x)),
                      value = NA)
  
  for(i in x) {
    print(paste0(which(i == x), "/", length(x)))
    for(j in y) {
      
      condition <- i == prior$x & j == prior$y
      prior[condition, "value"] <- dinvgamma(i, x1, x2)*dnorm(j, y1, y2)
      
      
    }
  }
  
  
  return(prior)
  
  
}

#### DUMP ####

# Alternative generator
ou_simulator <- function (T, mu, lambda, sigma, x0 = NULL, seed=1) {
  seed <- seed
  x <- c()
  
  # Initial value
  if (is.null(x0)) {
    x[[1]] <- mu + rnorm(1) * sqrt((sigma^2)/(2*lambda))
  } else {
    x[[1]] <- x0
  }
  
  # Consecutive values
  for (t in 2:T) {
    x[[t]] <- mu - (mu - x[[t-1]]) * exp(-lambda) + rnorm(1) * sqrt((sigma^2)/(2*lambda) * (1 - exp(-2 * lambda)))
  }
  
  #gompertz assumptions, poisson meanl 
  list(observations = x, T=T, n_series=1, samples_per_series=as.array(T), time=1:T)
  
  
}

# Generate gaussian multivariate normal data
generate_mvrnormal <- function(n_series, mu, lambda, sigma, intervals) {
  
  n <- length(intervals)
  
  # covariance matrix
  kernel <- matrix(NA, nrow=n, ncol=n)
  for(i in 1:n) {
    for(j in 1:n) {
      kernel[i,j] <- (sigma^2/(2*lambda))*exp(-lambda*abs(intervals[i]-intervals[j]))
    }
  }
  
  # observations
  obs <- list()
  for(i in 1:n_series) {
    obs[[i]] <- mvrnorm(1, mu = rep(mu, n), Sigma = kernel) %>% as.matrix()
  }
  
  # make matrix
  obs <- do.call(cbind, obs) %>% t()
  
  if(n_series == 1) {
    obs <- obs %>% t() %>% as.vector()
  }
  
  return(list(Y=obs, N=n_series, time=intervals, T=n))
}




# original generator
generateStanData <- function(kappa,
                             lambda,
                             mu,
                             intervals,
                             t.df = Inf,
                             seed = 1){
  set.seed(seed)
  N <- length(intervals)
  # lv.variates <- rnorm(N)
  dt <- outer(intervals,intervals,function(x,y) abs(x-y))
  x <- kappa * exp(-lambda*dt)
  L <- chol(x)
  
  scale <- if(is.finite(t.df)) rep(sqrt(rgamma(1,t.df/2,(t.df-2)/2)),each=N) else 1
  out.data <- list()
  out.data$latent_value <- as.vector(t(L) %*% (rnorm(N) * scale)) + mu
  
  out.data$value <- rpois(N,exp(out.data$latent_value))
  
  
  out.data$time <- intervals
  out.data$replicates <- 1L
  out.data$replicate_samples <- array(length(intervals))
  out.data$NSMPL <- length(intervals)
  out.data
}

# original generator w different parameter names
generate_a_series <- function(kappa=NULL, sigma=NULL, lambda, mu, intervals, t.df = Inf, seed = 1){
  set.seed(seed)
  N <- length(intervals)
  lv.variates <- rnorm(N)
  dt <- outer(intervals,intervals,function(x,y) abs(x-y))
  
  if(!is.null(kappa)) {
    x <- kappa * exp(-lambda*dt)
  } else {
    x <- (sigma^2)/(2*lambda) * exp(-lambda*dt)
  }
  
  
  L <- chol(x)
  
  scale <- if(is.finite(t.df)) rep(sqrt(rgamma(1,t.df/2,(t.df-2)/2)),each=N) else 1
  out.data <- list()
  out.data$Y <- as.vector(t(L) %*% (rnorm(N) * scale)) + mu
  # out.data$observations <- rpois(N,exp(as.vector(t(L) %*% (rnorm(N) * scale))+mu))
  out.data$time <- intervals
  out.data$n_series <- 1L
  out.data$samples_per_series <- array(length(intervals))
  out.data$T <- length(intervals)
  out.data$student_df <- t.df
  out.data
}

generate_n_series <- function(n, kappa=NULL, sigma=NULL, lambda, mu, intervals, t.df = Inf, seed = 1, fix_mu=NULL, fix_kappa_log=NULL) {
  
  s_list <- list()
  
  for(i in 1:n) {
    s_list[[i]] <- generate_a_series(kappa=kappa, sigma=sigma, lambda=lambda, mu=mu, intervals=intervals, t.df = t.df, seed = seed + i)
    
    if(!is.null(fix_mu)) {
      s_list[[i]][["mu"]] <- as.array(fix_mu)
    }
    if(!is.null(fix_kappa_log)) {
      s_list[[i]][["kappa_log"]] <- as.array(fix_kappa_log)
    }
    
  }
  
  if(length(s_list) == 1) {
    return(s_list[[1]])
  }
  s_list
  
}

# Generate one long for the original
generate_data_original <- function(kappa=0.1, lambda=0.5, mu=5, intervals, t.df = Inf, seed = 1){
  set.seed(seed)
  n_obs <- length(intervals)
  lv.variates <- rnorm(n_obs)
  dt <- outer(intervals,intervals,function(x,y) abs(x-y))
  x <- kappa * exp(-lambda*dt)
  L <- chol(x)
  
  scale <- if(is.finite(t.df)) rep(sqrt(rgamma(1,t.df/2,(t.df-2)/2)),each=n_obs) else 1
  
  data <- list()
  Y <- matrix(NA, nrow=n_series, ncol=length(intervals))
  
  
    X_latent <- as.vector(t(L) %*% (rnorm(n_obs) * scale)) + mu
    Y <- rpois(n_obs,exp(X_latent))
  
  
  data[["n_series"]] <- 1           # number of time series
  data[["samples_per_series"]] <- c(length(intervals)) # number of time series
  data[["T"]] <- length(intervals)  # number of time points
  data[["observations"]] <- Y                  # observations
  data[["time"]] <- intervals       # observation times
  
  return(data)
}

# Concatenate several short series

generate_data_set2 <- function(n_series, kappa, lambda, mu, intervals, t.df = Inf, seed = 1){
  
  setlist <- list()
  
  for(i in 1:n_series) {
    setlist[[i]] <- generateStanData(kappa,
                                     lambda,
                                     mu,
                                     intervals,
                                     t.df = Inf,
                                     seed = i)
  }
  
  T <- n_series*length(intervals)
  samples_per_series <- rep(length(intervals), n_series)
  observations <- c()
  time <- c()
  for(i in 1:length(setlist)) {
    observations <- c(observations, setlist[[i]][["value"]])
    time <- c(time, setlist[[i]][["time"]])
  }
  
  
  return(list(T=T, samples_per_series=samples_per_series, observations=observations, time=time, n_series=n_series))
  
}

# Generate a data set
generate_data_set <- function(kappa=0.1, lambda=0.5, mu=5, intervals, n_series,
                              t.df = Inf, seed = 1){
  
  set.seed(seed)
  n_obs <- length(intervals)
  lv.variates <- rnorm(n_obs)
  dt <- outer(intervals,intervals,function(x,y) abs(x-y))
  x <- kappa * exp(-lambda*dt)
  L <- chol(x)
  
  scale <- if(is.finite(t.df)) rep(sqrt(rgamma(1,t.df/2,(t.df-2)/2)),each=n_obs) else 1
  
  data <- list()
  Y <- matrix(NA, nrow=n_series, ncol=length(intervals))
  
  for(i in 1:n_series) {
    X_latent <- as.vector(t(L) %*% (rnorm(n_obs) * scale)) + mu
    Y[i, ] <- rpois(n_obs,exp(X_latent))
  }
  
  data[["N"]] <- n_series           # number of time series
  data[["T"]] <- length(intervals)  # number of time points
  data[["Y"]] <- Y                  # observations
  data[["time"]] <- intervals       # observation times
  
  return(data)
}


# concatenate series
concatenate_series <- function(s_list) {
  
  len <- length(s_list)
  
  if(len == 1) {
    return(s_list[[1]])
  }
  
  obs <- c()
  t <- c()
  s <- c()
  for(i in 1:length(s_list)) {
    obs <- c(obs, s_list[[i]][["observations"]])
    t <- c(t, s_list[[i]][["time"]])
    s <- c(s, s_list[[i]][["samples_per_series"]])
  }
  
  kappa_log <- rep(s_list[[1]][["kappa_log"]], len)
  mu <- rep(s_list[[1]][["mu"]], len)
  
  T <- length(obs)
  
  
list(observations=obs, time=t, n_series=len, samples_per_series=s, T=T, kappa_log=kappa_log, mu=mu)
  
}



# HPDI
success_rate_HPDI <- function(stan_fit, parameter="lambda", real_value=1, alpha=0.95) {
  
  pos <- rstan::extract(stan_fit)[[parameter]]
  
  success <- 0
  for(i in 1:ncol(pos)) {
    hpdi <- HPDI(pos[,i], prob=alpha)
    
    if(hpdi[1] < real_value && real_value < hpdi[2]) {
      success <- success + 1
    }
    
  }
  success/ncol(pos)
}

# Percentiles
success_rate_quantiles <- function(stan_fit, parameter="lambda", real_value=1, success_limit=0.95) {
  
  pos <- rstan::extract(stan_fit)[[parameter]]
  
  success <- 0
  for(i in 1:ncol(pos)) {
    q <- quantile(pos[,i], probs = c(1-success_limit, 0.5, success_limit))
    
    if(q[1] < real_value && real_value < q[3]) {
      success <- success + 1
    }
    
  }
  success/ncol(pos)
}





# Get df with lenght, mean and 50% intervals 
# samples_df <- function(sample_list) {
#   res_df <- matrix(NA, length(sample_list), 4)
#   
#   colnames(res_df) <- c("Observations", "lower25", "mean", "upper75")
#   
#   j <- 1
#   for(i in sample_list) {
#     res <- summary(i)$summary[grep("lambda\\[", rownames(summary(i)$summary)),c("25%", "50%", "75%")]
#     
#     res_df[j, c("lower25", "mean", "upper75")] <- res
#     j <- j + 1
#   }
#   
#   res_df[, "Observations"] <- as.numeric(names(sample_list))
#   res_df <- res_df %>% as.data.frame()
#   
#   res_df
# }

samples_df2 <- function(sample_list) {
  
  # make data frame for results
  res_df <- matrix(NA, 2*length(sample_list), 5) %>% set_colnames(c("Series", "Observations", "lower25", "mean", "upper75"))
  
  
  j <- 1
  for(i in sample_list) {
    res <- summary(i)$summary[grep("lambda\\[", rownames(summary(i)$summary)),c("25%", "50%", "75%")]
    
    res_df[c(j, j+1), c("lower25", "mean", "upper75")] <- res
    j <- j + 2
  }
  
  
  res_df[, "Observations"] <- rep(as.numeric(names(sample_list)), each=2)
  res_df <- res_df %>% as.data.frame()
  res_df[, "Series"] <- rep(c("A","B"), length(sample_list))
  
  res_df
}

# plot the posterior intervals
plot_posteriors <- function(res_df, sim_value, par) {
  p <- ggplot(res_df, aes(x=Observations, y=mean)) + 
    geom_errorbar(aes(ymin=lower25, ymax=upper75), width=1) +
    geom_line() +
    geom_point() + geom_hline(yintercept = sim_value, linetype="dashed") + labs(y="Estimate", x="Length") + theme_bw() 
  # + scale_y_continuous(limits = c(0, 6.5))
  
  p
} 

# hierarchial plotter
plot_hier_posteriors <- function(res_df, sim_value) {
  p <- ggplot(res_df, aes(x=Observations, y=mean, group=Series, color=Series)) + 
    geom_errorbar(aes(ymin=lower25, ymax=upper75), width=1) +
    geom_line() +
    geom_point()  + geom_hline(yintercept = sim_value, linetype="dashed", color="black") + labs(y="Estimate", x="Length", title="Lambda estimates with 50% error bars vs. series length") + theme_bw() + scale_color_manual(values=c("#f1a340", "#998ec3")) + scale_y_continuous(limits = c(0, max(res_df$upper)))
  
  p
}




samples_df_par <- function(sample_list, par) {
  
  res_df <- matrix(NA, length(sample_list), 4)
  
  colnames(res_df) <- c("length", "lower25", "mean", "upper75")
  
  j <- 1
  for(i in sample_list) {
    res <- summary(i)$summary[par,c("25%", "50%", "75%")]
    
    res_df[j, c("lower25", "mean", "upper75")] <- res
    j <- j + 1
  }
  
  res_df[, "length"] <- as.numeric(names(sample_list))
  res_df <- res_df %>% as.data.frame()
  
  res_df
  
}

samples_df_par2 <- function(sample_list, par) {
  
  res_df <- matrix(NA, 1, 4)
  
  colnames(res_df) <- c("length", "lower25", "mean", "upper75")
  

  for(i in sample_list) {
    
    s <- summary(i)$summary
    
    res <- s[grepl(paste0(par, "\\["), rownames(s)),c("25%", "50%", "75%")]
    
    if(class(res) == "numeric") {
      res_df <- res_df %>% rbind(c(1, res))
    } else {
      res_df <- res_df %>% rbind(cbind(nrow(res), res))
    }
   
  }
  
  res_df <- res_df[-1,]
  res_df <- res_df %>% as_tibble()
  
  res_df
  
}




# sigma generator
generate_a_sigma_series <- function(kappa,
                              lambda,
                              mu,
                              intervals,
                              t.df = Inf,
                              seed = 1){
  set.seed(seed)
  N <- length(intervals)
  lv.variates <- rnorm(N)
  dt <- outer(intervals,intervals,function(x,y) abs(x-y))
  x <- (.5^2)/lambda * exp(-lambda*dt)
  L <- chol(x)
  
  scale <- if(is.finite(t.df)) rep(sqrt(rgamma(1,t.df/2,(t.df-2)/2)),each=N) else 1
  out.data <- list()
  out.data$observations <- as.vector(t(L) %*% (rnorm(N) * scale)) + mu
  # out.data$observations <- rpois(N,exp(as.vector(t(L) %*% (rnorm(N) * scale))+mu))
  out.data$time <- intervals
  out.data$n_series <- 1L
  out.data$samples_per_series <- array(length(intervals))
  out.data$T <- length(intervals)
  out.data
}





# get plots form short series results
plot_short_series_posteriors <- function(res_df) {
  
  breaks <- res_df$length %>% unique()

  p <- ggplot(res_df, aes(x=length, y=mean)) + geom_point() + theme_bw() + labs(y="Posterior mean",x="Number of series" ) + scale_x_continuous(breaks = breaks)
  
  p
}


concatenate_series_matrix <- function(s_list) {
  
  obs <- s_list[[1]][["observations"]]
  
  for(i in 2:length(s_list)) {
    obs <- rbind(obs, s_list[[i]][["observations"]])
  }
  
  return(list(Y=obs, time=s_list[[1]][["time"]], N=length(s_list), T=length(s_list[[1]][["time"]])))
  
}



oup_invG_lambda <- function(n=compare_n_series, shape=2, scale=0.5) {
  vec <- c()
  for(i in 1:n) {
    
    x <- 1
    while(x >= 1) {
      x <- MCMCpack::rinvgamma(1, shape = shape, scale = scale)
    }
    
    vec[i] <- x
  }
  
  return(vec)
}

r_non_negative_normal <- function(n, mean, sd) {
  vec <- c()
  for(i in 1:n) {
    
    x <- -1
    while(x < 0) {
      x <- rnorm(1, mean = mean, sd = sd)
    }
    
    vec[i] <- x
  }
  
  return(vec)
}

oup_G_lambda <- function(n=compare_n_series, shape=2, scale=4) {
  vec <- c()
  for(i in 1:n) {
    
    x <- 1
    while(x >= 1) {
      x <- rgamma(1, shape = shape, scale = scale)
    }
    
    vec[i] <- x
  }
  
  return(vec)
}




plot_pos_sequence <- function(stan_list, par="lambda") {

  
  # get samples
  pos <- lapply(names(stan_list), function(x) rstan::extract(stan_list[[x]], par)[[1]]) %>% set_names(names(stan_list))
  
  pos <- do.call(cbind, pos) %>% set_colnames(names(stan_list)) %>% melt %>% mutate(observations=as.factor(Var2))
  
  # generate plot
  p <- ggplot(pos, aes(x=value, color=observations)) + stat_density(aes(color=observations), geom="line", position = "identity") + labs(y="Density", x=par %>% str_to_title())
  # + scale_color_manual(type = "seq", palette = 4)
  
  return(p)
}


plot_pos <- function(stan_obj, par="lambda") {
  # get samples
  pos <- rstan::extract(stan_obj, par)[[1]] %>% as_tibble %>% melt
  # generate plot
  p <- ggplot(pos, aes(x=value, color=variable)) + stat_density(geom="line", position = "identity") + labs(y="Density", x=par %>% str_to_title()) + xlab(TeX(paste0("$\\", par))) + scale_color_manual(name="Sample size", values="black", labels="100")
  return(p)
}
