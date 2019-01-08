# Functions

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
    

    return(c(parameter = parameter, index = i, simulation_value = parameter_values[i], IQR = interval))
    
  }) %>%
    do.call(rbind, .) %>% 
    as.data.frame() %>% 
    cbind(.,
          n_series = n_series,
          n_observations = n_observations)
  
  IQR_df$IQR <- IQR_df$IQR %>% as.character() %>% as.numeric
  IQR_df$simulation_value <- IQR_df$simulation_value %>% as.character() %>% as.numeric
  
  IQR_df <- IQR_df %>% mutate(mean_IQR = mean(IQR))
  
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
