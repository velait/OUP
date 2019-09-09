set.seed(1)
times <- seq(from = 0, to = 250, by = 0.1)
N_series <- 100


# Simulate data
ews_set <- lapply(1:N_series, function(j) {
  
  r <- 1
  K <- 10
  cs <- 1
  h <- 1
  sigma <- 0.07
  
  
  y <- rep(NA, length(times))
  y[1] <- 8
  
  cs <- 1 + (2.7 - 1)*(times/max(times))
  
  
  for(i in 2:length(times)) {
    
    
    # c <- cs[i]
    
    dt <- times[(i-1):i]
    
    y[i] <- ews_generator(y[i-1], dt, c = cs[i], milstein = T)[2]
    
    
  }
  
  
  return(log(y + 1))
  
  
}) %>%
  do.call(cbind, .) %>%
  cbind(x = times) %>%
  as.data.frame()


ews_set <- ews_set[1:(which((cs < 2.604)) %>% max()), ]

# ews_set <- ews_set %>% thin(modulo = )


# ews_set[, 89:101] %>% as.data.frame() %>% melt(id.vars = "x") %>% ggplot(aes(x = x, y = value)) + geom_line() + facet_wrap(~variable)





prediction_likelihoods <- lapply(1:N_series, function(i) {
  ou.lik <- function(x) {
    function(theta1,theta2,theta3) {
      n <- length(x)
      # dt <- deltat(x)
      dt <- .1
      -sum(dcOU(x=x[2:n], Dt=dt, x0=x[1:(n-1)],
                theta=c(theta1,theta2,theta3), log=TRUE))
    }
  }
  print(i)
  my_data <- ews_set[, i]
  
  
  # Proportion of time series used for a single window
  window_prop <- .5
  window_length <- (window_prop*length(my_data) - 1) %>% round
  N_window <- length(my_data) - window_length + 1
  
  my_res <- lapply(1:N_window, function(j) {
    
    # AR(1), get coefficients
    # arse <- ar(my_data[j:(j + window_length - 1)], order.max = 1, aic = FALSE)
    # coefs <- c(lambda = arse$ar, var = arse$var.pred, mu = arse$x.mean)
    
    
    
    
    ou.fit <- mle(ou.lik(my_data[j:(j + window_length - 1)]),
                  start=list(theta1=1,theta2=0.5,theta3=0.2),
                  method="L-BFGS-B",lower=c(0,1e-5,1e-3), upper=c(1,1,1))
    ou.coe <- coef(ou.fit)
    coefs <- c(ou.coe[2], ou.coe[1]/ou.coe[2], ou.coe[3]) %>% set_names(c("lambda", "mu", "sigma"))
    
    
    # Next time point likelihood
    next_m <- coefs["mu"] - (coefs["mu"] - my_data[j + window_length - 1])*exp(-coefs["lambda"])
    next_sd <- ((coefs["sigma"])/(2*coefs["sigma"]))*(1-exp(-2*coefs["lambda"]))
    LL <- dnorm(my_data[j + window_length],
          mean = next_m,
          sd = next_sd) 
    
    
    x <- c(log_likelihood = LL, window = j, series = i, from = j, to = (j + window_length - 1), pred_mean = next_m, pred_sd = next_sd, coefs)
    x["lambda_per_sigma"] <-  coefs["lambda"]/coefs["sigma"]
       
    x
    
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  
  my_res
  
  
}) %>% do.call(rbind, .)



# prediction_likelihoods %>% ggplot(aes(x = as.factor(window), y = log(log_likelihood), group = as.factor(series)))  + geom_line()




## Plot

# Likelihood
i <- 89

ggplot() +
  geom_line(data = ews_set[, i:ncol(ews_set)] %>%
              melt(id.vars = "x") %>%
              mutate(series = gsub("V", "", variable)) , aes(y = value, x = x)) + 
  geom_errorbar(data = prediction_likelihoods %>% filter(series >= i), 
              aes(x = (to + 1)/10, ymin = pred_mean.mu - pred_sd.sigma, ymax = pred_mean.mu + pred_sd.sigma), color = "royalblue", alpha = .1) + 
  facet_wrap(~series, scales = "free")  +
  geom_line(data =  prediction_likelihoods %>% filter(series >= i),
            aes(x = to/10, y = (log_likelihood)), color = "red") +
  guides(color=FALSE) +
  labs(x = "x", y = "y", subtitle = "Black = time series; Red = Log likelihood of the next observations; Blue = mean +- 1sd of the transition density")

  
  
# Resilience

# Standardize estimates
prediction_likelihoods$lambda_std <- (prediction_likelihoods$lambda - mean(prediction_likelihoods$lambda))/sd(prediction_likelihoods$lambda)
prediction_likelihoods$mu_std <- (prediction_likelihoods$mu - mean(prediction_likelihoods$mu))/sd(prediction_likelihoods$mu)
prediction_likelihoods$sigma_std <- (prediction_likelihoods$sigma - mean(prediction_likelihoods$sigma))/sd(prediction_likelihoods$sigma)
prediction_likelihoods$lambda_per_sigma_std <- (prediction_likelihoods$lambda_per_sigma - mean(prediction_likelihoods$lambda_per_sigma))/sd(prediction_likelihoods$lambda_per_sigma)


prediction_likelihoods$lambda_per_sigma <- prediction_likelihoods$lambda/(prediction_likelihoods$sigma^2)

# Standardize log_likelihood
for(i in 1:N_series) {
  
  lh <- prediction_likelihoods %>% filter(series == i) %>% pull(log_likelihood)
  
  prediction_likelihoods[prediction_likelihoods$series == i, "log_likelihood_std"] <- (lh - mean(lh, na.rm=T))/sd(lh, na.rm = T)
  
}


i <- 89

ggplot() +
  geom_line(data = ews_set[, i:ncol(ews_set)] %>%
              melt(id.vars = "x") %>%
              mutate(series = gsub("V", "", variable)) , aes(y = value, x = x)) + 
  geom_line(data = prediction_likelihoods[, c("to", "lambda_std", "sigma_std", "mu", "lambda_per_sigma_std", "series")] %>%
              melt(id.vars = c("to", "series")) %>% 
              filter(series >= i), aes(x = to/10, y = value, color = variable)) + 
  geom_line(data =  prediction_likelihoods %>% filter(series >= i),
            aes(x = to/10, y = (log_likelihood_std)), color = "red") + 
  facet_wrap(~series, scales = "free")







ggplot() +
  geom_line(data = ews_set[, i:ncol(ews_set)] %>%
              melt(id.vars = "x") %>%
              mutate(series = gsub("V", "", variable)) , aes(y = value, x = x)) + 
  geom_line(data = prediction_likelihoods[, c("to", "mu", "series")] %>%
              melt(id.vars = c("to", "series")) %>% 
              filter(series >= i), aes(x = to/10, y = value, color = variable)) + 
  # geom_line(data =  prediction_likelihoods %>% filter(series >= i),
  #           aes(x = to/10, y = (log_likelihood_std)), color = "red") + 
  facet_wrap(~series, scales = "free")
