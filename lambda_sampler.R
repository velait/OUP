# Compute completely pooled models posterior = likelihood*prior by hand
# Assume delta_t = 1 for simplicity

# OUP generator
generate_oup <- function(x0 = NULL, n = 100, lambda = .5, mu = 0, sigma = .25, seed = 1) {
  set.seed(seed)
  kappa <- (sigma^2)/(2*lambda)
  
  # First observation
  if(is.null(x0)) {
    x <- rnorm(1, mu, kappa)
  } else {
    x = x0
  }
  
  
  # Consecutive observations
  for(i in 2:n) {
    
    step <- rnorm(1,
                  mean = mu - (mu - x[i-1])*exp(-lambda),
                  sd = kappa*(1 - exp(-2*lambda)))
    
    x <- c(x, step)
    
  }
  
  return(x)
}


# Likelihood (normal, not t-distribution)
oup_log_likelihood <- function(x, lambda = .5, mu = 0, sigma = .25) {
  
  kappa <- (sigma^2)/(2*lambda)
  
  
    # First observation is sampled from stationary distribution
    lh <- dnorm(x[1], mean = mu, sd = sqrt(kappa))
    
    # lh_next <- 
    
    # Consecutive observations given the previous
    for(i in 2:length(x)) {
      
      step <- dnorm(x[i],
                  mean = mu - (mu - x[i-1])*exp(-lambda),
                  sd = kappa*(1 - exp(-2*lambda)))
      
      lh <- c(lh, step)
      
    }
  
  llh <- prod(lh) %>% log
    
  return(llh)
}


## Sample


# Priors
# mu <- rnorm(n = samples, mean = 0, sd = 1)
# sigma <- restricted_rnorm(n = samples , mean = .5, sd = 1, lower = 0)



resolution <- 100



## One series ************************************************* ####
 
series <- generate_oup(n = 50)
# marginalize over mu and sigma
lambda_marginal <- sapply(seq(0.0001, 1, length.out = resolution), function(l) {
  print(l)
  
  minus_mu <- sapply(seq(-1, 1, length.out = resolution), function(m) {
    
    minus_sigma <- sapply(seq(0.001, 1, length.out = resolution), function(s) {
      oup_log_likelihood(x = series, lambda = l, mu = m, sigma = s)
    }
    )
    
    return(minus_sigma %>% exp %>% sum %>% log)
    
  }
  )
  
  return(minus_mu %>% exp %>% sum %>% log)
  
})

## test likelihood for varying lambda and fixed mu, sigma ***** ####

# grow length of series
lh_grid <- lapply(10*12:20, function(i) {
  series <- generate_oup(n = i)
  
  single_series_lh <- sapply(seq(0.0001, 1, length.out = resolution), function(x) {
    oup_log_likelihood(series, lambda = x, sigma = .25, mu = 0)
  })
  
}) %>%
  do.call(cbind, .) %>% 
  as.data.frame()

# normalize
lh_grid <- apply(lh_grid, 2, FUN = function(i) {
  i <- exp(i)
  i/sum(i)
})

lh_grid_plot <- lh_grid %>%
  melt %>%
  cbind(x = rep(seq(0.0001, 1, length.out = resolution), 9)) %>%
  ggplot(aes(x = x, y = value, color = Var2)) +
  geom_line() +
  scale_color_brewer(type = "seq")

## Multiple series ******************************************** ####

many_series <- lapply(1:20, function(i) generate_oup(n = 25)) %>% 
  do.call(rbind, .)

# Hold sigma and mu fixed
lambda_marginal_many <- lapply(1:20, function(x) {
  print(x)
  many_series <- lapply(1:x, function(i) generate_oup(n = 25)) %>% 
    do.call(rbind, .)
  
  samples <- sapply(seq(0.0001, 1, length.out = resolution), function(l) {
    print(l)
    
    minus_mu <- sapply(0, function(m) {
      
      
      
      
      minus_sigma <- sapply(.25, function(s) {
        
        sapply(1:nrow(many_series), FUN = function(i) {
          oup_log_likelihood(x = many_series[i, ], lambda = l, mu = m, sigma = s)
        }) %>% sum
        
        
      }
      )
      
      
      
      return(minus_sigma)
      
    }
    )
    
    return(minus_mu)
    
  })
  
  samples/sum(samples)
})


lambda_marginal_many_df <- lapply(1:length(lambda_marginal_many), function(x) {
  
  df <- data.frame(n_series = x,
                   y = lambda_marginal_many[[x]],
                   x = seq(0.0001, 1, length.out = resolution))
  
}) %>% do.call(rbind, .)

lambda_marginal_many_df %>%
  filter((n_series %% 4) == 0) %>% 
  ggplot(aes(x = x, y = y, color = as.factor(n_series))) +
  geom_line() +
  scale_color_brewer(type = "seq", palette = 2)
 


lambda_marginal_many <- sapply(seq(0.0001, 1, length.out = resolution), function(l) {
  print(l)
  
  
  minus_mu <- sapply(seq(-1, 1, length.out = resolution), function(m) {
    
    
    
    
    minus_sigma <- sapply(seq(0.001, 1, length.out = resolution), function(s) {
      
      sapply(1:nrow(many_series), FUN = function(i) {
        oup_log_likelihood(x = many_series[i, ], lambda = l, mu = m, sigma = s)
      }) %>% sum
      
      
    }
    )
    
    
    
    return(minus_sigma %>% exp %>% sum %>% log)
    
  }
  )
  
  return(minus_mu %>% exp %>% sum %>% log)
  
})


## test likelihood many series ******************************** ####

# grow amount if series, fix length at 10
lh_grid <- lapply(5:10, function(i) {
  
  # generate data
  series <- lapply(1:i, function(x) {
    generate_oup(n = 50, seed = x)
  }
  ) %>% do.call(cbind, .) %>% as.data.frame()
  
  
  # lambda grid likelihoods
  single_series_lh <- sapply(seq(0.0001, 1, length.out = resolution), function(x) {
    
    apply(series, 2, FUN = function(y) {
      oup_log_likelihood(y, lambda = x, sigma = .25, mu = 0)
    }) %>% sum
    
  })
  
}) %>%
  do.call(cbind, .) %>% 
  as.data.frame()

# normalize
lh_grid <- apply(lh_grid, 2, FUN = function(i) {
  i <- exp(i)
  i/sum(i)
})

lh_grid_plot <- lh_grid %>%
  melt %>%
  cbind(x = rep(seq(0.0001, 1, length.out = resolution), 6)) %>%
  ggplot(aes(x = x, y = value, color = Var2)) +
  geom_line() +
  scale_color_brewer(type = "seq")


lh_grid %>%
  melt %>%
  cbind(x = rep(seq(0.0001, 1, length.out = resolution), 6)) %>%
  ggplot(aes(x = x, y = value, color = variable)) +
  geom_line() +
  scale_color_brewer(type = "seq")
