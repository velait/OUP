seed <- 2
times <- seq(1, 100, by = 0.01)
c_pars <- c(-.3, 2, 0, 1, 1)
init <- 0

drift <- function(x, c_pars) {
  
  alpha <- c_pars[1]
  beta <- c_pars[2]
  lambda <- c_pars[3]
  epsilon <- c_pars[4]
  r <- c_pars[5]
  
  r*(alpha + beta*(x - lambda) - (x - lambda)^3)
  
}
dispersion <- function(x, c_pars) {
  
  alpha <- c_pars[1]
  beta <- c_pars[2]
  lambda <- c_pars[3]
  epsilon <- c_pars[4]
  r <- c_pars[5]
  
  # epsilon

  # epsilon*x^.5
    
  epsilon + dnorm(x, lambda, epsilon)
  # epsilon + dnorm(x, lambda, epsilon) + .5/(1 + exp(-(x-lambda)))
  
}




## Gaussian ********************** ####
## Sigmoid *********************** ####
## Constant ********************** ####
## Power ************************* ####
## Gaussian plus sigmoid  ******** ####

## Simulate with Euler-Maruyama


em_generator <- function(init, times, c_pars, seed = NULL, milstein = FALSE, restrict_pos = FALSE) {
  
  if(is.null(seed)) {
    set.seed(sample(1:10000, 1))
  } else {
    set.seed(seed)
  }
  
  alpha <- c_pars[1]
  beta <- c_pars[2]
  lambda <- c_pars[3]
  epsilon <- c_pars[4]
  r <- c_pars[5]
  
  y <- init
  
  for(i in 2:length(times)) {
    
    y_prev <- y[i-1]
    delta_t <- times[i] - times[i-1]
    
    a <- y_prev + drift(y_prev, c_pars)*delta_t + dispersion(y_prev, c_pars)*rnorm(1, 0, sqrt(delta_t))
    
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


process <- em_generator(init, times, c_pars, seed = 2)


process_df <- data.frame(y = process, x = times) %>% 
  thin(modulo = .05)


process_df %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_line()


# constant_plus_normal_logistic <- process
# 
# 
# data.frame(constant, constant_plus_normal, constant_plus_normal_logistic, x = times) %>% 
#   thin() %>% 
#   melt(id.vars = "x") %>% 
#   mutate(dispersion = variable) %>% 
#   ggplot(aes(x = x, y = value, color = dispersion)) +
#   geom_line(size = 1, alpha = .75)
#   
#   
# 
# data.frame(x = (-50:50)/10, 
#            constant = epsilon,
#            constant_plus_normal = 0,
#            constant_plus_normal_logistic = 0, 
#            stationary = 0) %>% 
#   mutate(constant_plus_normal = epsilon + dnorm(x, lambda, epsilon), 
#          constant_plus_normal_logistic = epsilon + dnorm(x, lambda, epsilon) + .5/(1 + exp(-(x-lambda))),
#          stationary  = cc_density(x, 1, -.3, 2, 0, 1)) %>% 
#   melt(id.vars = "x") %>% 
#   mutate(dispersion = variable) %>% 
#   ggplot(aes(x = x, y = value, color = dispersion)) + 
#   geom_line(size = 1, alpha = .75) +
#   coord_cartesian(ylim = c(0, 1.75))




## Stan 

dispersion_test_model <- stan_model("dispersion_tests/em_cusp_dispersion.stan")


dispersion_samples <- sampling(dispersion_test_model, list(N_obs = nrow(process_df),
                                                           x = process_df$x, 
                                                           y = process_df$y))
