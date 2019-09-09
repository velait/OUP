## CKLS 





## OZAKI with LAMPERTI *********************************** ####

ckls_drift <- function(x, alpha, beta) {
  
  alpha*(beta - x)
  
}

ckls_dispersion <- function(x, epsilon, gamma) {
  
  epsilon*x^gamma
  
}

ckls_transformation_transition_density <- function(x, dt, alpha, beta, epsilon, gamma) {
  
  drift <- ((alpha*beta)/epsilon)*x^(-gamma) - (alpha/epsilon)*x^(1 - gamma) - .5*epsilon*gamma*x^(gamma - 1)
  
  L <- (-(alpha*beta*gamma)/epsilon)*x^(-gamma - 1) + 
    (alpha*(gamma - 1)/epsilon)*x^(-gamma) + 
    .5*epsilon*gamma*(1-gamma)*x^(gamma - 2)
    
  M <- .5*( (alpha*beta*gamma*(1+gamma))/epsilon*x^(-gamma - 2) +
              (alpha*gamma*(1-gamma))/epsilon*x^(-gamma - 1) +
              .5*epsilon*gamma*(1-gamma)*(gamma-2)*x^(gamma - 3)
    )
    
  
  exp_L_dt <- exp(L*dt)
  
  Ax <-  x + drift*(exp_L_dt - 1)/L + M*(exp_L_dt - 1 - L*dt)/(L^2)
  B <- sqrt((exp_L_dt^2 - 1)/(2*L))
  
  rnorm(1, Ax, B)
  
}

ckls_transformation <- function(x, alpha, beta, epsilon, gamma) {
  
  z <- (1/(epsilon*(1-gamma)))*x^(1-gamma)
  
  return(z)
  
}

ckls_inv_transformation <- function(z, alpha, beta, epsilon, gamma) {
  
  x <- (epsilon*(1-gamma)*z)^(1/(1-gamma))
  
  return(x)
  
}

ckls_series <- function(x0 = 10, times = 1:1000, alpha, beta, epsilon, gamma) {
 
  dt <- times[2:length(times)] - times[1:(length(times) - 1)]
  
  x <- rep(0, length(times))
  z <- rep(0, length(times))
  
  
  x[1] <- x0
  z[1] <- ckls_transformation(x[1], alpha, beta, epsilon, gamma)
  
  for(i in 2:length(times)) {
    
    z[i] <- ckls_transformation_transition_density(z[i-1], dt[i-1], alpha, beta, epsilon, gamma)
    
    
    # Prevent from getting absorbed into zero
    # z[i] <- ifelse(z[i] >= 0, z[i], .1)
    
  }
  
  ckls_inv_transformation(z, alpha, beta, epsilon, gamma)
  
}



times <- seq(from = 1, to = 100, by = .01)


beta <- 20
init <- beta
alpha <- 1
epsilon <- .5
gamma <- .1


z <- lapply((1:9)/10 , function(i) ckls_series(x0 = init, times = times, alpha, beta, epsilon, gamma = i)) %>%
  do.call(cbind, .) %>%
  cbind(x = times) %>% 
  as.data.frame() %>% 
  set_colnames(c(as.character((1:9)/10), "x"))

z <- thin(z)

z %>%
  melt(id.vars = "x") %>% 
  filter(variable == "0.9") %>% 
  ggplot(aes(x = x, y = value, color = variable)) +
  geom_line()

## Stan ************************************************** ####




beta <- 10
init <- beta
alpha <- 2
epsilon <- sqrt(.5)
gamma <- .1



# Data
ckls_df <- ckls_series(x0 = init, times = times, alpha = alpha, beta = beta, epsilon = epsilon, gamma = gamma) %>% 
  cbind(y = ., x = times) %>% 
  thin() %>% 
  as.data.frame()


ckls_df %>% ggplot(aes(x = x, y = y)) + geom_line()


ckls_data <- list(y = ckls_df$y, 
                  x = ckls_df$x, 
                  N_obs = nrow(ckls_df))


# Model
ckls_model <- stan_model("lamperti/CKLS_shoji.stan")



ckls_samples <- sampling(ckls_model, ckls_data, 
                         iter = 1000, 
                         chains = 1)

## Euler-Maruyma ************************************* ####


em_ckls <- function(y0, grid, alpha, beta, epsilon, gamma, seed = NULL, milstein = FALSE, restrict_pos = FALSE) {
  
  if(is.null(seed)) {
    set.seed(sample(1:10000, 1))
  } else {
    set.seed(seed)
  }
  
  
  y <- y0
  for(i in 2:length(grid)) {
    
    y_prev <- y[i-1]
    delta_t <- grid[i] - grid[i-1]
    
    a <- y_prev + ckls_drift(y_prev, alpha, beta)*delta_t + ckls_dispersion(y_prev, epsilon, gamma)*rnorm(1, 0, sqrt(delta_t))
    
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




z <- lapply(1 + (1:9)/10 , function(i) em_ckls(y0 = init, grid = times, alpha, beta, epsilon, gamma = i, restrict_pos = T)) %>%
  do.call(cbind, .) %>%
  cbind(x = times) %>% 
  as.data.frame() %>% 
  set_colnames(c(as.character((1:9)/10), "x"))



z <- thin(z)

z %>%
  melt(id.vars = "x") %>% 
  ggplot(aes(x = x, y = value, color = variable)) +
  geom_line()


## Stan 
# Data
# 
beta <- 10
init <- beta
alpha <- 1
epsilon <- .5
gamma <- 0.5


ckls_df <- em_ckls(y0 = init, grid = times, alpha = alpha, beta = beta, epsilon = epsilon, gamma = gamma) %>% 
  cbind(y = ., x = times) %>% 
  thin(modulo = 1) %>%
  as.data.frame()


ckls_df %>% ggplot(aes(x = x, y = y)) + geom_line()


ckls_data <- list(y = ckls_df$y, 
                  x = ckls_df$x, 
                  N_obs = nrow(ckls_df))


# Model
ckls_euler_model <- stan_model("lamperti/CKLS_euler.stan")



ckls_samples <- sampling(ckls_euler_model, ckls_data, 
                         iter = 1000, 
                         chains = 1)
