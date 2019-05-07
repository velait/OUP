library(mvtnorm)


## MODELS ************************************************ ####

## LOTKA - VOLTERRA 
A <- matrix(c(0, -0.4, 0.1, 0), 2, 2, byrow = T)
r <- c(1.1, -0.4)
error_sq <- 0.01*diag(dim(A))

lotka_drift <- function(x, A, r) {

  x_new <- x*(r + A%*%x)
  
  return(x_new)
}



## TOY 

k <- c(1, 1, 1)
b <- c(1, 0.95, 1.05)
K <- matrix(c(0, 0.1, 0.1,
              0.1, 0, 0.1,
              0.1, 0.1, 0), 3, 3, byrow = T)

x <- c(.1, .2, .3)

toy_drift <- function(x, K = matrix(c(0, 0.1, 0.1,
                                      0.1, 0, 0.1,
                                      0.1, 0.1, 0), 3, 3, byrow = T),
                      b = c(1, 0.95, 1.05),
                      k = c(1, 1, 1),
                      hill_n = 2) {
  
  f <- sapply(1:length(x), function(i) {
    
    prod((K[i, -i]^hill_n)/(K[i, -i]^hill_n + x[-i]^hill_n))
    
  })
  
  x_new <- x*(b*f - k*x)
  
  return(x_new)
}

dispersion <- function(x) {
  sqrt(error)
}



## Euler - Maruyama ************************************** ####
x <- c(.1, .2, .3)

## Models
multiv_em <- function(X0, delta_t = 0.1, n, drift = "toy_drift", error = 0.01) {
  
  drift <- match.fun(drift)
  
  dim <- length(X0)
  
  time_points <- seq(from = 0, to = n, by = delta_t)
  
  X <- matrix(NA, length(time_points), dim)
  X[1, ] <- X0
  
  for(i in 2:length(time_points)){
    
    mu <- X[i-1, ] + drift(X[i-1, ])*delta_t
    sigma <- delta_t*error*diag(3)
    
    X[i, ] <- rmvnorm(1, mu, sigma)
      
  }
  
  X <- X %>% 
    cbind(time_points, .) %>% 
    as.data.frame() %>% 
    set_colnames(c("time", "X", "Y"))
  
  return(X)
  
}




Y <- multiv_em(X0 = c(.1, .2, .3), delta_t = 0.001, n=50, drift = "toy_drift", error = 0.001)



last_non_inf <- (sapply(1:nrow(Y), function(i) {
  abs(sum(Y[i, ])) > 10
}) %>% which())[1] - 1

Y[1:last_non_inf, ] %>% 
  melt(id.vars ="time") %>% 
  ggplot(aes(x = time, y= value, color = variable)) + 
  geom_line()




Y %>% head(n = 500) %>% 
  ggplot(aes(x = X, y = Y, color = time)) +
  geom_path()

