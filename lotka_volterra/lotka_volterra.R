## Competitive Lotka Volterra model; n species


## Parameters ************************* ####

n_species <- 4
stochastic <- TRUE

A <- matrix(c(1, 1.09, 1.52, 0,
              0, 1, .44, 1.36, 
              2.33, 0, 1, .47,
              1.21, .51, .55, 1),
            4, 4, byrow = TRUE)

r <- c(1, 0.72, 1.53, 1.27)
K <- rep(1, ncol(A))


## Model
lotka <- function(times, y, parms) {
  
  dy <- r*y*(1 - (A %*% y)/K) %>% 
    as.vector
  
  list(dy)
}

# Error function, adds stochasticity
errorer <- function(x) {
  
  sqrt(abs(x))*rnorm(1, 0, 0.05)
  
}

## Initial values ********************* ####
times <- seq(from = 0, to = 100, by = 0.1)
x <- matrix(NA, length(times), n_species + 1)
x[1, ] <- c(0, runif(4, 0, 1))


## Model ****************************** ####
for(i in 2:length(times)) {
  
  out <- ode(times = times[(i-1):i],
             y = x[i-1, -1],
             func = lotka,
             parms = NULL,
             method = "euler")
  
  x_new <- out[2, 2:ncol(out)] 
  
  # x_new <- x_new + rnorm(3, 0, 0.001)
  
  # x_new <- x_new
  
  if(stochastic) {
    x_new <- x_new + sapply(x_new + 0.01, errorer)
  }
  
  x_new <- ifelse(x_new < 0, 0, x_new)
  
  x[i, ] <- c(out[2, 1], x_new)
  
}


## Plot ******************************* ####
x %>% 
  as.data.frame() %>% 
  set_colnames(c("time", paste("V", 1:(n_species)))) %>% 
  melt(id.vars = "time") %>%
  ggplot(aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs("Competetive lotka volterra")


