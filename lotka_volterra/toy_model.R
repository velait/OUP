library(deSolve)

theme_set(theme_bw(10))

## Parameters ************************* ####
n_species <- 5
stochastic <- TRUE

## Toy model
# k <- sample(c(1, 0.95, 1.05), n_species, replace = T)
k <- rep(1, n_species)
b <- sample(c(1, 0.95, 1.05), n_species, replace = T)
# b <- c(1, .95, 1.05)
# K <- matrix(runif(n_species, 0, 1) %>% round(2),
#             n_species, n_species, byrow = T)
K <- rep(0.1, n_species^2) %>% matrix(n_species, n_species, byrow = T)
K <- K - K*diag(n_species)




## Functions ************************** ####

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

## Toy model from "Multistability and the.... "; n species
toy_model <- function(times, y, parms) {
  
  f <- f_hill(y = y, K)
  
  dy <- y*(b*f - k*y)
  
  list(dy)
}



## Initial values ********************* ####

times <- seq(from = 0, to = 100, by = 0.1)
x <- matrix(NA, length(times), n_species + 1)
x[1, ] <- c(0, sample(1:10/10, n_species, replace = T))

## Model ****************************** ####


for(i in 2:length(times)) {
  
  out <- ode(times = times[(i-1):i],
             y = x[i-1, -1],
             func = toy_model,
             parms = NULL)
  
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

x <- x %>% as.data.frame() %>% set_colnames(c("time", paste0("V", 1:(n_species))))

x %>% 
  melt(id.vars = "time") %>%
  ggplot(aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(title = "Toy model") +
  theme(legend.position="none")



x[, 6] %>% hist


## DUMP ******************************* ####
# ODE; 3 species
# rigidode <- function(times, y, parms, p = 2) {
#   
#   f <- f_step(y, K, p = p)
#   
#   dy1 <- y[1]*(b[1]*f[1] - k[1]*y[1])
#   
#   dy2 <- y[2]*(b[2]*f[2] - k[2]*y[2])
#   
#   dy3 <- y[3]*(b[3]*f[3] - k[3]*y[3])
#   
#   dy <- c(dy1, dy2, dy3)
#   
#   list(dy)
# }




## The RigidODE problem



