# Shoji-Ozaki for the cusp catastrophy model


r <- 1
alpha <- -1
beta <- 1
lambda <- 2
epsilon <- 1

shoji_cusp_transition_density <- function(x, dt, r, alpha, beta, lambda, epsilon) {
  drift <- r*(alpha + beta*(x - lambda) - (x - lambda)^3)
  
  L <- r*(beta - 3*(x - lambda)^2)
  M <- r*((epsilon/2)*(- 6*(x - lambda)))
  
  Ax <-  x + drift*(exp(L*dt) - 1)/L + M*(exp(L*dt) - 1 - L*dt)/L^2
  B <- sqrt(epsilon)*sqrt((exp(2*L*dt) - 1)/(2*L))
  
  rnorm(1, Ax, B)
  
}

shoji_generator <- function(y0, times, r, alpha, beta, lambda, epsilon, seed = 1) {
  set.seed(seed)
  y <- y0
  
  for(i in 2:length(times)) {
    dt <- times[i] - times[i-1]
    
    y <- c(y, shoji_cusp_transition_density(y[i-1], dt = dt
                                            , r = r, alpha = alpha, beta = beta, lambda = lambda, epsilon = epsilon))
  }
  
  y
}

times <- seq(from = 0, to = 10, by = 0.01)

shoji_df <- shoji_generator(y0 = 5, times = times,
                            r, alpha, beta, lambda, epsilon, seed = 123)

em_df <- em_generator(y0 = 5, grid = times,
                      seed = 123, restrict_pos = FALSE)



cbind(shoji_df, em_df, x = times) %>%
  as.data.frame() %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value, color = variable)) +
  geom_line() +
  labs(title = "Euler vs. Shoji; dt = 0.15")
  

