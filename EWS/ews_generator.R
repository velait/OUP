# Simulate early warning data

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


r <- 1
K <- 10
cs <- 1
h <- 1
sigma <- 0.07

times <- seq(from = 0, to = 250, by = 0.1)


y <- rep(NA, length(times))
y[1] <- 8

# inflow <- rep(NA, length(times))
# inflow[1] <- 0.07*rnorm(1, 0, 1)


cs <- 1 + (2.7 - 1)*(times/max(times))


for(i in 2:length(times)) {
  
  
  c <- cs[i]
  # c <- cs[.575*length(cs)]
  # c <- 1 + (2.6771 - 1)*.75
  # c <- 2.6
    
  dt <- times[(i-1):i]
  
  y[i] <- ews_generator(y[i-1], dt, milstein = T)[2]
  
  
  # Add inflow
  # inflow[i] <- ((1 - 1/20)*inflow[i-1] + 0.07*rnorm(1, 0, 1))*y[i-1]
  
  # y[i] <- y[i] + 0.07*rnorm(1, 0, 1)*y[i-1]
  
}

# Measurement error
# y <- y + 0.1*rnorm(length(y), 0, 1)


df <- data.frame(y = log(1 + y), x = times)

df %>% 
  # mutate(y = (y-mean(y))/sd(y)) %>%
  # thin(modulo = 100) %>%
  ggplot(aes(x = x, y = y)) + 
  geom_vline(xintercept = df[which((cs < 2.604)) %>% max(), "x"]) + 
  geom_line()


df_data <- df[which((cs < 2.604)), ]

# log(y + 1) %>% density %>% plot
