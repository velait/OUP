# bimodality ala Forman


rho <- 5
alpha <- 1
sigma <- 1
n_observations <- 100
n_total_samples <- 1001
x_total <- 1000 * (0:(n_total_samples - 1)) / (n_total_samples - 1)
seed <- sample(1:100, 1)

# Simulate data (including hold out data)
single_data_set <- lapply(rho, function(r) {
  
  res <- lapply(alpha, function(a) {
    
    gp_simulator_data <- list(N = n_total_samples,
                              x = x_total,
                              alpha = a,
                              rho = r,
                              sigma =  sigma)
    
    samples <- sampling(oup_simulator,
                        data = gp_simulator_data, 
                        iter=1,
                        chains=1,
                        seed=seed,
                        algorithm="Fixed_param")
    
    samples
    
  }) %>% set_names(alpha_grid %>% as.character())
  
  
  res
  
}) %>% set_names(rho_grid %>% as.character())



processes <- data.frame(x =x_total, y = summary(single_data_set[[1]][[1]])$summary[ grepl("f", rownames(summary(single_data_set[[1]][[1]])$summary)), "mean"]
) %>% 
  `rownames<-`(NULL)


x_grid <- -500:500/100
distributions <- data.frame(x = x_grid, y = dnorm(x_grid, 0, 1), bi = bimodal_distribution(x_grid))



transformator <- function(x) {
  
  qmixnorm(pnorm(x, 0, alpha),mean = c(3, 5), sd = c(0.1, 1), pro = c(.9, .1))
  
}


processes <- processes %>% 
  mutate(bi = transformator(y))


processes %>% 
  ggplot() +
  geom_line(aes(x = x, y = y),color = "red", alpha = .5) +
  geom_line(aes(x = x, y = bi))




# bimodal_quantile <- function(x) {
#   
#   # sapply(1:length(x), function(i) {
#   #   
#   #   return(sum(x[1:i]))
#   #   
#   # })
#   
#   
#   
# }
# 
# 
# 
# bimodal_distribution <- function(x) {
#   
#   y <- dnorm(x, -1, .5) + dnorm(x, 1, .5)
#   
#   return(y/sum(y))
#   
# }
