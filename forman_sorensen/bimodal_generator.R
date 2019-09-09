# bimodality ala Forman
library(KScorrect)


rho <- 5
sigma <- 0.01
alpha <- 1
n_total_samples <- 1001
x_total <- seq(from = 0, to = 1000,length.out = n_total_samples)
# seed <- sample(1:100, 1)
# 


mix_mean = c(1, 3)
mix_sd = c(.5, .75)
mix_pro = c(.25, .75)

transformator <- function(x, alpha, mix_mean, mix_pro) {
  
  qmixnorm(pnorm(x, 0, alpha),
           mean = mix_mean, 
           sd = mix_sd,
           pro = mix_pro)
  
}

marginal_distributions <- data.frame(x = -500:500/100) %>% 
  mutate(oup = dnorm(x, 0, alpha),
         bimodal = dmixnorm(x, mix_mean, mix_sd, mix_pro))

marginal_distributions %>% melt(id.vars = "x") %>% ggplot(aes(x = x, y=value, color = variable)) + geom_line() 

## Stan models ******************** ####

# Simulator
oup_simulator <- stan_model("fiddling/stan_models/simulate_gp_OUP.stan")


## Data *************************** ####
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
    
  }) 
  
  
  res
  
})



processes <- data.frame(x = x_total,
                        y = summary(single_data_set[[1]][[1]])$summary[ grepl("f", rownames(summary(single_data_set[[1]][[1]])$summary)), "mean"]
) %>% 
  `rownames<-`(NULL)



processes <- processes %>% 
  mutate(bi = transformator(y, alpha, mix_mean, mix_pro)) %>% 
  set_colnames(c("x", "oup", "bimodal"))


# plot process
processes[1:250, ] %>% 
  ggplot() +
  geom_line(aes(x = x, y = oup),color = "red", alpha = .5) +
  geom_line(aes(x = x, y = bimodal)) +
  labs(subtitle = "Red = OUP; Black = transformed diffusion",
       y = "Value", x = "Time") +
  theme_bw(15)

# 
# 
# # Plot distributions
# processes %>%
#   melt(id.var = "x") %>%
#   ggplot(aes(x = value, fill = variable)) +
#   geom_histogram(alpha = .5, color = "black", position = "identity") +
#   scale_fill_tron() +
#   geom_line(data = marginal_distributions %>% melt(id.var = "x"), aes(x = x, y = value*400, color = variable)) +
#   scale_y_continuous(sec.axis = sec_axis(~./400, name = "Density"))
# 

