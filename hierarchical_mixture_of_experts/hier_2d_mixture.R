# Mixture of experts: two OUPs

## Settings ********* ####
# x values
n_obs <- 101
x_val <- seq(from = 0, to = 100, length.out = n_obs)

## Functions

# sigmoid function as the latent process
sigmoid_function <- function(x) {
  
  return(1/(1 + exp(-.5*x + 25)))
}

# Step function
stepper <- function(x) {
  
  ifelse(x > 33 & x < 66, 0, 1)
  
}

normie <- function(x) {
  1 - (dnorm(x, 50, 5)/dnorm(50, 50, 5))
}


line <- function(x) {
  
  x/100
  
}

func <- sigmoid_function


## Parameters ******* ####
mu <- list(c(0, 0), c(4, 2), c(2, 4))
B <- list(matrix(c(1, 0,
                   0, 1), ncol=2, nrow=2, byrow = T),
          matrix(c(.5, 0,
                   0, .5), ncol=2, nrow=2, byrow = T),
          matrix(c(.75, 0,
                   0, .75), ncol=2, nrow=2, byrow = T))

G <- list(matrix(c(.5, 0,
                   0, .5), ncol=2, nrow=2, byrow = T),
          matrix(c(.5, 0,
                   0, .5), ncol=2, nrow=2, byrow = T),
          matrix(c(.5, 0,
                   0, .5), ncol=2, nrow=2, byrow = T))



## Generate data **** ####
oup_2d_mixture_list() <- function(weight, Y, mu, B, G) {
  
  
  if(class(Y) == "numeric") {
    Y <- Y %>% as.matrix()
  }
  
  exp_mB1 <- expm(-B1) %>% as.matrix()
  Mean1 <- mu1 + exp_mB1%*%(Y - mu1)
  Cov1 <- G1 - exp_mB1%*%G1%*%exp_mB1 %>%
    as.matrix()
  
  
  exp_mB2 <- expm(-B2) %>% as.matrix()
  Mean2 <- mu2 + exp_mB2%*%(Y - mu2)
  Cov2 <- G2 - exp_mB2%*%G2%*%exp_mB2 %>%
    as.matrix()
  
  oup <- rbind(rmvnorm(n = 1, mean = Mean1, sigma = Cov1), 
               rmvnorm(n = 1, mean = Mean2, sigma = Cov2))
  
  choose <- sample(1:2, 1, prob = c(weight, 1 - weight))
  
  return(oup[choose, ])
  
}



seed <- sample(1:1000, 1)
set.seed(seed)
y <- matrix(NA, nrow = length(x_val), ncol = 2)
y[1, ] <- c(1, 1)

for(x in 2:length(x_val)) {
  print(x)
  delta_t <- x_val[x] - x_val[x-1]
  
  obs <- oup_2d_mixture(weight = func(x_val[x]),Y =  y[x-1, ], mu1, mu2, B1, B2, G1, G2)
  
  y[x, ] <- obs
}


df <- data.frame(x = y[, 1], y = y[, 2], time = x_val)


## Plot ************* ####
ggplot(df, aes(x = x, y = y, color = time)) + 
  geom_point() + 
  geom_path() +
  scale_color_gradient()



## Stan ************* ####
fit_2d_mixture_oup <- stan_model("2d_oup/2d_mixture_fitter.stan")


samples <- sampling(fit_2d_mixture_oup, 
                    list(N = 101, x = df$time, y = df[, 1:2]), 
                    iter = 1000, 
                    chains = 1)


summary(samples)$summary[, c("25%", "75%")] %>% 
  as.data.frame() %>% 
  set_colnames(c("lower", "upper")) %>% 
  rownames_to_column("parameter") %>% 
  filter(parameter != "lp__") %>% 
  filter(!grepl("mix", parameter)) %>% 
  filter(!grepl("_G", parameter)) %>% 
  mutate(true_value = c(mu1[1], mu1[2], mu2[1], mu2[2],
                        B1[1, 1], B1[1, 2], B1[2, 1], B1[2, 2],
                        B2[1, 1], B2[1, 2], B2[2, 1], B2[2, 2],
                        G1[1, 1], G1[1, 2], G1[2, 1], G1[2, 2],
                        G2[1, 1], G2[1, 2], G2[2, 1], G2[2, 2])) %>% 
  ggplot() +
  geom_point(aes(x = parameter, y = true_value), color = "red") +
  geom_errorbar(aes(x = parameter, ymin =lower, ymax =upper), width = 0.1)




summary(samples)$summary[, c("mean", "25%", "75%")] %>% 
  as.data.frame() %>% 
  set_colnames(c("mean", "lower", "upper")) %>% 
  rownames_to_column("parameter") %>% 
  filter(grepl("mix_", parameter)) %>% 
  mutate(true_value = func(x_val), x = x_val) %>% 
  ggplot() +
  geom_line(aes(x = x, y = true_value), color = "red") +
  geom_ribbon(aes(x = x, ymin =lower, ymax = upper), fill = "grey", alpha = .25) +
  geom_line(aes(x = x, y = mean))
