# cups catastrophe model
# see Iacus p.57

# seed <- 1
seed = sample(1:1000, 1)

# 
r <- 1
alpha <- -1
beta <- 1.5
lambda <- 0
epsilon <- 1




C <- 27*alpha^2 - 4*beta^3



x <- seq(-5, 5, length.out = 100)

inv_density_plot <- data.frame(x = x, y = cc_density(x, r, alpha, beta, lambda, epsilon)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line() +
  labs(title = "Population level density",
       x = "Abundance", y = "Density")

inv_density_plot




## Shoji *********************************** ####
grid <- seq(from = 0, to = 100, by = 0.01)

inits <- c(0)

df <- lapply(inits, function(y0) shoji_generator(y0 = y0, times = grid,
                                                 seed = seed,
                                                 r, alpha, beta, lambda, epsilon)) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  cbind(time = grid) %>%
  set_colnames(c("y", "x"))
  # set_colnames(c("y1", "y2", "x"))

thin_df <- (df)

diff_plot <- thin_df %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value, color = variable)) +
  geom_line() +
  # scale_color_tron() +
  labs(title = "", x = "Time", y = "Abundance") +
  theme(legend.position = "none") 

diff_plot



plot_grid(inv_density_plot, diff_plot, rel_widths = c(1, 1.5)) 



## Cusp fitter 1 *************************** ####

Euler_Maruyama <- stan_model(file = "euler_maruyama/euler_maruyama.stan")
Euler_Maruyama_samples <- sampling(Euler_Maruyama,
                                   data = list(N = length(thin_df$x),
                                               x = thin_df$x,
                                               y = thin_df$y),
                                   iter = iter,
                                   chains = chains,
                                   init = 1)

euler_parameters <- c("alpha", "beta", "lambda", "r", "epsilon")

euler_res <- summary(Euler_Maruyama_samples)$summary[grep(paste0(euler_parameters, collapse = "|"), rownames(summary(Euler_Maruyama_samples)$summary)), c("50%", "25%", "75%")] %>% 
  as.data.frame() %>% 
  set_colnames(c("mode", "lower25", "upper75"))


euler_res["r", 1]
euler_res["alpha", 1]
euler_res["beta", 1]
euler_res["lambda", 1]
euler_res["epsilon", 1]


x_grid <- seq(from = -5, to = 5, length.out = 250)
data.frame(x = x_grid,
           y_est = cc_drift(x_grid,
                            euler_res["r", 1],
                            euler_res["alpha", 1],
                            euler_res["beta", 1],
                            euler_res["lambda", 1],
                            euler_res["epsilon", 1]),
           y_sim = cc_drift(x_grid,
                            r, 
                            alpha, 
                            beta, 
                            lambda, 
                            epsilon)) %>% 
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value, color = variable)) + 
  geom_line()

euler_plot <- ggplot() +
  geom_line(data = df, aes(x = x, y = y)) +
  geom_hline(yintercept = c(euler_res[1:3, "mode"]),
             linetype = "dashed")

euler_plot

## Cusp fitter 2 *************************** ####
invariant_fitter <- stan_model(file = "cusp/cusp_fitter.stan")


invariant_samples <- sampling(invariant_fitter,
                                   data = list(N = length(thin_df$y),
                                               y = thin_df$y),
                                   iter = 500,
                                   chains = 1)

invariant_samples






## Randomly generate from cusp density ***** ####

sample_inv_density <- function(n, a, b) {
  
  max_y <- max(cc_density(seq(a, b, length.out = 1000), r, alpha, beta, lambda, epsilon))
  x <- runif(n, a, b)
  y <- runif(n, 0, max_y)
  ok <- cc_density(x, r, alpha, beta, lambda, epsilon) >= y
  
  return(x[ok])
}


xx <- rnorm(1000, 0, 1)

invariant_samples <- sampling(invariant_fitter,
                              data = list(N = length(xx),
                                          y = xx),
                              iter = iter,
                              chains = 2)

invariant_samples
## Fit with R function ********************* ####


r <- 1
alpha <- 0
beta <- 1
lambda <- 0
epsilon <- 1

grid <- seq(from = 0, to = 250, by = 0.01)
inits <- c(0)



rep_cusp <- replicate(2, 
{df <- lapply(inits, function(y0) em_generator(y0 = y0, grid = grid, seed = NULL, restrict_pos = FALSE)) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  cbind(time = grid) %>%
  set_colnames(c(paste0("y", 1:length(inits)), "x"))

thin_df <- df[(1:(nrow(df)) %% 100) == 1, 1]
thin_df}
)


ab_coefs <- apply(rep_cusp, 2, FUN = function(i) {
  
  fit <- cusp(y ~ i, alpha ~ 1, beta ~ 1)
  coefs <- fit$coefficients[1:2] %>% set_names(c("a", "b"))
  coefs
  
}) %>% cbind() %>% t %>% as.data.frame()


ab_coefs %>% 
  melt() %>% 
  ggplot(aes(x = variable, y = value)) + 
  geom_boxplot() + 
  geom_hline(yintercept = c(0, 2), linetype = "dashed")
  
## Pipe ************************************ ####

### 1. Simulate many cusp series

r <- 1
alpha <- 1
beta <- 0
lambda <- 0
epsilon <- 1


# Use reasonable sample sizes: 100 observations per 100 subjects
grid <- seq(from = 0, to = 100, by = 0.01)

rep_cusp <- replicate(100, 
                      {df <- lapply(runif(1, -1, 1), function(y0) em_generator(y0 = y0, grid = grid, seed = NULL, restrict_pos = FALSE)) %>%
                        do.call(cbind, .) %>%
                        as.data.frame() %>%
                        cbind(time = grid) %>%
                        set_colnames(c(paste0("y", 1:length(inits)), "x"))
                      
                      thin_df <- df[(1:(nrow(df)) %% 100) == 1, 1]
                      thin_df}
)


### 2. Fit the population density to these --> a, b
rep_cusp_vector <- rep_cusp %>% as.vector()
ab_coefs <- cusp(y ~ rep_cusp_vector, alpha ~ 1, beta ~ 1)
empirical_coefs <- ab_coefs$coefficients[1:2] %>% set_names(c("a", "b"))

boot_coef <- replicate(100, expr = {
  
  boot_vec <- sample(rep_cusp_vector, length(rep_cusp_vector), replace = T)
  
  coefs <- cusp(y ~ boot_vec, alpha ~ 1, beta ~ 1)
  coefs$coefficients[1:2] %>% set_names(c("a", "b"))
  
}) %>% 
  t %>% 
  as.data.frame()
  

boot_coef %>% 
  ggplot() +
  geom_point(aes(x = a, y = b)) + 
  geom_point(aes(x = empirical_coefs["a"], y = empirical_coefs["b"]), color = "red", size = 3) + 
  geom_point(aes(x = alpha, y = beta), color = "blue", size = 3)


# Use bootstrap means and variances as priors
empirical_coefs
apply(boot_coef, 2, FUN = sd)

### 3. Fit Euler-Maruyama fit priors based on 2. 
Euler_Maruyama <- stan_model(file = "euler_maruyama/euler_maruyama.stan")

EM_samples <- apply(rep_cusp, 2, FUN = function(i) {
  ss <- sampling(Euler_Maruyama,
           data = list(N = length(i),
                       x = 1:length(i),
                       y = i),
           iter = 1000,
           chains = 1)
  
  summary(ss)$summary[c("alpha", "beta", "lambda", "epsilon"), "mean"]
  
})
  

### 4. Results

EM_samples %>%
  t %>%
  cbind(., n = 1:ncol(rep_cusp)) %>%
  as.data.frame() %>% 
  melt(id.vars = "n") %>% 
  ggplot(aes(x = variable, y = value)) + 
  geom_boxplot()








