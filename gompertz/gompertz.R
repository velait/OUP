# Cusp plus gompertz



## Data ******************************** ####
N <- 9

cusp_parameters <- c("alpha", "beta", "lambda", "epsilon", "r" )

r <- runif(N, .5, 1.5)
alpha <- c(0, 0, 0, 0, 0,
           2, -0.5, -1, 1)
# beta <- c(2, -2, 1.5, -3, -1,
#           1, -1,  0, 0)
beta <- rep(2, N)
# alpha <- rnorm(N, 0, .5)
# beta <- rnorm(N, 0, 2)
lambda <- runif(N, 3, 5)
epsilon <- runif(N, 1, 2)


n_points <- 100

grid <- seq(from = 0, to = n_points, by = 0.1)
seed = sample(1:n_points, 1)


inits <- rnorm(N, lambda[1], 0.1)

df_shoji <- lapply(1:N, function(i) shoji_generator(y0 = inits[i], times = grid,
                                                    r = r, alpha = alpha[i], beta = beta[i],
                                                    lambda = lambda[i], epsilon = epsilon[i],
                                                    seed = sample(1:1000, 1))) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  cbind(time = grid) %>%
  set_colnames(c(paste0("y", 1:N), "x"))


df_shoji <- thin(df_shoji)

df_shoji %>% 
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value)) +
  geom_line() +
    facet_wrap(~variable)


## Exponentiate

gomp_cusp <- df_shoji %>% 
  select(-x) %>% 
  apply(MARGIN = 2, FUN = function(i) exp(i)) %>% 
  cbind(x = df_shoji$x)

gomp_cusp_pois <- lapply(1:N, function(i) {
  
  lapply(1:nrow(gomp_cusp), function(t) {
    
    rpois(1, gomp_cusp[t, i])
    
  }) %>% unlist
  
}) %>% do.call(cbind, .) %>% as.data.frame() %>% cbind(x = gomp_cusp[, "x"])

gomp_cusp_pois %>% 
  as.data.frame() %>% 
  melt(id.vars = "x") %>% 
  ggplot() +
  geom_line(aes(x = x, y = value)) + 
  facet_wrap(~variable, scales = "free")







## STAN ******************************** ####

stan_gompertz <- stan_model("gompertz/gompertz_cusp.stan")

gompertz_samples <- lapply(1:N, function(i) {
  
  gomp_data <- list(N_obs = nrow(gomp_cusp_pois),
                    N_OTUs = 1,
                    x = gomp_cusp_pois$x,
                    y = (select(gomp_cusp_pois, -x))[, i] %>% as.matrix() %>% t)
  
  
   sampling(stan_gompertz, 
            gomp_data,
            iter = 1000, 
            chains = 1, 
            control = list(
            adapt_delta = 0.8, 
            max_treedepth = 10
            )
  )
  
  
})




results_with_r_df <- lapply(1:5, function(x) get_stan_results(gompertz_samples, paste0("theta\\[.+,", x,"\\]"))) %>% do.call(rbind, .) %>% 
  cbind(parameter = rep(c("alpha", "beta", "lambda", "epsilon", "r"), each = N),
        index = factor(rep(1:N, 5)),
        true = c(alpha, beta, lambda, epsilon, r))



results_with_r_df %>% 
  ggplot() + 
  geom_point(aes(x = index, y = true), shape = 4) + 
  geom_errorbar(aes(x = index, ymax = upper97.5, ymin = lower2.5), width = .2) +
  facet_wrap(~parameter, scales = "free")



## get latent states

latent_df <-  lapply(1:N, function(x) {
  df <- get_stan_results(gompertz_samples, paste0("latent_cusp\\[", x,",.+\\]"))
  
  df <- df %>%
    as.data.frame() %>% 
    mutate(time = 1:nrow(df), series = x, true = df_shoji[, x]) 
  
  df
  
} ) %>% do.call(rbind, .)


latent_df %>% 
  filter(time != 1) %>% 
  ggplot() + 
  geom_line(aes(x = time, y = true)) + 
  geom_line(aes(x = time, y = mode), color = "red") + 
  geom_ribbon(aes(x = time, ymin = lower2.5, ymax = upper97.5), alpha = .25) +
  facet_wrap(~series, scales = "free")


## Simulate data

# map_simulated_densities_with_r <- lapply(1, function(i) {
#   
#   
#   # simulate series
#   series_sim <- lapply(1, function(s) {
#     lapply(1:N, function(j) {
#       
#       map_pars <- results_with_r_df %>%
#         filter(index == j) %>%
#         select(mode, parameter)
#       
#       shoji_generator(y0 = rnorm(1, 0, 1), times = seq(from=0, to=100, by = 0.1),
#                       # r = 1,
#                       map_pars[5, 1],
#                       alpha = map_pars[1, 1],
#                       beta = map_pars[2, 1],
#                       lambda = map_pars[3, 1],
#                       epsilon = map_pars[4, 1],
#                       seed = sample(1:1000, 1))
#     }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame() %>%
#       cbind(time = seq(from=0, to=100, by = 0.1)) %>%
#       set_colnames(c(paste0("y", 1:N), "x")) %>% 
#       mutate(sample = s)
#   }) %>% 
#     do.call(rbind, .)
#   
#   
#   
#   # thin: dt = 1
#   series_sim <- thin(series_sim)
#   series_sim_no_time <- series_sim %>% select(-c(x, sample))
#   
#   series_sim_no_time <- series_sim_no_time %>% 
#     apply(MARGIN = 2, FUN = function(i) {
#       cols <- exp(i)
#       
#       
#       rpois(length(cols), cols)
#       
#       })
#   
#   
# })






gomp_oup_pois <- rnorm(100, exp(oup_generator(mu = 5, y0 = 5)), 5)
gomp_oup_pois <- cbind(gomp_oup_pois, x= 1:100) %>% as.data.frame() %>% set_colnames(c("y1", "x"))


oup_gompertz <- stan_model("gompertz/gompertz_oup.stan")


oup_samples <- lapply(1, function(i) {
  
  gomp_data <- list(N_obs = nrow(gomp_oup_pois),
                    N_OTUs = 1,
                    x = gomp_oup_pois$x,
                    y = (select(gomp_oup_pois, -x))[, i] %>% as.matrix() %>% t)
  
  
  
  
  
  
  
  sampling(oup_gompertz, 
           gomp_data,
           iter = 1000, 
           chains = 1, 
           control = list(
             adapt_delta = 0.8, 
             max_treedepth = 10
           )
  )
  
  
})


## OUP gompertz ************************ ####
oup_gompertz <- oup_generator(N = 100, times = 1:100, y0 = 5, mu = 5, lambda = 0.5, sigma = .25)


oup_gompertz <- exp(oup_gompertz)

oup_gompertz <- lapply(1:100, function(i) {
  
  rpois(1, oup_gompertz[i])
  
}) %>% unlist()


oup_gompertz_stan_data <- list(N_obs = 100, N_OTUs = 1, x = 1:100, y = oup_gompertz %>% as.array %>% t)


oup_gompertz_stan_model <- stan_model("gompertz/gompertz_oup.stan")



oup_gompertz_samples <- sampling(oup_gompertz_stan_model, oup_gompertz_stan_data, 
                                 chains = 1, iter = 1000)



# Compare



oup_gompertz_generator <- function(N = 100, times = 1:100, y0 = 5, mu = 5, lambda = 0.5, sigma = .25) {
  
  x <- oup_generator(N = N, times = times, y0 = y0, mu = mu, lambda = lambda, sigma = sigma)
  
  
  x <- exp(x)
  
  x <- lapply(1:100, function(i) {
    
    rpois(1, x[i])
    
  }) %>% unlist()
  
  x
}






data.frame(res = oup_gompertz_generator(100, 1:100, 5, 4.77, .5, .81), 
      data = oup_gompertz, 
      x = 1:100) %>% 
  melt(id.vars = "x") %>% 
  ggplot(aes(x = x, y = value, color = variable)) + 
  geom_line()
  
