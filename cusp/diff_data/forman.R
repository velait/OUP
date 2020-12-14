## Simulate data from forman, use cusp to infere
library(KScorrect)

## Data ********************************************** ####
N <- 9


cusp_parameters <- c("alpha", "beta", "lambda", "epsilon", "r" )

oup_simulator <- stan_model("fiddling/stan_models/simulate_gp_OUP.stan")



rho <- 5
sigma <- 0.01
alpha <- 1
n_total_samples <- 101
x_total <- seq(from = 0, to = 100,length.out = n_total_samples)
# seed <- sample(1:100, 1)
 
 

mix_mean = list(c(1, 3), c(1, 3), c(1, 3), c(1, 3),
                c(0, 3), c(0, 3), c(0, 3), c(0, 3), c(0, 3))

mix_sd = list(c(.5, .5), c(.5, .75), c(.5, 1), c(.25, .75), 
              c(.5, .5), c(.5, .75), c(.5, 1), c(.25, .75), c(.25, .5))


mix_pro = list(c(.5, .5), c(.5, .5), c(.5, .5), c(.5, .5), c(.5, .5),
               c(.5, .5), c(.5, .5), c(.5, .5), c(.5, .5))



# marginal_distributions <- data.frame(x = -500:500/100) %>% 
#   mutate(oup = dnorm(x, 0, alpha),
#          bimodal = dmixnorm(x, mix_mean, mix_sd, mix_pro))
# 
# marginal_distributions %>% melt(id.vars = "x") %>% ggplot(aes(x = x, y=value, color = variable)) + geom_line() 

## Stan models ******************** ####

# Simulator
oup_simulator <- stan_model("fiddling/stan_models/simulate_gp_OUP.stan")


## Data *************************** ####
# Simulate data (OUP)
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
                        seed=sample(1:100, 1),
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
  mutate(bi = transformator(y, alpha, mix_mean, mix_sd, mix_pro)) %>% 
  set_colnames(c("x", "oup", "bimodal"))


forman <- lapply(1:N, function(i) {
    transformator(processes$y, alpha, mix_mean[[i]], mix_sd[[i]], mix_pro[[i]])
}) %>% do.call(cbind, .) %>% 
  as.data.frame() %>%
  set_colnames(paste0("y", 1:N))






## Inference ***************************************** ####

shoji_cusp_with_r <- stan_model("cusp/shoji_cusp_with_r.stan")
# shoji_cusp <- stan_model("cusp/shoji_cusp.stan")


forman_samples <- lapply(1:N, function(i) {
  
  # Stan data
  
  forman_data <- list(N_obs = n_total_samples,
                      N_OTUs = 1,
                      x = 1:n_total_samples, 
                      y = forman[1:101, 1] %>% as.matrix %>% t)
  
  
  samples <- sampling(shoji_cusp_with_r,
           forman_data,
           iter = 1000,
           chains = 1,
           control = list(adapt_delta = 0.8))
  
  
  samples
})
  




## Results ******************************************* ####
forman_results_df <- lapply(1:N, function(i) {
  
  lapply(1:length(cusp_parameters),
         function(x) get_stan_results(forman_samples[[i]],
                                      paste0("theta\\[.+,", x,"\\]"))) %>%
    do.call(rbind, .) %>% 
    cbind(parameter = cusp_parameters, 
          series = i)
  
}) %>% do.call(rbind, .)


## Plot estimates ************************************ ####


forman_results_df %>% 
  ggplot() +
  geom_errorbar(aes(x = as.factor(series), ymin = lower2.5, ymax = upper97.5)) + 
  facet_wrap(~parameter)



## Draw samples, simulate, stationary density ******** ####

n_sim <- 25

forman_cusp_simulated_densities <- lapply(1:N, function(i) {
  
  # number of samples
  n_samples <- forman_samples[[i]]@stan_args[[1]]$iter - forman_samples[[i]]@stan_args[[1]]$warmup
  
  # samples indices
  sample_ind <- sample(1:n_samples, n_sim)
  
  par_samples <- rstan::extract(forman_samples[[i]], "theta")[[1]][sample_ind, , ]
  
  
  
  # simulate series
  ticks <- 0.1
  series_sim <- lapply(1:n_sim, function(s) {
    
    shoji_generator(y0 = rnorm(1, par_samples[s, 1], 1), times = seq(from=0, to=n_total_samples, by = ticks),
                    # r = 1,
                    r = par_samples[s, 5],
                    alpha = par_samples[s, 1],
                    beta = par_samples[s, 2],
                    lambda = par_samples[s, 3],
                    epsilon = par_samples[s, 4],
                    seed = sample(1:1000, 1)) %>% 
      as.data.frame() %>%
      cbind(time = seq(from=0, to=n_total_samples, by = ticks)) %>%
      mutate(sample = s)
  }) %>% 
    do.call(rbind, .) %>% 
    cbind(series = i) %>% 
    set_colnames(c("value", "time", "sample", "series"))
  
  
  
  # thin: dt = 1
  series_sim <- thin(series_sim, time = "time")
  
  return(series_sim)
  
}) %>% do.call(rbind, .)
forman_map_densities <- lapply(1:N, function(i) {
  
  # simulate series
  ticks <- 0.1
  
  par_samples <- forman_results_df %>% filter(series == i)
    
  series_sim <- shoji_generator(y0 = rnorm(1, par_samples[3, "mode"], 1), times = seq(from=0, to=n_total_samples, by = ticks),
                    # r = 1,
                    r = par_samples[5, "mode"],
                    alpha = par_samples[1, "mode"],
                    beta = par_samples[2, "mode"],
                    lambda = par_samples[3, "mode"],
                    epsilon = par_samples[4, "mode"],
                    seed = sample(1:1000, 1)) %>% 
      as.data.frame() %>%
      cbind(time = seq(from=0, to=n_total_samples, by = ticks)) %>%
      mutate(sample = s) %>% 
    cbind(series = i) %>% 
    set_colnames(c("value", "time", "sample", "series"))
  
  
  
  # thin: dt = 1
  series_sim <- thin(series_sim, time = "time")
  
  return(series_sim)
  
}) %>% do.call(rbind, .)


 ggplot() +
  stat_density(data = forman_cusp_simulated_densities %>% mutate(series = paste0("y", series)),
               aes(x = value, group = sample),
               geom = "line",
               position = "identity",
               color = "red",
               alpha = 0.25) +
  stat_density(data =  cbind(forman[1:n_total_samples, ], 1:n_total_samples) %>%
                 set_colnames(c(paste0("y", 1:N), "time")) %>%
                 melt(id.vars = "time") %>% mutate(series = variable),
               aes(x = value),
               geom = "line", position = "identity", size = 1) +
  
   stat_density(data =forman_map_densities %>% mutate(series = paste0("y", series)), aes(x = value),
                geom = "line", position = "identity", color = "blue", size = 1) +
 
  facet_wrap(~series, scales = "free") +
  labs(subtitle = "OUP transformation and posterior samples; 100 time points") 



## PPC *********************************************** ####

ggplot() +
   geom_line(data = forman_cusp_simulated_densities %>% mutate(series = paste0("y", series)),
                aes(x = time, y = value, group = sample),
                position = "identity",
                color = "red",
                alpha = 0.25)  +
   geom_line(data =  cbind(forman[1:n_total_samples, ], 1:n_total_samples) %>%
                  set_colnames(c(paste0("y", 1:N), "time")) %>% 
               melt(id.vars = "time") %>% mutate(series = variable), 
                aes(x = time, y = value),
                position = "identity") +
   geom_line(data =  forman_map_densities %>% mutate(series = paste0("y", series)), 
             aes(x = time, y = value),
             position = "identity", color = "blue") +
   facet_wrap(~series)
 
 
 
 
 
 
 
 
 
 
 
 
 
 