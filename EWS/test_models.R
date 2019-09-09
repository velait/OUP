## Data *********************


#### Regular change point inference
# Get data
source("EWS/ews_generator.R")

ews_data <- list(Y = df %>% thin(modulo = 1) %>% pull(y),
                 time = df %>% thin(modulo = 1) %>% pull(x), 
                 T = nrow(df %>% thin(modulo = 1)))

## Basic OUP **************** ####

oup_parameters <- c("lambda", "sigma", "mu")
oup_model <- stan_model("stan_models/oup_single.stan")
oup_samples <- sampling(oup_model, ews_data, 
                        chains = 1,
                        iter = 1000)


# Sliding window inference
window_length <- 50
N_windows <- length(ews_data[["Y"]]) - window_length

sliding_oup_res <- lapply(1:N_windows, function(i) {
  
  data <- list()
  
  data[["Y"]] <- ews_data[["Y"]][i:(window_length -1 + i)]
  data[["time"]] <- ews_data[["time"]][i:(window_length -1 + i)]
  data[["T"]] <- window_length
  
  
  oup_samples <- sampling(oup_model, data, 
                          chains = 1,
                          iter = 500)
  
  df <- summary(oup_samples)$summary[oup_parameters, c("mean", "2.5%", "25%", "50%", "75%", "97.5%")] %>% 
    set_colnames(c("mean", "lower2.5", "lower25", "mode", "upper75", "upper97.5")) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "parameter") %>% 
    mutate(window = i, from = i, to = (window_length -1 + i)) 
    
  return(df)
  
}) %>% do.call(rbind, .)


## Results ****************

ggplot() + 
  geom_line(data = df, aes(x = df$x, y = df$y)) + 
  geom_vline(xintercept = df[which((cs < 2.604)) %>% max(), "x"], linetype = "dashed") + 
  geom_line(data = sliding_oup_res, aes(x = to, y = mode, color = parameter)) + 
  geom_line(aes(x = 50:250, y = .4*(sliding_oup_res %>% filter(parameter == "lambda") %>% pull(mode))/(sliding_oup_res %>% filter(parameter == "sigma") %>% pull(mode))), color = "darkgoldenrod") +
  geom_ribbon(data = sliding_oup_res, aes(x = to, ymin = lower2.5, ymax = upper97.5, fill = parameter), alpha = 0.25) + 
  labs(x = "x", y = "y", subtitle = "OUP fit to sliding windows,  dl = 50")







#### Flickering data set
source("EWS/flickerer.R")

flicker_data <- list(Y = flicker_df %>% pull(y),
                 time =flicker_df %>% pull(x), 
                 T = nrow(flicker_df))


N_windows <- length(flicker_data[["Y"]]) - 50

flicker_sliding_oup_res <- lapply(1:N_windows, function(i) {
  
  data <- list()
  
  data[["Y"]] <- flicker_data[["Y"]][i:(window_length -1 + i)]
  data[["time"]] <- flicker_data[["time"]][i:(window_length -1 + i)]
  data[["T"]] <- window_length
  
  
  oup_samples <- sampling(oup_model, data, 
                          chains = 1,
                          iter = 500)
  
  df <- summary(oup_samples)$summary[oup_parameters, c("mean", "2.5%", "25%", "50%", "75%", "97.5%")] %>% 
    set_colnames(c("mean", "lower2.5", "lower25", "mode", "upper75", "upper97.5")) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "parameter") %>% 
    mutate(window = i, from = i, to = (window_length -1 + i)) 
  
  return(df)
  
}) %>% do.call(rbind, .)


ggplot() + 
  geom_line(data = flicker_df, aes(x = flicker_df$x, y = flicker_df$y)) + 
  geom_line(data = flicker_sliding_oup_res, aes(x = to, y = mode, color = parameter)) + 
  geom_line(aes(x = 51:1000, y = .4*(flicker_sliding_oup_res %>% filter(parameter == "lambda") %>% pull(mode))/(flicker_sliding_oup_res %>% filter(parameter == "sigma") %>% pull(mode))), color = "darkgoldenrod") +
  geom_ribbon(data = flicker_sliding_oup_res, aes(x = to, ymin = lower2.5, ymax = upper97.5, fill = parameter), alpha = 0.25) + 
  labs(x = "x", y = "y", subtitle = "OUP fit to sliding windows,  dl = 50")






## HMM ********************** ####
 


## Mixture of experts ******* ####
expert_mixture <- stan_model("mixture_of_experts/fit_oup_mixture.stan")

# Test on flikcering data set

data <- list()

data[["y"]] <- flicker_data[["Y"]]
data[["x"]] <- flicker_data[["time"]]
data[["N"]] <- flicker_data[["Y"]] %>% length


expert_mixture_samples <- sampling(expert_mixture, data, 
                        chains = 1,
                        iter = 1000)

df <- summary(expert_mixture_samples)$summary[1:(nrow(summary(expert_mixture_samples)$summary) - 1), c("mean", "2.5%", "25%", "50%", "75%", "97.5%")] %>% 
  set_colnames(c("mean", "lower2.5", "lower25", "mode", "upper75", "upper97.5")) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "parameter")
  
  
  
df <- df %>% 
  separate(col = parameter, sep = "\\[", into = c("parameter", "index")) %>% 
  mutate(index = gsub("\\]", "", index) %>% as.numeric)


ggplot() + 
  geom_line(data = flicker_df, aes(x = flicker_df$x, y = flicker_df$y)) + 
  geom_line(data = df %>% filter(parameter == "mix_weight"), aes(x = index, y = mode), color = "red") +
  geom_ribbon(data = df %>% filter(parameter == "mix_weight"), aes(x = index, ymin = lower2.5, ymax = upper97.5, fill = parameter), alpha = 0.25) + 
  labs(x = "x", y = "y", subtitle = "")


## Cusp model *************** ####

cusp_model <- stan_model("cusp/shoji_cusp_with_r.stan")


cusp_parameters <- c("alpha", "beta", "lambda", "epsilon", "r")
window_length <- 50
N_windows <- length(ews_data[["Y"]]) - window_length

sliding_cusp_res <- lapply(1:N_windows, function(i) {
  
  data <- list()
  
  data[["y"]] <- ews_data[["Y"]][i:(window_length -1 + i)] %>% as.matrix() %>% t
  data[["x"]] <- ews_data[["time"]][i:(window_length -1 + i)]
  data[["N_obs"]] <- window_length
  data[["N_OTUs"]] <- 1
  
  
  cusp_samples <- sampling(cusp_model,
                          data, 
                          chains = 1,
                          iter = 500)
  
  df <- get_stan_results(cusp_samples, cusp_parameters, regex = T) %>% 
    rownames_to_column(var = "parameter") %>% 
    mutate(parameter = parameter %>% gsub("\\[[0-9]\\]", "", .)) %>% 
    mutate(window = i, from = i, to = (window_length -1 + i)) 
  
  return(df)
  
}) %>% do.call(rbind, .)


## Results ****************

ggplot() + 
  geom_line(aes(x = ews_data$time, y = ews_data$Y)) + 
  geom_line(data = sliding_cusp_res, aes(x = to, y = mode, color = parameter)) + 
  # geom_ribbon(data = sliding_cusp_res, aes(x = to, ymin = lower2.5, ymax = upper97.5, fill = parameter), alpha = 0.25) + 
  labs(x = "x", y = "y", subtitle = "Cusp fit to sliding windows,  dl = 50") + 
  coord_cartesian(ylim = c(-1, 3))
  

## EWS package ************** ####

library(earlywarnings)

ews <- earlywarnings::generic_ews(ews_data$Y, winsize = 20)







